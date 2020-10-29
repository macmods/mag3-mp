function envt = envt_les_v9(envt,farm,dir_LES)
% Load velocity deficits from LES files
% INPUT: envt, farm, and les directory where fields are saved
% OUTPUT: 3-d velocity deficit fields + diffusivity for three kelp scenarios
%   u_def_a; u_def_b; u_def_c
%   v_def_a; v_def_b; v_def_c
%   w_def_a; w_def_b; w_def_c
%       a = "subsurface" and deep mixed layer
%       b = "subsurface" and shallow mixed layer
%       c = "max canopy" and deep mixed layer
%       d = "max canopy" and shallow mixed layer
%       These are normalized velocity deficits and will be scaled with flow
%       from ROMS conditions
%   K_a; K_b;
%       normalized vertical diffusivity
%       K = K_a * u-star * cultivation depth
%       u-star = 0.0061 cm/s in LES case and will be scaled with mean flow
%       of ROMS conditions

%% File Directory
    lesfiles = dir(fullfile(dir_LES,'*.csv'));
    
    % There are four files
    % 1. full canopy; deep mixed layer
    % 2. full canopy; shallow mixed layer
    % 3. subsurface; deep mixed layer
    % 4. subsurface; shallow mixed layer
    
    for ff = 1:length(lesfiles)
        filename = strcat(lesfiles(ff).folder,'/',lesfiles(ff).name);
        les = readtable(filename);
        
    %% Define LES grid
    % Information from LES group; note that dx/dy/dz not determined dynamically

    lesgrid.x = max(unique(les.x)); lesgrid.dx = 2;
    lesgrid.y = max(unique(les.y)); lesgrid.dy = 2;
    lesgrid.z = max(abs(unique(les.z))); lesgrid.dz = 0.5;
    
    % Translate single column grid to matrix
    % Order of shared LES files are along x,then y, then z
    lesgrid.u = reshape(les.u,length(unique(les.x)),length(unique(les.y)),length(unique(les.z)));
    lesgrid.v = reshape(les.v,length(unique(les.x)),length(unique(les.y)),length(unique(les.z)));
    lesgrid.w = reshape(les.w,length(unique(les.x)),length(unique(les.y)),length(unique(les.z)));
    lesgrid.K = reshape(les.Kdiff,length(unique(les.x)),length(unique(les.y)),length(unique(les.z)));
    
    % u, v, and w are normalized velocity deficits (nondimensional)
    % Kdiff is normalized as well
    
    %% Average across repeating backbones
    % Issue of infinite in cross-flow LES domain versus finite cross-flow MAG
    % domain. Solution applied here is to take the average y-domain across many
    % backbones

    % where backbones repeat on LES grid
    b_repeats = [0:farm.bline_spacing/lesgrid.dy:lesgrid.y/lesgrid.dy]';

    % take the mean in y dimension for same position along growth line
    % (inbetween backbones)
    for b_count = 1:farm.bline_spacing/lesgrid.dy % 13 grid points between backbones
        lesgrid.u_bavg(:,b_count,:) = nanmean(lesgrid.u(:,b_repeats+b_count,:),2);
        lesgrid.v_bavg(:,b_count,:) = nanmean(lesgrid.v(:,b_repeats+b_count,:),2);
        lesgrid.w_bavg(:,b_count,:) = nanmean(lesgrid.w(:,b_repeats+b_count,:),2);
        lesgrid.K_bavg(:,b_count,:) = nanmean(lesgrid.K(:,b_repeats+b_count,:),2);
    end
    clear b_count b_repeats
    
%% LES to MAG - per backbone
% LES grid is different from MAG grid
% interpolate between the two

    % MAG GRID - repeating backbones
    % set with farm; this is a special subset grid of farm.gridx with the
    % mesh-y being just the backbone strip; also include where the kelp
    % starts in the MAG grid (e.g, replace buffer with negative values)
    
        [farm.magmeshx, farm.magmeshb, farm.magmeshz] = ndgrid(0-farm.buffer_horz:farm.dx:farm.x-farm.buffer_horz,1:farm.dy:farm.bline_spacing-1,1:farm.dz:farm.z);
        
    % LES GRID
    % incorporate where the farm starts in LES so that we can align MAG and
    % LES
        lesgrid.les_farmstart = 150; % farm starts at x = 150; y = 0;
        [lesgrid.lesmeshx, lesgrid.lesmeshy, lesgrid.lesmeshz] = ndgrid(unique(les.x)-lesgrid.les_farmstart,0:lesgrid.dy:farm.bline_spacing-lesgrid.dy,flip(abs(unique(les.z))));
        
    % LES -> MAG

        % "bdef" means it is velocity deficit for the "average" backbone
        u_bdef = interpn(lesgrid.lesmeshx,lesgrid.lesmeshy,lesgrid.lesmeshz,lesgrid.u_bavg,farm.magmeshx,farm.magmeshb,farm.magmeshz,'spline');
        v_bdef = interpn(lesgrid.lesmeshx,lesgrid.lesmeshy,lesgrid.lesmeshz,lesgrid.v_bavg,farm.magmeshx,farm.magmeshb,farm.magmeshz,'spline');
        w_bdef = interpn(lesgrid.lesmeshx,lesgrid.lesmeshy,lesgrid.lesmeshz,lesgrid.w_bavg,farm.magmeshx,farm.magmeshb,farm.magmeshz,'spline');
        K_b    = interpn(lesgrid.lesmeshx,lesgrid.lesmeshy,lesgrid.lesmeshz,lesgrid.K_bavg,farm.magmeshx,farm.magmeshb,farm.magmeshz,'spline');
            % just in case interp extrapolates a negative value -> if it
            % does will cause BIG problems in transport function
            K_b(K_b < 0) = min(min(min(lesgrid.K_bavg)));
         
    % Concatenate for number of MAG backbones (repeat the "average" backbone)
    
        % doesn't include buffer in the y-dimension because LES doesn't contain
        % that information (issue of translating between infinite and finite
        % space between LES and MAG y-dimension)
    
        u_def = u_bdef;
        v_def = v_bdef;
        w_def = w_bdef;
        K     = K_b;
        
        for repeat = 1:length(farm.bline_grid)-1
            
            % This is the data product to be applied in transport function
            u_def = cat(2,u_def,u_bdef);
            v_def = cat(2,v_def,v_bdef);
            w_def = cat(2,w_def,w_bdef);
            K     = cat(2,K,K_b);
            
        end
        clear repeat
        clear u_bdef v_bdef w_bdef K_b
        
        
        % Save into the structure envt and identify whether run a
        % ("subsurface) or run b ("max canopy")
        
        % ff=1 full canopy; deep ml
        % ff=2 full canopy; shallow ml
        % ff=3 subsurface; deep ml
        % ff=4 subsurface; shallow ml
        
        if ff == 3 
            
            % preallocate space; zeros because there is a buffer in the
            % cross-shore where information is unavailable from LES so set
            % to zero so that value remains ROMS flow condition
            envt.u_def_a = zeros(size(farm.gridx));
            envt.v_def_a = zeros(size(farm.gridx));
            envt.w_def_a = zeros(size(farm.gridx));
            envt.K_a = zeros(size(farm.gridx));
            
            envt.u_def_a(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = u_def;
            envt.v_def_a(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = v_def;
            envt.w_def_a(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = w_def;
            envt.K_a(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = K;
            
            % fill buffer region with nearest data ...
            for yy = 1:farm.buffer_horz-1
            envt.u_def_a(:,yy,:) = envt.u_def_a(:,farm.buffer_horz,:);
            envt.v_def_a(:,yy,:) = envt.v_def_a(:,farm.buffer_horz,:);
            envt.w_def_a(:,yy,:) = envt.w_def_a(:,farm.buffer_horz,:);
            envt.K_a(:,yy,:) = envt.K_a(:,farm.buffer_horz,:);
            end
            for yy = farm.y:-1:farm.y-farm.buffer_horz
            envt.u_def_a(:,yy,:) = envt.u_def_a(:,farm.y-farm.buffer_horz-1,:);
            envt.v_def_a(:,yy,:) = envt.v_def_a(:,farm.y-farm.buffer_horz-1,:);
            envt.w_def_a(:,yy,:) = envt.w_def_a(:,farm.y-farm.buffer_horz-1,:);
            envt.K_a(:,yy,:) = envt.K_a(:,farm.y-farm.buffer_horz-1,:);
            end
        
            % interp onto transport grid
            envt.u_def_a = interpn(farm.gridx,farm.gridy,farm.gridz,envt.u_def_a,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.v_def_a = interpn(farm.gridx,farm.gridy,farm.gridz,envt.v_def_a,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.w_def_a = interpn(farm.gridx,farm.gridy,farm.gridz,envt.w_def_a,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.K_a = interpn(farm.gridx,farm.gridy,farm.gridz,envt.K_a,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
        
        elseif ff == 4
            
            % preallocate space; zeros because there is a buffer in the
            % cross-shore where information is unavailable from LES so set
            % to zero so that value remains ROMS flow condition
            envt.u_def_b = zeros(size(farm.gridx));
            envt.v_def_b = zeros(size(farm.gridx));
            envt.w_def_b = zeros(size(farm.gridx));
            envt.K_b = zeros(size(farm.gridx));
            
            envt.u_def_b(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = u_def;
            envt.v_def_b(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = v_def;
            envt.w_def_b(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = w_def;
            envt.K_b(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = K;
            
            % fill buffer region with nearest data ...
            for yy = 1:farm.buffer_horz-1
            envt.u_def_b(:,yy,:) = envt.u_def_b(:,farm.buffer_horz,:);
            envt.v_def_b(:,yy,:) = envt.v_def_b(:,farm.buffer_horz,:);
            envt.w_def_b(:,yy,:) = envt.w_def_b(:,farm.buffer_horz,:);
            envt.K_b(:,yy,:) = envt.K_b(:,farm.buffer_horz,:);
            end
            for yy = farm.y:-1:farm.y-farm.buffer_horz
            envt.u_def_b(:,yy,:) = envt.u_def_b(:,farm.y-farm.buffer_horz-1,:);
            envt.v_def_b(:,yy,:) = envt.v_def_b(:,farm.y-farm.buffer_horz-1,:);
            envt.w_def_b(:,yy,:) = envt.w_def_b(:,farm.y-farm.buffer_horz-1,:);
            envt.K_b(:,yy,:) = envt.K_b(:,farm.y-farm.buffer_horz-1,:);
            end    
            
            % interp onto transport grid
            envt.u_def_b = interpn(farm.gridx,farm.gridy,farm.gridz,envt.u_def_b,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.v_def_b = interpn(farm.gridx,farm.gridy,farm.gridz,envt.v_def_b,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.w_def_b = interpn(farm.gridx,farm.gridy,farm.gridz,envt.w_def_b,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.K_b = interpn(farm.gridx,farm.gridy,farm.gridz,envt.K_b,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
        
        elseif ff == 1
            
            % preallocate space; zeros because there is a buffer in the
            % cross-shore where information is unavailable from LES so set
            % to zero so that value remains ROMS flow condition
            envt.u_def_c = zeros(size(farm.gridx));
            envt.v_def_c = zeros(size(farm.gridx));
            envt.w_def_c = zeros(size(farm.gridx));
            envt.K_c = zeros(size(farm.gridx));
            
            envt.u_def_c(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = u_def;
            envt.v_def_c(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = v_def;
            envt.w_def_c(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = w_def;
            envt.K_c(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = K;
            
            % fill buffer region with nearest data ...
            for yy = 1:farm.buffer_horz-1
            envt.u_def_c(:,yy,:) = envt.u_def_c(:,farm.buffer_horz,:);
            envt.v_def_c(:,yy,:) = envt.v_def_c(:,farm.buffer_horz,:);
            envt.w_def_c(:,yy,:) = envt.w_def_c(:,farm.buffer_horz,:);
            envt.K_c(:,yy,:) = envt.K_c(:,farm.buffer_horz,:);
            end
            for yy = farm.y:-1:farm.y-farm.buffer_horz
            envt.u_def_c(:,yy,:) = envt.u_def_c(:,farm.y-farm.buffer_horz-1,:);
            envt.v_def_c(:,yy,:) = envt.v_def_c(:,farm.y-farm.buffer_horz-1,:);
            envt.w_def_c(:,yy,:) = envt.w_def_c(:,farm.y-farm.buffer_horz-1,:);
            envt.K_c(:,yy,:) = envt.K_c(:,farm.y-farm.buffer_horz-1,:);
            end    
            
            % interp onto transport grid
            envt.u_def_c = interpn(farm.gridx,farm.gridy,farm.gridz,envt.u_def_c,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.v_def_c = interpn(farm.gridx,farm.gridy,farm.gridz,envt.v_def_c,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.w_def_c = interpn(farm.gridx,farm.gridy,farm.gridz,envt.w_def_c,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.K_c = interpn(farm.gridx,farm.gridy,farm.gridz,envt.K_c,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
        
        elseif ff == 2
            
            % preallocate space; zeros because there is a buffer in the
            % cross-shore where information is unavailable from LES so set
            % to zero so that value remains ROMS flow condition
            envt.u_def_d = zeros(size(farm.gridx));
            envt.v_def_d = zeros(size(farm.gridx));
            envt.w_def_d = zeros(size(farm.gridx));
            envt.K_d = zeros(size(farm.gridx));
            
            envt.u_def_d(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = u_def;
            envt.v_def_d(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = v_def;
            envt.w_def_d(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = w_def;
            envt.K_d(:,farm.buffer_horz:farm.y-farm.buffer_horz-1,:) = K;
            
            % fill buffer region with nearest data ...
            for yy = 1:farm.buffer_horz-1
            envt.u_def_d(:,yy,:) = envt.u_def_d(:,farm.buffer_horz,:);
            envt.v_def_d(:,yy,:) = envt.v_def_d(:,farm.buffer_horz,:);
            envt.w_def_d(:,yy,:) = envt.w_def_d(:,farm.buffer_horz,:);
            envt.K_d(:,yy,:) = envt.K_d(:,farm.buffer_horz,:);
            end
            for yy = farm.y:-1:farm.y-farm.buffer_horz
            envt.u_def_d(:,yy,:) = envt.u_def_d(:,farm.y-farm.buffer_horz-1,:);
            envt.v_def_d(:,yy,:) = envt.v_def_d(:,farm.y-farm.buffer_horz-1,:);
            envt.w_def_d(:,yy,:) = envt.w_def_d(:,farm.y-farm.buffer_horz-1,:);
            envt.K_d(:,yy,:) = envt.K_d(:,farm.y-farm.buffer_horz-1,:);
            end    
 
            % interp onto transport grid
            envt.u_def_d = interpn(farm.gridx,farm.gridy,farm.gridz,envt.u_def_d,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.v_def_d = interpn(farm.gridx,farm.gridy,farm.gridz,envt.v_def_d,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.w_def_d = interpn(farm.gridx,farm.gridy,farm.gridz,envt.w_def_d,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
            envt.K_d = interpn(farm.gridx,farm.gridy,farm.gridz,envt.K_d,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'spline');
        
        end
        
    end
    
%% Estimate mixed layer depth from ROMS as selection criteria for LES case a, b or c
% CALCULATE Depth of Upper Boundary Layer
% = shallowest depth were AKv is 10% of mean AKv of upper 5 m

    % smooth vertical profile across time series
    for dd=1:size(envt.Dv,2)
        Dvsm(:,dd) = smooth(envt.Dv(:,dd));
    end
    % mean value of AKv in upper 5 metters
    surfDv = nanmean(Dvsm(1:5,:),1);
    bl_logical = Dvsm < surfDv./10;
    for dd=1:size(bl_logical,2)
         dbl_t = find(bl_logical(:,dd), 1, 'first');
         if ~isempty(dbl_t)
             dbl(1,dd) = dbl_t;
         else
             dbl(1,dd) = NaN;
         end
    end
    dbl(dbl<5) = NaN;
    % smooth result in time - two-week running mean
    envt.mld = smooth(dbl,7*6);

end
