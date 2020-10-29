function [Tr] = transport_setup_v18(envt,farm,time,envt_step,arg_in1,arg_in2)
% Generates necessary variables to carry-out transport equations
%
% INPUT:  envt, farm, time, step (ENVT), arg_in1, arg_in2, arg_in3
%         arg_in1 = N species; 'NO3' 'NH4' 'DON' arg_in2 = kelp scenario:
%         case_0 (no drag); case_a ('subsurface'); case_b ('max canopy;
%         deep mixed layer'); case_c ('max canopy; shallow mixed layer')
% OUTPUT: Tr.Nutrient
%         used by transport function


%% Advection [m/h], Diffusion [m2/h]
% Select advection and diffusion fields across farm based on kelp scenario
% (arg_in2)
        
    % Advection
    
        ux = permute(repmat(envt.u(:,envt_step),1,size(farm.gridx_tr,1),size(farm.gridx_tr,2)),[2 3 1]);
        vy = permute(repmat(envt.v(:,envt_step),1,size(farm.gridx_tr,1),size(farm.gridx_tr,2)),[2 3 1]);
        wz = permute(repmat(envt.w(:,envt_step),1,size(farm.gridx_tr,1),size(farm.gridx_tr,2)),[2 3 1]);
              
    % Diffusion
    % envt.Dh = horizontal diffusion
    % Dv = vertical diffusion
    
        % constant; same for LES and ROMS
        Dh = repmat(envt.Dh,size(farm.gridx_tr,1),size(farm.gridx_tr,2),size(farm.gridx_tr,3));
        
        if strcmp(arg_in2,'case_0')
            
            %Don't modify flows
            
            %use ROMS Dv; vertically averaged
            Dv = repmat(nanmean(envt.Dv(:,envt_step)),size(farm.gridx_tr,1),size(farm.gridx_tr,2),size(farm.gridx_tr,3));
            
        elseif strcmp(arg_in2,'case_a')
            
            if envt.mld(envt_step) >= farm.z_cult
            %modify ROMS flow with LES "subsurface" velocity deficits
            ux = envt.u_def_a .* ux + ux;
            vy = envt.v_def_a .* vy + vy;
            wz = envt.w_def_a .* abs(nanmean(envt.u(:,envt_step)));
            
            % u-star = 0.0061 m/s * 60 * 60 = 21.96 m/h
            % u-start is surface friction velocity -> this should probably
            % be dynamic based on ROMS B.C.
            Dv = envt.K_a .* (21.96 .* farm.z_cult); % normalized by u* and zi [m2/s];
            
                % scale to ROMS Dv
                ROMS_to_LES = nanmean(envt.Dv(:,envt_step)) ./ nanmean(nanmean(nanmean(Dv)));
                Dv = Dv .*  ROMS_to_LES;
            
            elseif envt.mld(envt_step) < farm.z_cult
            ux = envt.u_def_b .* ux + ux;
            vy = envt.v_def_b .* vy + vy;
            wz = envt.w_def_b .* abs(nanmean(envt.u(:,envt_step)));
            
            % u-star = 0.0061 m2/s * 60 * 60 = 21.96 m2/h
            Dv = envt.K_b .* (21.96 .* farm.z_cult); % normalized by u* and zi [m2/s];
            
                % scale to ROMS Dv
                ROMS_to_LES = nanmean(envt.Dv(:,envt_step)) ./ nanmean(nanmean(nanmean(Dv)));
                Dv = Dv .*  ROMS_to_LES;
            end
            
        elseif strcmp(arg_in2,'case_b')
            
            if envt.mld(envt_step) >= farm.z_cult
            %modify ROMS flow with LES "max canopy" velocity deficits
            ux = envt.u_def_c .* ux + ux;
            vy = envt.v_def_c .* vy + vy;
            wz = envt.w_def_c .* abs(nanmean(envt.u(:,envt_step)));
            
            % u-star = 0.0061 m2/s * 60 * 60 = 21.96 m2/h
            Dv = envt.K_c .* (21.96 .* farm.z_cult); % normalized by u* and zi [m2/s];
            
                % scale to ROMS Dv
                ROMS_to_LES = nanmean(envt.Dv(:,envt_step)) ./ nanmean(nanmean(nanmean(Dv)));
                Dv = Dv .*  ROMS_to_LES;
                   
            elseif envt.mld(envt_step) < farm.z_cult
            %modify ROMS flow with LES "max canopy" velocity deficits
            ux = envt.u_def_d .* ux + ux;
            vy = envt.v_def_d .* vy + vy;
            wz = envt.w_def_d .* abs(nanmean(envt.u(:,envt_step)));
            
            % u-star = 0.0061 m2/s * 60 * 60 = 21.96 m2/h
            Dv = envt.K_d .* (21.96 .* farm.z_cult); % normalized by u* and zi [m2/s];
            
                % scale to ROMS Dv
                ROMS_to_LES = nanmean(envt.Dv(:,envt_step)) ./ nanmean(nanmean(nanmean(Dv)));
                Dv = Dv .*  ROMS_to_LES;
            end
                                
        end

        % flip LES fields if transport coming in from opposite end of farm
        
            % Direction of advection
            direction = nanmean(atan2d(envt.v(:,envt_step),envt.u(:,envt_step)));

            if direction > 90 || direction < -90
                ux = flip(ux,1);
                vy = flip(vy,1);
                wz = flip(wz,1);
                Dv = flip(Dv,1);
            end
            
        %diffusion plus 1/2 step grid; minus 1/2 step grid
        stepx = 1/2 * farm.dx_tr;
        stepy = 1/2 * farm.dy_tr;
        stepz = 1/2 * farm.dz_tr;
        
        [gridDxPlus,~,~] = ndgrid(0+stepx:farm.dx_tr:farm.x+stepx,0:farm.dy_tr:farm.y,1:farm.dz_tr:farm.z);
        [gridDxMinus,~,~] = ndgrid(0-stepx:farm.dx_tr:farm.x-stepx,0:farm.dy_tr:farm.y,1:farm.dz_tr:farm.z);
        [~,gridDyPlus,~] = ndgrid(0:farm.dx_tr:farm.x,0+stepy:farm.dy_tr:farm.y+stepy,1:farm.dz_tr:farm.z);
        [~,gridDyMinus,~] = ndgrid(0:farm.dx_tr:farm.x,0-stepy:farm.dy_tr:farm.y-stepy,1:farm.dz_tr:farm.z);
        [~,~,gridDzPlus] = ndgrid(0:farm.dx_tr:farm.x,0:farm.dy_tr:farm.y,1+stepz:farm.dz_tr:farm.z+stepz);
        [~,~,gridDzMinus] = ndgrid(0:farm.dx_tr:farm.x,0:farm.dy_tr:farm.y,1-stepz:farm.dz_tr:farm.z-stepz);
        
        clear stepx stepy stepz

           

%% TDMA: Generate tridiagonal matrices for LHS of ADI
% Constants -> used by two separate functions: TDMA and ADI
        
    % Advection

        Ax = ux .* time.dt_Tr ./ (2 .* farm.dx_tr);
        Ay = vy .* time.dt_Tr ./ (2 .* farm.dy_tr);
        Az = wz .* time.dt_Tr ./ (2 .* farm.dz_tr);

        % Because applying upwind scheme, as stored in above variables, set
        % velocities to be positive regardless of direction

        Tr.absAx = abs(Ax);
        Tr.absAy = abs(Ay);
        Tr.absAz = abs(Az);

        Tr.absAx_pos = double(Ax > 0) .* Tr.absAx; % upwind
        Tr.absAx_neg = double(Ax < 0) .* Tr.absAx; % downwind
        Tr.absAy_pos = double(Ay > 0) .* Tr.absAy; % upwind
        Tr.absAy_neg = double(Ay < 0) .* Tr.absAy; % downwind
        Tr.absAz_pos = double(Az > 0) .* Tr.absAz; % upwind
        Tr.absAz_neg = double(Az < 0) .* Tr.absAz; % downwind

    % Diffusion

        % original grid
        Tr.Bx = Dh .* time.dt_Tr ./ (farm.dx_tr .^ 2);
        Tr.By = Dh .* time.dt_Tr ./ (farm.dy_tr .^ 2);
        Tr.Bz = Dv .* time.dt_Tr ./ (farm.dz_tr .^ 2);
        % plus 1/2 step
        Tr.BxPlus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dh,gridDxPlus,farm.gridy_tr,farm.gridz_tr,'spline') .* time.dt_Tr ./ (farm.dx_tr .^ 2);
        Tr.ByPlus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dh,farm.gridx_tr,gridDyPlus,farm.gridz_tr,'spline') .* time.dt_Tr ./ (farm.dy_tr .^ 2);
        Tr.BzPlus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dv,farm.gridx_tr,farm.gridy_tr,gridDzPlus,'spline') .* time.dt_Tr ./ (farm.dz_tr .^ 2);
        % minus 1/2 step
        Tr.BxMinus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dh,gridDxMinus,farm.gridy_tr,farm.gridz_tr,'spline') .* time.dt_Tr ./ (farm.dx_tr .^ 2);
        Tr.ByMinus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dh,farm.gridx_tr,gridDyMinus,farm.gridz_tr,'spline') .* time.dt_Tr ./ (farm.dy_tr .^ 2);
        Tr.BzMinus = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,Dv,farm.gridx_tr,farm.gridy_tr,gridDzMinus,'spline') .* time.dt_Tr ./ (farm.dz_tr .^ 2);
            % if Bz is < 0; replace with 0; as diffusion approaches
            % surface diffusivity approaches zero but this interpn function
            % will extrapolate less than zero at z(-1/2) at the surface bin
            Tr.BzMinus(Tr.BzMinus < min(min(min(Tr.Bz)))) = min(min(min(Tr.Bz)));
            Tr.BzPlus(Tr.BzPlus < min(min(min(Tr.Bz)))) = min(min(min(Tr.Bz)));

            
%% ADI Solver
% Requires ghost matrices with boundary conditions
% Extra row at beginning and end of x and y-dimensions and z=end

%preallocate space
Tr.Nt = NaN(size(farm.gridx_tr,1)+2,size(farm.gridy_tr,2)+2,size(farm.gridz_tr,3)+1);


     if strcmp(arg_in1,'NO3')

         BC = envt.NO3(:,envt_step);
         % Generate tridiagonal matrices for ADI solver; only has
         % to be done once per boundary condition

         [Tr.A1, Tr.A2, Tr.A3, Tr.D1, Tr.D2, Tr.D3] = TDMA_v4(Tr, BC, farm); 

         % t* -> this is the first ghost matrix

         Tr.Nt(:,:,1:farm.z) = permute(repmat(BC,1,size(farm.gridx_tr,1)+2,size(farm.gridx_tr,2)+2),[2 3 1]);
         Tr.Nt(:,:,farm.z/farm.dz_tr+1) = (Tr.Nt(:,:,farm.z/farm.dz)-Tr.Nt(:,:,farm.z-1/farm.dz)) + Tr.Nt(:,:,farm.z/farm.dz);


     elseif strcmp(arg_in1,'NH4')

         BC = envt.NH4(:,envt_step);

         % Generate tridiagonal matrices for ADI solver; only has
         % to be done once per boundary condition

         [Tr.A1, Tr.A2, Tr.A3, Tr.D1, Tr.D2, Tr.D3] = TDMA_v4(Tr, BC, farm); 

         % t*
         Tr.Nt(:,:,1:farm.z) = permute(repmat(BC,1,size(farm.gridx_tr,1)+2,size(farm.gridx_tr,2)+2),[2 3 1]);
         Tr.Nt(:,:,farm.z/farm.dz_tr+1) = (Tr.Nt(:,:,farm.z/farm.dz)-Tr.Nt(:,:,farm.z-1/farm.dz)) + Tr.Nt(:,:,farm.z/farm.dz);


      elseif strcmp(arg_in1,'DON')

         BC = envt.DON(:,envt_step);

         % Generate tridiagonal matrices for ADI solver; only has
         % to be done once per boundary condition

         [Tr.A1, Tr.A2, Tr.A3, Tr.D1, Tr.D2, Tr.D3] = TDMA_v4(Tr, BC, farm); 

         % t*
         Tr.Nt(:,:,1:farm.z) = permute(repmat(BC,1,size(farm.gridx_tr,1)+2,size(farm.gridx_tr,2)+2),[2 3 1]);
         Tr.Nt(:,:,farm.z/farm.dz_tr+1) = (Tr.Nt(:,:,farm.z/farm.dz)-Tr.Nt(:,:,farm.z-1/farm.dz)) + Tr.Nt(:,:,farm.z/farm.dz);

     end

     % t**; 2nd ghost matrix for ADI step 2
     Tr.N_LHS1t = Tr.Nt;

     % t***; 3rd ghost matrix for ADI step 3
     Tr.N_LHS2t = Tr.Nt;

end