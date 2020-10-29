function farm = farm_v6
% Farm Properties
% Output: farm.()
%   buffer; space on each side of farm (including below cult.), [m]
%   z_cult; depth of cultivation, [m]
%   x,y,z; dimension of farm [m]
%   dx,dy,dz; bin size [m] 
%   dx_tr,dy_tr,dz_tr; bin size for transport function [m]
%   bline_spacing; spacing between backbones [m]
%   gline_spacing; spacing between growth lines [m]
%   gline_length; length of growth lines [m]
%   rotation; rotation of farm -> used to align dominant flow direction
%   parallel to backbone
%   frondcount = used to preallocate space for tracking fronds

%% farm dimensions and grid

        farm.buffer_horz = 20; % buffer for canopy smoothing
        farm.buffer_vert = 10; % buffer below farm structure
        
    farm.z_cult = 20; %[m below surface]
    farm.x      = 400+farm.buffer_horz*2; % [m]
    farm.y      = 400+farm.buffer_horz*2; % [m]
    farm.z      = farm.z_cult+farm.buffer_vert; % [m]
    farm.dx     = 1; % [m]
    farm.dy     = 1; % [m]
    farm.dz     = 1; % [m]
    farm.dx_tr  = farm.x ./ (0.1 .* farm.x); % [m]
    farm.dy_tr  = farm.y ./ (0.1 .* farm.y); % [m]
    farm.dz_tr  = farm.z ./ (1.0 .* farm.z); % [m]
        
 
%% farm grid

    %farm grid -> %transport grid
    
        [farm.gridx,farm.gridy,farm.gridz] = ndgrid(0:farm.dx:farm.x,0:farm.dy:farm.y,1:farm.dz:farm.z);
        [farm.gridx_tr,farm.gridy_tr,farm.gridz_tr] = ndgrid(0:farm.dx_tr:farm.x,0:farm.dy_tr:farm.y,1:farm.dz_tr:farm.z);

%% plant spacing

    farm.bline_spacing = 26;
    farm.gline_spacing = 2;
    farm.gline_length  = 8; % [m]
    
    farm.bline_grid = farm.buffer_horz:farm.bline_spacing:farm.y-farm.buffer_horz; % [m]
    farm.gline_grid = farm.buffer_horz:farm.gline_spacing:farm.x-farm.buffer_horz; % [m]
    
%% farm rotation

    % rotate farm to be parallel/perpendicular with dominant flow direction
    % -> based on pre-analysis with ROMS data for a given location; this is
    % important to align with LES
    % positive rotation to match an alignment of left-to-right
    % visualization
    
    farm.rotation = 133; % degrees

%% track fronds (preallocate and initialize)

    farm.frondcount = 60;

end
