function Height = height_v1(Nf,farm)
% Height, allometric conversion from Nf to height [meters]
% Calculation based on allometric relationships between Nf_capacity and
% length (e.g., mg Nf per meter). Type of frond dependent (subsurface,
% canopy)
%
% Output: Height per frond per m

global param
% preallocate space
Height = NaN(size(Nf)); % preallocate variable space


%% Height
% Approximate height; all bins with Nf

    appH = sum(~isnan(Nf(:,:)),2);

% FRONDS < farm.dz

    Height(appH == 1,farm.z_cult) = Nf(appH == 1,farm.z_cult) ./ param.Nf_capacity_subsurface .* farm.dz;
    
    
% SUBSURFACE

    for ss = 2:farm.z_cult-1
        if any(appH == ss)
        % below growing tip; length is length of z-bin
        Height(appH == ss,farm.z_cult-ss+2:farm.z_cult) = farm.dz;
        % at growing tip; length is based on amount of Nf in z-bin
        Height(appH == ss,farm.z_cult-ss+1) = Nf(appH == ss,farm.z_cult-ss+1) ./ param.Nf_capacity_subsurface .* farm.dz;
        end
    end
    clear ss
    
    
% CANOPY

    % subsurface portion
    Height(appH == farm.z_cult,2:farm.z_cult) = farm.dz;
    % canopy portion
    Height(appH == farm.z_cult,1) = Nf(appH == farm.z_cult,1) ./ param.Nf_capacity_canopy;
    clear appH  


end



          