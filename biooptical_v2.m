function PAR_field = biooptical_v2(Nf,envt,farm,step)
% Bio-optical model, calculates the attenuation of light due to water,
% chl-a (from ROMS) and nitrogen-specific shading. Incoming PAR from ROMS.
% 
% Input: Nf (per m3; not per frond) * already smoothed at canopy
%        envt, farm, step (ENVT data)
% 
% Output: PAR as a function of depth across the entire farm regardless of
% whether there is kelp present in a given area or not.
%

global param


%% PAR, incoming
PARo = envt.PAR(1,step);

% preallocate space
K=NaN(size(farm.gridx));
PARz=NaN(size(farm.gridx));

% replacement of NaN with zero so that sum of K works in loop    
Nf(isnan(Nf)) = 0; 
     

%% Attenuation of PAR with depth
% Calculate attenuation coefficents and resulting PAR with depth loop

    for z = 1:farm.z_cult
            
        if z==1 
           % no attenuation at surface
           PARz(:,:,z) = PARo;
           
        else
           % attenuate with sum of three contributions
           % 
           K(:,:,z) = param.PAR_Ksw   .* farm.dz...
                    + param.PAR_Kchla .* envt.chla(z,step)...
                    + param.PAR_KNf   .* Nf(:,:,z-1);
           PARz(:,:,z) = PARz(:,:,z-1) .* (exp(-K(:,:,z)));
           
        end
            
    end
       
% Output
PAR_field = PARz;
clear z K PARz

end