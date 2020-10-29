function [kelp_fr, kelp_ar] = biolchar_v4(kelp_fr,kelp_ar,farm)
% biolchar_vX 
% Calculate biological characteristics from Nf, Ns, and Age (known)
% function dependency: type_vX.m, height_vX.m
%
% OUTPUT:  
%   kelp_fr.Q
%   kelp_fr.Biomass
%   kelp_fr.Type
%   kelp_fr.Height, Height_tot
%   kelp_fr.frBlade, fractional biomass that is blade
%   kelp_ar.Nf, assemble data from per frond to per area and smooth in
%   canopy
%
% NOTES:
% Q is integrated across depth, Rassweiler et al. (2018) found that %N does
% not vary with depth. This can be interpreted as translocation occuring on
% the scale of hours (Parker 1965, 1966, Schmitz and Lobban 1976). Side
% note: Ns is redistributed after uptake as a function of fractional Nf.
% This is how the model "translocates" along a gradient of high-to-low
% uptake. Mathematically, this keeps Q constant with depth. -> There may be
% more recent work coming out of SBC LTER that indicates %N varies along
% the frond, particularly in the canopy in a predictable age-like matter
% (e.g., what tissues are doing most of the photosynthesis) (T. Bell pers.
% comm.)
%
% Biomass calculation from Hadley et al. (2015), Table 4 [g(dry)/frond/dz]
%                
% Type of frond: following SBC LTER designation 
%   1 = subsruface
%   2 = canopy
%   3 = senescing
%            
% Blade-to-stipe ratio derived from Nyman et al. 1993 Table 2
  

global param
% preallocate space to sum from per frond to per area
avg_frlength = NaN(numel(fieldnames(kelp_fr)),1);
max_frlength = NaN(numel(fieldnames(kelp_fr)),1);    

%% Calculate DERIVED variables on a per frond basis

kelplist = fieldnames(kelp_fr);
for iField = 1:numel(fieldnames(kelp_fr))

    kelpid = char(kelplist(iField));
    % KNOWN STATE VARIABLES
    % Ns, Nf, Age known 
    % Create temporary variables
    
        Ns = kelp_fr.(kelpid).Ns;
        Nf = kelp_fr.(kelpid).Nf;
        Age= kelp_fr.(kelpid).Age;


    % DERIVED VARIABLES
    
        kelp_fr.(kelpid).Q = param.Qmin .* (1 + nansum(Ns,2) ./ nansum(Nf,2));
        kelp_fr.(kelpid).B = Nf ./ param.Qmin;
        kelp_fr.(kelpid).Type = type_v1(Nf,Age,farm);
        kelp_fr.(kelpid).Height = height_v1(Nf,farm);
        kelp_fr.(kelpid).Height_tot = nansum(kelp_fr.(kelpid).Height,2);
        kelp_fr.(kelpid).Height_tot(kelp_fr.(kelpid).Height_tot == 0) = NaN;

        % Blade to Stipe for blade-specific parameters
        fHeight = cumsum(kelp_fr.(kelpid).Height./kelp_fr.(kelpid).Height_tot,2,'reverse'); % fractional frond height (0 at base; 1 at tip)
        BtoS = param.Blade_stipe(1) - param.Blade_stipe(2) .* fHeight + param.Blade_stipe(3) .* fHeight .^ 2;
        kelp_fr.(kelpid).frBlade = BtoS ./ (BtoS + 1);

        % preallocate
        kelp_fr.(kelpid).UptakeNavg = zeros(farm.frondcount,farm.z_cult/farm.dz);
          
    % PER AREA
    
        kelp_ar.Nf(kelp_fr.(kelpid).kelploc(1),kelp_fr.(kelpid).kelploc(2),1:farm.z_cult) = ...
            nansum(kelp_fr.(kelpid).Nf(:,:),1);

        % Average length of canopy frond for smoothing function
        avg_frlength(iField) = nanmean(kelp_fr.(kelpid).Height(:,1));
        
        % Max frond length for deciding LES velocity deficits
        max_frlength(iField) = max(kelp_fr.(kelpid).Height_tot);
        

end
clear Ns Nf Age fHeight BtoS kelpid kelplist iField
 

%% Redistribute the canopy in the horizontal using smooth2a

    % Rectangular redistribution uniform across farm and based on the
    % average length of canopy-portion within each location

    avg_frlength(isnan(avg_frlength)) = 0;
    kelp_ar.Nf(:,:,1) = smooth2a(kelp_ar.Nf(:,:,1),ceil(nanmean(avg_frlength)));

    % In case MAG is set to 1-D mode
    if farm.x == 1
        if avg_frlength > 1
        kelp_ar.Nf(:,:,1) = param.Nf_capacity_canopy + 0.2 .* kelp_ar.Nf(:,:,1);
        end
    end
    clear avg_frlength
    
%% Identify "Kelp Scenario" for LES velocity deficits

    %
    
    if nanmean(max_frlength) < farm.z_cult/2
        kelp_ar.scenario = char('case_0');
    elseif nanmean(max_frlength) >= farm.z_cult/2 && ...
         nanmean(max_frlength) < farm.z_cult + farm.z_cult * 0.2
        kelp_ar.scenario = char('case_a');
    elseif nanmean(max_frlength) >= farm.z_cult + farm.z_cult * 0.2
        kelp_ar.scenario = char('case_b');
    end
        
end