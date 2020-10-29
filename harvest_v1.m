function [harvest_Nf, kelp_fr] = harvest_v1(kelp_fr,farm)
% Calculate harvest per frond and sum across area
% Output: harvest structure
%   Nf harvest per area (x, y)

% preallcoate harvest Nf
harvest_Nf = NaN(farm.x,farm.y);

kelplist = fieldnames(kelp_fr);
for iField = 1:numel(fieldnames(kelp_fr))

    kelpid = char(kelplist(iField));
    
    % sum all Nf in canopy per location
    harvest_Nf(kelp_fr.(kelpid).kelploc(1),kelp_fr.(kelpid).kelploc(2)) = ...
            nansum(kelp_fr.(kelpid).Nf(:,1));
        
    % identify which fronds were in canopy and replace age with 9999
    cut_i = kelp_fr.(kelpid).Nf(:,1) > 0;
    kelp_fr.(kelpid).Age(cut_i) = 9999;
    clear cut_i
    
    % replace cut with NaN in canopy
    kelp_fr.(kelpid).Nf(:,1) = NaN;
    kelp_fr.(kelpid).Ns(:,1) = NaN;
    
end
clear kelplist iField kelpid


end