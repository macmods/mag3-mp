function [kelp_fr, kelp_ar] = growth_setup_v4(kelp_fr,kelp_ar,envt,farm,time,growth_step,envt_step)

%% Calculate growth per frond
kelplist = fieldnames(kelp_fr);
for iField = 1:numel(fieldnames(kelp_fr))

    kelpid = char(kelplist(iField));
    
    % Calculate growth on a per frond basis
    [kelp_fr.(kelpid), kelp_ar] = growth_v5(kelp_fr.(kelpid),kelp_ar,envt,farm,time,envt_step,growth_step);
    
end
clear kelplist iField kelpid


end