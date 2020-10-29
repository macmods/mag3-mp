function kelp_fr = uptakesum_v1(kelp_fr,duration,duration_dt)
% average out uptake ...

kelplist = fieldnames(kelp_fr);
for iField = 1:numel(fieldnames(kelp_fr))

    kelpid = char(kelplist(iField));
    kelp_fr.(kelpid).UptakeNavg = ...
        ((kelp_fr.(kelpid).UptakeNavg * duration) + (kelp_fr.(kelpid).UptakeN * duration_dt)) ./ (duration + duration_dt);
    
end
       

end