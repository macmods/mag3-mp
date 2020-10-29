function [] = save_netcdf_v1(kelp_fr,time,save_step,filename)

%% Save INITIAL Conditions

    % Time step
    ncwrite(filename,'Time',time.timevec_Gr(save_step),save_step)
    
    % generate biomass matrix
    kelplist = fieldnames(kelp_fr);
    
    Btot = NaN(length(kelplist),1); % preallocate
    Bcan = NaN(length(kelplist),1); % preallocate
    
    for iField = 1:numel(fieldnames(kelp_fr))

        kelpid = char(kelplist(iField));
        Btot(iField,:) = nansum(nansum(kelp_fr.(kelpid).B));
        Bcan(iField,:) = nansum(kelp_fr.(kelpid).B(:,1));

    end
    clear iField kelplist kelpid    


    % write to file
    ncwrite(filename,'Btot',Btot,[1 save_step])
    ncwrite(filename,'Bcan',Bcan,[1 save_step])
  
   
end
