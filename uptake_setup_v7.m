function [kelp_fr, uptake_m3] = uptake_setup_v7(kelp_fr,envt,farm,envt_step,arg_in1,arg_in2)
% setup for uptake calculation

if strcmp(arg_in1,'transportON')
% pre-allocate space to Uptake variables
uptake_m3.NO3 = zeros(size(farm.gridx));
uptake_m3.NH4 = zeros(size(farm.gridx));
uptake_m3.DON = zeros(size(farm.gridx));
avg_frlength = NaN(numel(fieldnames(kelp_fr)),1);
end

%% Calculate uptake rate (NO3 + NH4 + DON-Urea) per frond
kelplist = fieldnames(kelp_fr);
for iField = 1:numel(fieldnames(kelp_fr))

    kelpid = char(kelplist(iField));
    % Calculate uptake on a per frond basis
    % This is because Q within an area can vary among
    % fronds

        % Save amount of uptake per uptake_step to be transferred to growth
        % function
        if strcmp(arg_in2,'Michaelis')
            
            [kelp_fr.(kelpid).UptakeN, UptakeFactor] = ...
                uptake_v2(kelp_fr.(kelpid),envt,farm,envt_step,arg_in1);
            
        elseif strcmp(arg_in2,'Michaelis+DBL')
            
        [kelp_fr.(kelpid).UptakeN, UptakeFactor] = ...
            uptake_v1(kelp_fr.(kelpid),envt,farm,envt_step,arg_in1);
        
        end
        
        if strcmp(arg_in1,'transportON')
        % For transport function need to combine infro from per frond to
        % per area. Must be done for each nitrogen species; [umol N/m3/h]
        UptakeNO3 = UptakeFactor.Uptake_NO3_mass .* kelp_fr.(kelpid).B; 
            UptakeNO3(isnan(UptakeNO3)) = 0; % per frond
        UptakeNH4 = UptakeFactor.Uptake_NH4_mass .* kelp_fr.(kelpid).B; 
            UptakeNH4(isnan(UptakeNH4)) = 0; % per frond
        UptakeDON = UptakeFactor.Uptake_DON_mass .* kelp_fr.(kelpid).B; 
            UptakeDON(isnan(UptakeDON)) = 0; % per frond


        % Translate onto FARM grid space -> for use by transport function
        uptake_m3.NO3(kelp_fr.(kelpid).kelploc(1),kelp_fr.(kelpid).kelploc(2),1:farm.z_cult) = ...
            nansum(UptakeNO3,1);
        uptake_m3.NH4(kelp_fr.(kelpid).kelploc(1),kelp_fr.(kelpid).kelploc(2),1:farm.z_cult) = ...
            nansum(UptakeNH4,1);
        uptake_m3.DON(kelp_fr.(kelpid).kelploc(1),kelp_fr.(kelpid).kelploc(2),1:farm.z_cult) = ...
            nansum(UptakeDON,1);

        % Canopy length for smoothing below
        avg_frlength(iField) = nanmean(kelp_fr.(kelpid).Height(:,1));
        end
        
end
clear iField kelpid kelplist UptakeFactor UptakeNO3 UptakeNH4 UptakeDON

  
%% Redistribute the canopy in the horizontal using smooth2a

    % Rectangular redistribution uniform across farm and based on the
    % average length of canopy-portion within each location
    
    if strcmp(arg_in1,'transportON')
    avg_frlength(isnan(avg_frlength)) = 0;
    uptake_m3.NO3(:,:,1) = smooth2a(uptake_m3.NO3(:,:,1),ceil(nanmean(avg_frlength)));
    uptake_m3.NH4(:,:,1) = smooth2a(uptake_m3.NH4(:,:,1),ceil(nanmean(avg_frlength)));
    uptake_m3.DON(:,:,1) = smooth2a(uptake_m3.DON(:,:,1),ceil(nanmean(avg_frlength)));
    clear avg_frlength
    end
                  
                            
end                    