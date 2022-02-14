%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  20200805, Christina Frieder
% 
%  MACMODS
%   Full integration of ROMS, LES, MAG
% 
%  mag_sim - model of macroalgal growth in 3-D, a biogeochemical model
%  of nutrient flow through macroalgae. Nitrate, ammonium, urea uptake by
%  macroalgae in a two-step process: first into stored intracellular pools
%  (Ns) and then assimilated into fixed pools (Nf). Both occur at rates
%  dependent on envrionmental factors (nutrient concentrations,
%  hydrodynamics - waves period and current magnitude, temperature,
%  irradiance. Ns is continually lost via exudation to the dissolved
%  organic nitrogen pool (DON), and Nf is continually lost via mortality to
%  the particulate organic nitrogen pool (PON). Both are lost to wave
%  action - significant wave height. 
%
%  Model species: Macrocystic pyrifera 
%  State Variables:
%    NO3, Concentration of nitrate in seawater, [umol NO3/m3]
%    NH4, Concentration of ammonium in seawater, [umol NH4/m3]
%    Ns, macroalgal stored nitrogen, [mg N/m frond]
%    Nf, macroalgal fixed nitrogen, [mg N/m frond]
%    DON, dissolved organic nitrogen, [mmol N/m3]
%    PON, particulate organic nitrogen, [mg N/m3]
%  Farm Design:
%    x, farm lines
%    y, space between lines
%    z, cultivation depth
%  Environmental input:
%    ROMS, 1-km grid solution for SCB
%    NBDC for waves
%    LES for modified velocity fields and vertical diffusivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load MAG parameters
global param
param = param_v1;

%% Directories

    % HPC
    dir_ROMS   = '/data/homezvol2/friederc/mag/SBC_farm/SBCfarm_';
    dir_WAVE   = '/data/homezvol2/friederc/mag/SBC_farm/';
    dir_LES    = '/data/homezvol2/friederc/mag/LES_data/';
    dir_OUTPUT = '/dfs2/davis/friederc/hpc3_output/';
    % Local
%     dir_ROMS   = 'D:\Data\SBC_Farm\SBCfarm_';
%     dir_WAVE   = 'D:\Data\SBC_Farm\';
%     dir_LES    = 'D:\LES_sim\LES_data\';


for year = 2005
    
    % Time 
    
        time = time_v5([year 1 1; year 12 31]);

    % Farm dimensions, plant spacing, cultivation depth
        % There are two grids:
        % 1. Kelp grid
        % 2. Transport grid (coarser to speed up transport function)
    
        farm = farm_v6;
    
    % Environmental Boundary Conditions
        % loads data for duration of simulation (daily)
        % loads LES velocity deficit data for use by transport function
    
        envt = envt_sb_v5(farm,time,dir_ROMS,dir_WAVE); % Santa Barbara 
        envt = envt_les_v11(envt,farm,dir_LES); % load velocity deficit fields
        
    % Initialize kelp
        % Kelp tracked per frond -> kelp_fr
        % Kelp tracked per area -> kelp_ar; sum of kelp_fr per area
    
       [kelp_fr, kelp_ar] = init_v3(farm,time);
       
    % NetCDF -> save variables: frond ID, location, time, Btot, Bcan
    
        file_date = sprintf('MAG_18_05_h.nc');
        output_filename1  = strcat(dir_OUTPUT,file_date); clear file_date
        %file_date = sprintf('MAG_18_00_h_fronds.nc');
        %output_filename2  = strcat(dir_OUTPUT,file_date); clear file_date
        
        create_netcdf_v1(kelp_fr,time,output_filename1)
        %create_netcdf_v2(kelp_fr,time,output_filename2)


%% Run MAG
for growth_step = time.dt_Gr:time.dt_Gr:time.duration

    
    Gr_step = growth_step / time.dt_Gr; % growth counter
    ROMS_step = ceil(Gr_step*time.dt_Gr/time.dt_ROMS); % ROMS counter
        disp(Gr_step);
        disp(datetime('now'));
        
        if Gr_step == 99 % explore doy 99
           output_workspace = strcat(dir_OUTPUT,'MAG_18_05_h_DOY99.mat');
           save(output_workspace,'-v7.3')
        end
            
    
    %% DERIVED BIOLOGICAL CHARACTERISTICS

        [kelp_fr, kelp_ar] = biolchar_v4(kelp_fr,kelp_ar,farm);
        save_netcdf_v1(kelp_fr,time,Gr_step,output_filename1)
        %save_netcdf_v2(kelp_fr,time,Gr_step,output_filename2)
        
        % save frond level data (biomass per frond for validation, doesn't
        % need to be depth resolved)
        
    %% DERIVED ENVT
        % attenuation of light based on ROMS incoming PAR + sw + chla +
        % kelp self-shading
        envt.PAR_field  = biooptical_v2(kelp_ar.Nf,envt,farm,ROMS_step);
        
        % magnitude velocity = ROMS + LES
        envt.magu_field = magu_v9(envt,farm,ROMS_step,'LES',kelp_ar.scenario);
       
        % setup transport function
        % dynamic time step for transport
        time = deltat_v9(kelp_ar.scenario,envt,farm,time,ROMS_step);
        envt.Tr.NO3 = transport_setup_v20(envt,farm,time,ROMS_step,'NO3',kelp_ar.scenario);               
        envt.Tr.NH4 = transport_setup_v20(envt,farm,time,ROMS_step,'NH4',kelp_ar.scenario);               
        envt.Tr.DON = transport_setup_v20(envt,farm,time,ROMS_step,'DON',kelp_ar.scenario);               
    
        
    %% Calc UPTAKE
    
        dd = 0; % counter; time (hours) into the day (dynamic uptake look)
        while dd < 24
        [kelp_fr, uptake_m3] = uptake_setup_v7(kelp_fr,envt,farm,ROMS_step,'transportON','Michaelis+DBL');   
           
           %% TRANSPORT LOOP
           % transport every dt_transport until any nutrient value changes
           % by more than 0.2; then come back to recalculate uptake and
           % start transport again until 24 hours is reached
           [envt.NO3_field, envt.NH4_field, envt.DON_field, duration_dt] = ...
               transport_solver_v7(envt.NO3_field,envt.NH4_field,envt.DON_field,uptake_m3,envt.Tr,kelp_ar,farm,time,dd);
           
           %% calc average uptake
           % because uptake is recalculated every duration_dt -> keep track
           % of the average-weighted value for use by the growth function
           kelp_fr = uptakesum_v1(kelp_fr,dd,duration_dt);
           dd = dd + duration_dt
           
        end   
            

    %% GROWTH MODEL

       [kelp_fr, kelp_ar] = growth_setup_v4(kelp_fr,kelp_ar,envt,farm,time,Gr_step,ROMS_step);

 
    %% Informational output -> track progress
    output_table(Gr_step,:) = table(Gr_step,{kelp_ar.scenario},envt.mld(ROMS_step),time.dt_Tr,...
        [min(min(min(envt.NO3_field(:,:,1:farm.z_cult)))) nanmean(nanmean(nanmean(envt.NO3_field(:,:,1:farm.z_cult)))) max(max(max(envt.NO3_field(:,:,1:farm.z_cult))))],...
        [min(min(min(envt.NH4_field(:,:,1:farm.z_cult)))) nanmean(nanmean(nanmean(envt.NH4_field(:,:,1:farm.z_cult)))) max(max(max(envt.NH4_field(:,:,1:farm.z_cult))))],...
        [min(min(min(envt.DON_field(:,:,1:farm.z_cult)))) nanmean(nanmean(nanmean(envt.DON_field(:,:,1:farm.z_cult)))) max(max(max(envt.DON_field(:,:,1:farm.z_cult))))]);
    output_notes = strcat(dir_OUTPUT,'MAG_18_05_h_output.mat');
    save(output_notes,'output_table','-v7.3')
    
end % end MAG
end % end years of simulation

exit

