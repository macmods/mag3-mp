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
%  action - significant wave height. Inorganic nutrients and DON are
%  transported through the farm via a transport function that is informed
%  by boundary conditions (ROMS) and modified velocity and diffusivity
%  fields (LES).
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
%    x, oriented to be parallel to backbone
%    y, 
%    z, depth below sea furface, positive downwards
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
%     dir_ROMS   = '/data/users/friederc/hpc_runs/SBC_farm/SBCfarm_';
%     dir_WAVE   = '/data/users/friederc/hpc_runs/SBC_farm/';
%     dir_LES    = '/data/users/friederc/hpc_runs/LES_data/';
%     dir_OUTPUT = '/dfs2/davis/friederc/mag_v18_output/';

% Local
dir_ROMS   = 'D:\github\mag1-mp\envtl_data\SBCfarm_';
dir_WAVE   = 'D:\github\mag1-mp\envtl_data\';
dir_LES    = 'D:\github\mag1-mp\envtl_data\';
    
    % Time 
    
        % seeding date based on prior analysis of seeding
        % date maximizing output
        % end date of Dec 31
        time = time_v5([2001 11 18; 2002 12 31]);
    
        % add harvest schedule -> informed by Javier
        time.harvest = [datenum([2002 05 01])...
                        datenum([2002 08 01])...
                        datenum([2002 11 01])];
                    
    % Farm dimensions, plant spacing, cultivation depth
        % There are two grids:
        % 1. Kelp grid
        % 2. Transport grid (coarser to speed up transport function)
    
        farm = farm_v6;
    
    % Environmental Boundary Conditions
        % loads data for duration of simulation (daily)
        % loads LES velocity deficit data for use by transport function
    
        envt = envt_sb_v5(farm,time,dir_ROMS,dir_WAVE); % Santa Barbara 
        envt = envt_les_v9(envt,farm,dir_LES); % load velocity deficit fields
        
    % Initialize kelp
        % Kelp tracked per frond -> kelp_fr
        % Kelp tracked per area -> kelp_ar; sum of kelp_fr per area
    
       [kelp_fr, kelp_ar] = init_v3(farm,time);
       
    % NetCDF -> save variables: frond ID, location, time, Btot, Bcan
    
        file_date = 'MAG_18_02_d.nc';
        output_filename  = strcat(dir_OUTPUT,file_date); clear file_date
        create_netcdf_v1(kelp_fr,time,output_filename)

    %% Run MAG
    for growth_step = time.dt_Gr:time.dt_Gr:time.duration

        Gr_step = growth_step / time.dt_Gr; % growth counter
        ROMS_step = ceil(Gr_step*time.dt_Gr/time.dt_ROMS); % ROMS counter
            disp(datetime('now'));
            disp(Gr_step);

        %% DERIVED BIOLOGICAL CHARACTERISTICS

            [kelp_fr, kelp_ar] = biolchar_v4(kelp_fr,kelp_ar,farm);
            save_netcdf_v1(kelp_fr,time,Gr_step,output_filename)

        %% DERIVED ENVT
            % attenuation of light based on ROMS incoming PAR + sw + chla +
            % kelp self-shading
            envt.PAR_field  = biooptical_v2(kelp_ar.Nf,envt,farm,ROMS_step);

            % magnitude velocity = ROMS + LES
            envt.magu_field = magu_v7(envt,farm,ROMS_step,'LES',kelp_ar.scenario);

            % setup transport function
            % dynamic time step for transport
            time = deltat_v9(kelp_ar.scenario,envt,farm,time,ROMS_step);
            envt.Tr.NO3 = transport_setup_v18(envt,farm,time,ROMS_step,'NO3',kelp_ar.scenario);               
            envt.Tr.NH4 = transport_setup_v18(envt,farm,time,ROMS_step,'NH4',kelp_ar.scenario);               
            envt.Tr.DON = transport_setup_v18(envt,farm,time,ROMS_step,'DON',kelp_ar.scenario);               


        %% TRANSOPRT & UPTAKE; NESTED LOOPS
        % transport as long as delta N is < threshold, then recalc uptake and
        % do again as long as duration is less than time.dt_Gr

            dd = 0; % counter; time (hours) into the day (dynamic uptake look)
            while dd < time.dt_Gr
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
               dd = dd + duration_dt %#ok<NOPTS>

            end   


        %% GROWTH MODEL

           [kelp_fr, kelp_ar] = growth_setup_v4(kelp_fr,kelp_ar,envt,farm,time,Gr_step,ROMS_step);


        %% HARVEST
        % harvest is based on a prescribed schduled informed by Javier

            % does current date match any of the harvest dates. Yes -> harvest
           if sum(ismember(time.harvest,time.timevec_Gr(Gr_step))) == 1

               % which harvest number is it (for structure id); save date
               h_no = find(ismember(time.harvest,time.timevec_Gr(Gr_step)));
               h_no = sprintf('h%d',h_no);
               harvest.(h_no).date = time.timevec_Gr(Gr_step);

               % chop off the canopy and store it in harvest structure
               % update kelp_fr with
                    % 1. removed Nf and Ns at surface
                    % 2. update "age = 9999" so that cut fronds are identified
                    % as senescing
                    [harvest.(h_no).Nf, kelp_fr] = harvest_v1(kelp_fr,farm);

                    output_harvest = strcat(dir_OUTPUT,'MAG_18_02_d_harvest.mat');
                    save(output_harvest,'harvest','-v7.3')

            end
       
       
    %% Informational output -> track progress
    
       %output_table(Gr_step,:) = ...
       %output_notes = strcat(dir_OUTPUT,'MAG_18_02_d_output.mat');
       %save(output_notes,'output_table','-v7.3')
    
    end % end mag loop

exit

