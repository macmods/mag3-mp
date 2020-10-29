function envt = envt_sb_v5(farm,time,dir_ROMS,dir_WAVE)
% Environmental Parameter Input Data
% INPUT: farm, time, directories for ROMS and WAVE data
% OUTPUT: envt.()
%   NO3; daily ROMS Nitrate in seawater; [umol/m3]
%   NH4; daily ROMS Ammonium in seawater; [umol/m3]
%   DON; daily ROMS DON, dissolved organic nitrogen; [mmol N/m3]
%   T; daily ROMS Temperature; [Celcius]
%   u; daily ROMS Uo Seawater velocity, u x-direction; [m/h]
%   v; daily ROMS Vo Seawater velocity, v y-direction; [m/h]
%   w; daily ROMS Wo Seawater velocity, w z-direction; [m/h]
%   Dv; daily ROMS vertical diffusivity; [m2/h]
%   Dh; horizontal diffusivity; constant [m2/h]
%   PAR; daily ROMS PAR, photosynthetically active radiation; incoming PAR; [W/m2]
%   chla; adily ROMS chl-a: sum of DIAZ+DIAT+SP (small phytoplankton); [mg-chla/m3]
%   Tw; daily NDCP Wave period; [h]
%   Hs; daily NDCP Significant wave height; [m]
%
% The ROMS simulations have been pre-processed and saved as mat files for
% use by MAG.
% The NDCB data has been downloaded and saved as mat files for use by MAG.

%% File Directory
           
% SBC farm site; offshore Mohawk; 60-m water depth

    ROMS_start = datenum([1994 01 01 0 0 0]); % from ROMS folks
    
% Extract days

    filename = strcat(dir_ROMS,'NO3.mat'); NO3 = load(filename);
    ROMS_time = datevec(ROMS_start + NO3.ocean_time ./ (60*60*24)); clear NO3
    ROMS_time = datenum(ROMS_time(:,1:3));
    
    start = min(time.timevec_ROMS);
    stop  = max(time.timevec_ROMS);
    [~,idx_start]=min(abs(ROMS_time-start));
    [~,idx_stop] =min(abs(ROMS_time-stop));
        
    ROMS_extract  = idx_start:idx_stop;
    clear val start stop filename ROMS_time idx_start idx_stop
    
    
% Load ROMS data and create variable to store boundary condition.
%% Nitrate 
            
    filename = strcat(dir_ROMS,'NO3.mat');
    NO3 = load(filename);

    envt.NO3 = NO3.NO3(ROMS_extract,1:farm.z)';
    envt.NO3 = envt.NO3 .* 1e3; % umol/m3
    envt.NO3(envt.NO3 <= 0.01e3) = 0.01e3; % replace negatives

        clear NO3 filename

        
%% Ammonium
        
    filename = strcat(dir_ROMS,'NH4.mat');
    NH4 = load(filename);

    envt.NH4 = NH4.NH4(ROMS_extract,1:farm.z)';
    envt.NH4 = envt.NH4 .* 1e3; % umol/m3

        clear NH4 filename

        
%% DON
        
    filename = strcat(dir_ROMS,'DON.mat');
    DON = load(filename);
    
    envt.DON = DON.DON(ROMS_extract,1:farm.z)';
    
        clear DON filename


%% Temperature
    
    filename = strcat(dir_ROMS,'temp.mat');
    temp = load(filename);
    
    envt.T = temp.temp(ROMS_extract,1:farm.z)';
    
        clear temp filename

            
%% Seawater Velocity, u,v,w
% Seawater velocity in x
            
    filename = strcat(dir_ROMS,'u.mat');
    u = load(filename);
    
    envt.u= u.u(ROMS_extract,1:farm.z)';
    envt.u = envt.u .* 60 .* 60; % [m/h]
    
        clear u filename

% Seawater velocity in y
            
    filename = strcat(dir_ROMS,'v.mat');
    v = load(filename);
    
    envt.v = v.v(ROMS_extract,1:farm.z)';
    envt.v = envt.v .* 60 .* 60; % [m/h]
    
        clear v filename
            
% ROTATE velocity fields based on farm orientation

    R = [cosd(farm.rotation) -sind(farm.rotation); sind(farm.rotation) cosd(farm.rotation)];

    for step = 1:size(envt.u,2)
        magu = [envt.u(:,step) envt.v(:,step)];

        vR = magu * R;

        envt.u(:,step) = vR(:,1);
        envt.v(:,step) = vR(:,2);
    end
    clear R
    
% Seawater velocity in z
% This is informed by ROMS and based on conversation with Kristen on
% 20191015 going to set vertical velocities to zero

    filename = strcat(dir_ROMS,'w.mat');
    w = load(filename);
    
    envt.w = w.w(ROMS_extract,1:farm.z)';
    envt.w = envt.w .* 60 .* 60; % [m/h]
    envt.w = envt.w .* 0;
    
        clear w_BC filename


%% Diffusivity

    % Vertical diffusivity provided by ROMS
    
    filename = strcat(dir_ROMS,'AKv.mat');
    AKv = load(filename);
    
    envt.Dv = AKv.AKv(ROMS_extract,1:farm.z)';
    envt.Dv = envt.Dv .* 60 .* 60; % [m2/h]
    
    % Horizontal diffusivity -> still TBD; using placeholder as per
    % conversation with LES and ROMS
    
    % McWilliams: "In a general way one can expect values of 0.1-1 m^2/s on these
    % scales.  We can look at some ROMS dispersion analyses on the shelf to
    % try to be more specific.  The LES is not informative about this
    % number as it lacks the submesoscale and mesoscale influences.:

    envt.Dh = 0.02 .* 60 .* 60;
    
    
%% PAR

    % Load ROMS data of PAR and extract surface value only
    % Also load PARincoming which isn't all that different from PAR at
    % the surface. PARincoming is 0.45 * penetration of solar heat
    % PAR at surface is modified as = PARinc * (1 - exp(Kpar))
    % [W/m2]
        
    filename = strcat(dir_ROMS,'PAR.mat');
    PAR = load(filename);
    
    envt.PAR = PAR.PAR(ROMS_extract,1)';
    
        clear PAR filename
           
            
%% CHL-a

    filename1 = strcat(dir_ROMS,'DIATCHL.mat');
    filename2 = strcat(dir_ROMS,'DIAZCHL.mat');
    filename3 = strcat(dir_ROMS,'SPCHL.mat');
    DIATCHL = load(filename1);
    DIAZCHL = load(filename2);
    SPCHL = load(filename3);

    envt.chla = ...
          DIATCHL.DIATCHL(ROMS_extract,1:farm.z)'...
        + DIAZCHL.DIAZCHL(ROMS_extract,1:farm.z)'...
        + SPCHL.SPCHL(ROMS_extract,1:farm.z)';
    envt.chla(envt.chla < 0) = 0; % replace negatives
    
        clear DIATCHL DIAZCHL SPCHL filename1 filename2 filename3
            
            
%% Wave period, Significant wave height
NDCPfilename = strcat(dir_WAVE,'NDBC46053/NDCP46053.mat');

    % NDCP; mat file
    % Column 1 = matlab time vec
    % Column 2 = Hs
    % Column 3 = Tw
    wave = load(NDCPfilename);
     
    envt.Tw = interp1(wave.NDCP46053(:,1),wave.NDCP46053(:,3),time.timevec_ROMS);
    envt.Tw = envt.Tw ./ (60*60); % [h]
    envt.Hs = interp1(wave.NDCP46053(:,1),wave.NDCP46053(:,2),time.timevec_ROMS);

       clear wave
         
%% INITIALIZE NUTRIENT FIELDS

    envt.NO3_field = permute(repmat(envt.NO3(:,1),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
    envt.NH4_field = permute(repmat(envt.NH4(:,1),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
    envt.DON_field = permute(repmat(envt.DON(:,1),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
        
        clear ROMS_dir ROMS_start ROMS_extract 

end