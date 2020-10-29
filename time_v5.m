function time = time_v5(duration)
% Time Domain, [h]
% Output: time.()
%   start = date of start in date vec [yr,m,d]
%   duration = duration of simulation in hours
%       accounts for different years being different lengths (e.g., leap
%       years)
%   dt_xx = time step for 
%       Tr = transport
%       Up = uptake
%       Gr = growth (MAG model)
%       ROMS = environmental input
%   timevec = time domain in matlab time (for plotting)

    time.start      = duration(1,:);
    time.duration   = (datenum(duration(2,:))-datenum(duration(1,:)))*24; % days * hours/day
    time.dt_Tr      = []; % solve transport function every X hours
    time.dt_Up      = []; % solve uptake every X hours
    time.dt_Gr      = 24; % solve growth every X hours
    time.dt_ROMS    = 24; % based on extracted ROMS files
    time.timevec_Gr = datenum(time.start):time.dt_Gr/24:datenum(time.start)+time.duration/24-1;
    time.timevec_ROMS = datenum(time.start):time.dt_ROMS/24:datenum(time.start)+time.duration/24-1;
    
end  
