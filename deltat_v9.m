function [time] = deltat_v9(scenario,envt,farm,time,ROMS_step)
% Sets time step for uptake and transport loops
% Dynamic to improve efficiency
% Time steps need to be shorter when applying LES fields with lots of
% biomass on farm; and time steps can be longer when there isn't much kelp
% and only using ROMS input for transport function


%% Dynamic transport time step
% Courant number
c = 1;      


    if strcmp(scenario,'case_0')
        
        % u/dx + v/dy + w/dz + Dh/dx*dy + Dz/dz2
        time.dt_Tr = min(c ./ ...
            (abs(envt.u(:,ROMS_step)./farm.dx_tr) ...
            + abs(envt.v(:,ROMS_step)./farm.dy_tr) ...
            + abs(envt.w(:,ROMS_step)./farm.dz_tr) ...
            + envt.Dh./(farm.dx_tr*farm.dy_tr) ...
            + envt.Dv(:,ROMS_step)./(farm.dz_tr*farm.dz_tr)));
  
    else
        
        time.dt_Tr = 30/60/60; % 30 seconds -> constant
        
    end

    if ROMS_step == 159 || 167
        
        time.dt_Tr = 10/60/60; % 10 seconds -> constant
        
    end
end