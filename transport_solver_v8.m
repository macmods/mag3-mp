function [NO3_field, NH4_field, DON_field, duration_dt] = transport_solver_v8(NO3_field,NH4_field,DON_field,Sink,Tr,kelp_ar,farm,time,duration)
% transport_solver calculates changes to nutrient field from sources/sinks
% and then runs it through the ADI_solver
%
% Input: Nutrient field
%        Tr contains constants for ADI_solver
%        farm, time
%        arg_in = nutrient, 'NO3' 'NH4' 'DON
%
% Output:Nutrient field (modified)


global param

%%
  
    % Convert the nutrient field from FARM GRID To TRANSPORT GRID
        
            NO3_field = interpn(farm.gridx,farm.gridy,farm.gridz,NO3_field,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
            NH4_field = interpn(farm.gridx,farm.gridy,farm.gridz,NH4_field,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
            DON_field = interpn(farm.gridx,farm.gridy,farm.gridz,DON_field,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
            
            NO3_fieldi = NO3_field; % store initial field
            NH4_fieldi = NH4_field;
            DON_fieldi = DON_field;
    
    
            % Exudation + Mortality = kelp_ar.rDON
            rDON_sm = NaN(size(kelp_ar.rDON));
            
            for zz = 1:farm.z % can stop at cultivation depth because no sink term beyond
                rDON_sm(:,:,zz) = smooth2a(kelp_ar.rDON(:,:,zz),farm.dx_tr,farm.dy_tr);
            end; clear zz
            
            rDON_tr = interpn(farm.gridx,farm.gridy,farm.gridz,rDON_sm,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
            clear rDON_sm

        

    % SINK term needs to be smoothed first and then inteprolated to
    % TRANSPORT GRID
    
        up_NO3 = Sink.NO3;
        up_NH4 = Sink.NH4;
        up_DON = Sink.DON;
            
        up_NO3_sm = NaN(size(up_NO3));
        up_NH4_sm = NaN(size(up_NH4));
        up_DON_sm = NaN(size(up_DON));
        for zz = 1:farm.z % can stop at cultivation depth because no sink term beyond
            up_NO3_sm(:,:,zz) = smooth2a(up_NO3(:,:,zz),farm.dx_tr,farm.dy_tr);
            up_NH4_sm(:,:,zz) = smooth2a(up_NH4(:,:,zz),farm.dx_tr,farm.dy_tr);
            up_DON_sm(:,:,zz) = smooth2a(up_DON(:,:,zz),farm.dx_tr,farm.dy_tr);
        end
        clear zz
        
        up_NO3_tr = interpn(farm.gridx,farm.gridy,farm.gridz,up_NO3_sm,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
        up_NH4_tr = interpn(farm.gridx,farm.gridy,farm.gridz,up_NH4_sm,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
        up_DON_tr = interpn(farm.gridx,farm.gridy,farm.gridz,up_DON_sm,farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,'linear');
        clear up_NO3_sm up_NH4_sm up_DON_sm
    

         
%% TRANSPORT
delta=0;
duration_dt = 0;
while (duration < 24 && delta < 0.2) ... % keep going if delta is less than 20%
    ||(duration < 24 && delta > 0.2 && duration_dt < 10/60) % keep going even if delta is > 20% but duration_dt is less than 10 minutes (e.g., not necessary to calc uptake more often than every hour)

            % STEP 1: APPLY SINK AND SOURCE TERMS
            % NO3 = NO3 [umol N/m3] - Uptake [umol N/m3/h] + Remin[DON]

            NO3_field = NO3_field; ...
                      %- up_NO3_tr .* time.dt_Tr ...
                      %+ DON_field .* 1e3 .* param.remin .* time.dt_Tr;
                      % SINK: uptake * biomass
                      % SOURCE: DON in mmol converted to umol * remin rate

            % STEP 1: APPLY SINK AND SOURCE TERMS
            % NH4 = NH4 [umol N/m3] - Uptake [umol N/m3/h]

            NH4_field = NH4_field; ...
                      %- up_NH4_tr .* time.dt_Tr;
                      % SINK: uptake * biomass
                      % SOURCE: none

            % STEP 1: APPLY SINK AND SOURCE TERMS
            % DON = DON [mmol N/m3] - Uptake [umol N/m3/h] -> [mmol
            % N/m3/h] - Remin[DON] + (Exudation+Mortality)[Ns][mmol
            % N/m3/h]
   
            DON_field = DON_field; ...
                      %- up_DON_tr ./ 1e3 .* time.dt_Tr ...
                      %- DON_field .* param.remin .* time.dt_Tr ...
                      %+ rDON_tr .* time.dt_Tr;
                  % SINK: uptake * biomass and remin to NO3
                  % SOURCE: Exudation and mortality froms Ns

            % STEP 2: TRANSPORT

            NO3_field = ADI_3D_v5(NO3_field,Tr.NO3);
            NH4_field = ADI_3D_v5(NH4_field,Tr.NH4);
            DON_field = ADI_3D_v5(DON_field,Tr.DON);
            
            % if delta nutrient is more than 10% recalc uptake
            delta = max(...
                [max(max(max(abs((NO3_field(:,:,1:farm.z_cult) - NO3_fieldi(:,:,1:farm.z_cult)) ./ NO3_fieldi(:,:,1:farm.z_cult)))))...
                 max(max(max(abs((NH4_field(:,:,1:farm.z_cult) - NH4_fieldi(:,:,1:farm.z_cult)) ./ NH4_fieldi(:,:,1:farm.z_cult)))))...
                 max(max(max(abs((DON_field(:,:,1:farm.z_cult) - DON_fieldi(:,:,1:farm.z_cult)) ./ DON_fieldi(:,:,1:farm.z_cult)))))]);
             
            duration_dt = duration_dt + time.dt_Tr;
            duration = duration + time.dt_Tr;
end

   
%% Convert from TRANSPORT GRID to KELP GRID and OUTPUT
% check for negative values and report are they due to transport function
% or interpolation
    
    % NO3
    if any(any(any(NO3_field < 0)))
    disp('NO3 negative transport function');
    end
        
    NO3_field = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,NO3_field,farm.gridx,farm.gridy,farm.gridz,'linear');

    if any(any(any(NO3_field < 0)))
    disp('NO3 negative interpolation -> replaced values');
    end

    NO3_field(NO3_field <= 0.001e3) = 0.001e3;

    % NH4
    if any(any(any(NH4_field < 0)))
    disp('NH4 negative transport function');
    end

    NH4_field = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,NH4_field,farm.gridx,farm.gridy,farm.gridz,'linear');

    if any(any(any(NH4_field < 0)))
    disp('NH4 negative interpolation -> replaced values');
    end

    NH4_field(NH4_field <= 0.001e3) = 0.001e3;

    % DON
    if any(any(any(DON_field < 0)))
    disp('DON negative transport function');
    end

    DON_field = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,DON_field,farm.gridx,farm.gridy,farm.gridz,'linear');

    if any(any(any(DON_field < 0)))
    disp('DON negative interpolation -> replaced values');
    end

    DON_field(DON_field <= 0.01) = 0.01;

            
end            