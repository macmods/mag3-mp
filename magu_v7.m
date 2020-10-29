function [magu,ux,vy,wz] = magu_v7(envt,farm,envt_step,arg_in1,arg_in2)
% Calculate magnitude velocity across entire farm, used by uptake
% Can either be ROMS informed (homogenous fields) or LES informed
% (heterogeneous fields)
if strcmp(arg_in1,'ROMS')
    
    % Generate a homogenous field (*not LES informed*)
    % delta x needs to match dx for growth (full farm dimensions

        ux = permute(repmat(envt.u(:,envt_step),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
        vy = permute(repmat(envt.v(:,envt_step),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
       
            
elseif strcmp(arg_in1,'LES')
    
    % Generate a homogenous field (*not LES informed*)
    % delta x needs to match dx for growth (full farm dimensions

        ux = permute(repmat(envt.u(:,envt_step),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
        vy = permute(repmat(envt.v(:,envt_step),1,size(farm.gridx,1),size(farm.gridx,2)),[2 3 1]);
        
        
    % Now modify this field with LES based on kelp status
    
        if strcmp(arg_in2,'case_0')
            
            % do not modify ROMS flow with LES velocity deficits
            % create a dummy u_def
            u_def = zeros(size(envt.u_def_a));
            v_def = zeros(size(envt.u_def_a));

        elseif strcmp(arg_in2,'case_a')
            
            %modify ROMS flow with LES "subsurface" velocity deficits
            if envt.mld(envt_step) >= farm.z_cult
            u_def = envt.u_def_a;
            v_def = envt.v_def_a;
            
            elseif envt.mld(envt_step) < farm.z_cult
            u_def = envt.u_def_b;
            v_def = envt.v_def_b;
            end
            
        elseif strcmp(arg_in2,'case_b')
            
            if envt.mld(envt_step) >= farm.z_cult
            %modify ROMS flow with LES "max canopy" velocity deficits
            u_def = envt.u_def_c;
            v_def = envt.v_def_c;
             
            elseif envt.mld(envt_step) < farm.z_cult
            u_def = envt.u_def_d;
            v_def = envt.v_def_d;
            end
                
        end

        u_def = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,u_def,farm.gridx,farm.gridy,farm.gridz);
        v_def = interpn(farm.gridx_tr,farm.gridy_tr,farm.gridz_tr,v_def,farm.gridx,farm.gridy,farm.gridz);

        ux = u_def .* ux + ux;
        vy = v_def .* vy + vy;

end
            
% magnitude velocity -> to be used in uptake; needs to be on dx_gr
magu = sqrt(ux.^2 + vy.^2);

end
