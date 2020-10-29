function kelp = configkelp_v4(farm,time)
global param

        preFrondx = farm.frondcount;
            
    kelp.Ns = NaN(preFrondx,farm.z_cult/farm.dz);
    kelp.Nf = NaN(preFrondx,farm.z_cult/farm.dz);
    kelp.Nf_capacity = NaN(preFrondx,1);
    kelp.Age = NaN(preFrondx,1);
    kelp.ID = NaN(preFrondx,1);
    
        clear preFrondx
    
    
%% Initial Kelp Conditions
% Nf, Ns, AGE, NF_CAPACITY, ID

    kelp.Nf_capacity(1,1) = param.Nf_capacity_subsurface;
    kelp.Nf(1,farm.z_cult) = param.Nf_capacity_subsurface; % equivalent to a single 1 m frond; [mg N]
    kelp.Ns(1,farm.z_cult) = ((20-param.Qmin)*param.Nf_capacity_subsurface)/param.Qmin; % corresponds to a Q of 20
    kelp.Age(1) = time.dt_Gr;
    kelp.ID(1) = 1;
    kelp.lastFrond = time.dt_Gr;
    
end