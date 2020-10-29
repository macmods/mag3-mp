function param = param_v1
% structure with the named parameters of the model

    param = struct(...
        'VmaxNO3', 752,... % [umol N/m2/h]   
        'VmaxNH4', 739,... % [umol N/m2/h]    
        'VmaxDON', 108,... % [umol N/m2/h]   
        'KsNO3', 10.2e3,... % [umol N/m3]    
        'KsNH4', 5.31e3,... % [umol N/m3]     
        'KsDON', 7.755e3,... % [umol N/m3]     
        'Qmax', 40,... % [mg N/g(dry)]         
        'Qmin', 10,... % [mg N/g(dry)]        
        'umax', 0.2/24,... % [h-1]   
        'Tmin',14,... % [Celcius]
        'Tmax',20,... % [Celcius]
        'Tlim',23,... % [Celcius]
        'kPAR', -1/3,...
        'PARc', 1.7864,... % [W/m2]
        'PAR_Ksw', 0.0384,... % [m-1]
        'PAR_Kchla', 0.0138,... % [m2/mg chl-a]
        'PAR_KNf', 0.0001,... % [m2/mg N]
        'Hmax', 30,... % [m]
        'Frond_init', [0.042/10/24 0.17/24],... % [h-1]
        'Age_max', 150*24,... % [h]
        'd_dissolved', 0.002/24,... % [h-1]
        'd_wave_m',0.000455,... % [h-1]
        'd_blade', 0.009/24,... % [h-1]
        'd_frond', 0.1/24,... % [h-1]
        'dry_wet',0.094,... % [g(dry)/g(wet)]
        'Nf_capacity_subsurface', 109.98,... % [mg N/m]
        'Nf_capacity_watercolumn',98.7,... % [mg N/m]
        'Nf_capacity_canopy', 243.46,... % [mg N/m]
        'Biomass_surfacearea_subsurface', 32.2/1e4,... % [m2/g(wet)]
        'Biomass_surfacearea_watercolumn', 10.9/1e4,... % [m2/g(wet)]
        'Biomass_surfacearea_canopy', 58.7/1e4,... % [m2/g(wet)]
        'Blade_stipe', [2.32 13.36 23.39],...
        'remin', 0.01/24,... % [h-1]
        'MW_N', 14.006720,... % [g/mol]
        'MW_NO3', 62.0049,... % [g/mol]
        'MW_NH4', 18.038508); % [g/mol]
       
end


%    VmaxNO3 = maximum uptake rate based on Gerard (1982), Haines and
%    Wheeler 1978; [umol NO3/m2/h] important to have these units as uptake
%    is calculated as amount per area per time
%
%    VmaxNH4 = maximum uptake rate based on Haines and Wheeler (1974);
%    [umol NH4/m2/h]
%
%    VmaxDON = maximum uptake rate of DON-urea based on Smith et al.
%    (2018); [umol N/m2/h]
%
%    KsNO3 = half-saturation constant for NO3; Gerard (1982), Haines and
%    Wheeler (1978), Brzezinski et al. (2015); [umol NO3/m3]
%
%    KsNH4 = half-saturation constant for NH4, Haines and Wheeler (1974);
%    [umol NH4/m3]
%
%    KsDON = half-saturation constant for DON-urea, no empirical values
%    available; [umol N/m3]
%
%    Qmax = Maximum N quota, Brezinski et al. (2013), Rassweiler et al.
%    (2018); [mg N/g(dry)]
%
%    Qmin = Minimum N quota, Brezinski et al. (2013), Rassweiler et al.
%    (2018); [mg N/g(dry)]
%
%    umax = Maximum growth rate, informed by Rassweiler et al. (2018), max
%    observed by Wheeler and North (1980); [h-1]
%
%    gT1, gT2, gT3, temperature growth parameters, Zimmerman and Kremer
%    (1986), Schiel and Foster (2015); [Celcius]
%
%    kPAR = -1/3; rate constant for light-dependent growth; Dean and
%    Jacobsen (1984)
%
%    PARc = compensating PAR level; Dean and Jacobsen (1984); [W/m2]
%
%    Max_height; SBC LTER; [n]
%
%    Frond_init_m, Frond_init_b, the initation rate of fronds, linear
%    relationship with Q where Initiaiton = m * Q + b; Rodriguez et al.
%    (2013), Bell et al. (2018), pers. comm. Tom Bell; [d-1]
%
%    Mortality_age = rate of age-dependent mortality [h-1]; Rodriguez et
%    al. (2013)
%
%    Mortality_wave = rate of wave-dependent mortality; based on linear
%    relationship extracted from Rodriguez et al. (2013) where
%    Mortality_wave = m * Hs; [h-1]
%
%    Mortality_frond = rate of frond loss following senescence; no data
%    available on rate of senescence but 3-4 weeks from pers. comm. Tom
%    Bell, optimized to match a couple weeks; [h-1]
%
%    Mortality_blade = rate of mortality of blades from Rassweiler et al.
%    (2018); [h-1]
%
%    Mortality_dissolved = rate of dissolved loss from fronds from
%    Rassweiler et al. (2018); [h-1]
%
%    dry-wet = conversion of biomass from g(dry) to g(wet); Rassweiler et
%    al. (2018). This value is consistent in space and time and in blades
%    versus stipes
%
%    Nf_capacity = Biomass capacity factor; maximum allowable biomass
%    within a frond per m3; type dependent (subsurface, watercolumn,
%    canopy), allometric conversions based on Rassweiler et al. (2018)
%
%    Biomass_surfacearea = biomass to surface-area conversion; type
%    dependent (subsurface, watercolumn, canopy), allometric conversions
%    based on Haines and Wheeler (1978) and Fram et al. (2008)
%
%    Blade_stipe = blade to stipe weight conversion; based on Nyman et al.
%    (1993); = 2.32 - 13.36 * Fractional frond height + 23.39 * Fractional
%    frond height ^ 2
%
%    PAR_Ksw = 0.0384; % diffuse attenuation coefficient for water;
%    Lorenzen (1972); [m-1]
%
%    PAR_Kchla = 0.0138; % diffuse attenuation coefficient per mg chla;
%    Lorenzen (1972); [m2/mg chl-a]
%
%    PAR_KNf = 0.0001; % N-specific shading, Hadley et al. (2015); [m2/mg
%    Nf]
%
%    MW = Molecular weights of nitrogen, nitrate, and ammonium
%
%% Other constants embedded in functions
% configkelp_v1.m
%      preFrondx = 100 * TIME.duration/TIME.dt/365; amount of memory to
%      preallocate
%
% uptake_vX.m
%      visc = kinematic viscosity; [m2/s converted to m2/h] Dm = molecular
%      diffusion coefficient; [m2/s converted to m2/h] *** Dm = T *
%      3.65e-11 + 0.72e-10; determined for 18C (Li and Gregory 1974);
%      n_length = 25; number of iterations of wave-driven currents;
%      selected based on number of iterations attaining 95% of maximum
%      value for n_length = 1000
%
% growth_vX.m
%      gH1,gH2,gHk exponential dropoff in growth as frond height reaches
%      maximum height