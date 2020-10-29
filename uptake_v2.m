function [UptakeN, UptakeFactor] = uptake_v2(kelp_fr,envt,farm,envt_step,arg_in)
% Determination of uptake rate for nitrate, ammonium, and urea.
%
%   Input: (ENVT,Q,Type,farm)
%           NO3 = nitrate concentration in seawater [umol NO3/m3]
%           NH4 = ammonium concentration in seawater [umol NH4/m3]
%           DON = dissolved organic nitrogen concentration [mmol/m3]
%           mag_u = seawater velocity [m/h]; magnitude velocity of x,y,z
%           Tw = wave period [h]
%           Type = 'subsurface/canopy' for conversion from surface area to
%           g(dry)
%           FARM_DIM = dimensions of farm needed for allometric
%           equations
%
%    Output:
%           Uptake, [mg N/g(dry)/h]; converted from [umol NO3+NH4+N/m2/h]
%           UptakeFactor.UptakeFactor_NO3/NH4/DON, dimensionless; 0-1 
%           UptakeFactor.Uptake_NOS/NH4/DON_mass, [umol NO3;NH4;N /g(dry)/h]
%
%    Note -> Uptake is only calculated for Canopy and Subsurface fronds. If
%    fronds are senscing, uptake not determined (NaN).


global param


%% KELP INPUT

kelploc = kelp_fr.kelploc;
Q = kelp_fr.Q;
Type = kelp_fr.Type;

    subsurface = Type==1;
    canopy     = Type==2;

    
%% ENVT INPUT

if strcmp(arg_in,'transportON') 
NO3   = squeeze(envt.NO3_field(kelploc(1),kelploc(2),1:farm.z_cult));
NH4   = squeeze(envt.NH4_field(kelploc(1),kelploc(2),1:farm.z_cult));
DON   = squeeze(envt.DON_field(kelploc(1),kelploc(2),1:farm.z_cult));
else
NO3   = envt.NO3(1:farm.z_cult,envt_step);
NH4   = envt.NH4(1:farm.z_cult,envt_step);
DON   = envt.DON(1:farm.z_cult,envt_step);
end

% mag_u = squeeze(envt.magu_field(kelploc(1),kelploc(2),1:farm.z_cult));
% Tw = envt.Tw(1,envt_step);


%% Q                                        
% Quota-limited uptake: maximum uptake when Q is minimum and
% approaches zero as Q increases towards maximum; Possible that Q
% is greater than Qmax. Set any negative values to zero.

    % By only calculating uptake for subsurface and canopy fronds
    % ensures that dead and senescing fronds don't uptake
    UptakeFactor.vQ             = NaN(length(Q),1);
    UptakeFactor.vQ(subsurface) = (param.Qmax-Q(subsurface))./(param.Qmax-param.Qmin);
    UptakeFactor.vQ(canopy)     = (param.Qmax-Q(canopy))./(param.Qmax-param.Qmin);

    % Uptake factor between 0 and 1
    % negative and positive values could happen if Q exceeds Qmin or
    % Qmax
    UptakeFactor.vQ(UptakeFactor.vQ < 0) = 0;
    UptakeFactor.vQ(UptakeFactor.vQ > 1) = 1;


%% Michaelis-Menten FACTOR: Kinetically limited uptake only, 
% ranges from 0-1
% umol/m3 / umol/m3 + umol/m3

        UptakeFactor.vMichaelis_NO3 = NO3 ./ (param.KsNO3+NO3); 
        UptakeFactor.vMichaelis_NH4 = NH4 ./ (param.KsNH4+NH4);
        UptakeFactor.vMichaelis_DON = DON.*0.2.*1e3 ./ (param.KsDON+DON);

        
%% Calc Michaelis-Menten uptake rate
% umol/m2/h

        Uptake_NO3 = param.VmaxNO3 .* UptakeFactor.vMichaelis_NO3;
        Uptake_NH4 = param.VmaxNH4 .* UptakeFactor.vMichaelis_NH4;
        Uptake_DON = param.VmaxDON .* UptakeFactor.vMichaelis_DON;
        
        
%% Convert from surface area to g(dry) and to mg N 
% Based on allometric conversions from param that are dependent on kelp
% type (subsurface, canopy, watercolumn) 
% [mg N/g(dry)/h] 

        % Pre-allocate space because calculating per frond type (only
        % partially filling in matrix at a time
        Uptake_NO3_mass = NaN(length(Type),farm.z_cult/farm.dz);
        Uptake_NH4_mass = NaN(length(Type),farm.z_cult/farm.dz);
        Uptake_DON_mass = NaN(length(Type),farm.z_cult/farm.dz);


        Uptake_NO3_mass(subsurface,:) = repmat(Uptake_NO3',nnz(subsurface),1) ... % repeat for number of subsurface fronds
                                     .* UptakeFactor.vQ(subsurface,1) ... % Q-limited uptake
                                     .* param.Biomass_surfacearea_subsurface*2 ./ param.dry_wet ... % multiply by two because surface are on both sides
                                     .* param.MW_N ./ 1e3; % convert to mg N
        Uptake_NO3_mass(canopy,1)     = repmat(Uptake_NO3(1,:)',nnz(canopy),1) ...
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_canopy*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;
        Uptake_NO3_mass(canopy,2:end) = repmat(Uptake_NO3(2:farm.z_cult,:)',nnz(canopy),1) ...
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_watercolumn*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;

        Uptake_NH4_mass(subsurface,:) = repmat(Uptake_NH4',nnz(subsurface),1) ...
                                     .* UptakeFactor.vQ(subsurface,1) ...
                                     .* param.Biomass_surfacearea_subsurface*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;
        Uptake_NH4_mass(canopy,1)     = repmat(Uptake_NH4(1,:)',nnz(canopy),1) ...
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_canopy*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;
        Uptake_NH4_mass(canopy,2:end) = repmat(Uptake_NH4(2:farm.z_cult,:)',nnz(canopy),1) ...
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_watercolumn*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;

        Uptake_DON_mass(subsurface,:) = repmat(Uptake_DON',nnz(subsurface),1) ...
                                     .* UptakeFactor.vQ(subsurface,1) ...
                                     .* param.Biomass_surfacearea_subsurface*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;
        Uptake_DON_mass(canopy,1)     = repmat(Uptake_DON(1,:)',nnz(canopy),1) ...
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_canopy*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;
        Uptake_DON_mass(canopy,2:end) = repmat(Uptake_DON(2:farm.z_cult,:)',nnz(canopy),1) ... 
                                     .* UptakeFactor.vQ(canopy,1) ...
                                     .* param.Biomass_surfacearea_watercolumn*2 ./ param.dry_wet ...
                                     .* param.MW_N ./ 1e3;


%% TOTAL Uptake = Uptake NO3 + Uptake NH4 + Uptake DON
% [mg N/g(dry)/h]

UptakeN = Uptake_NO3_mass + Uptake_NH4_mass + Uptake_DON_mass;
clear subsurface canopy

end