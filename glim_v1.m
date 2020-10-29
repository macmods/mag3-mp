function [Growth, gQ_fr, gT_fr, gE_fr, gH_fr] = glim_v1(Q,Type,Height_tot,envt,farm,envt_step,kelploc)
% Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH
%
% Input: (Q,Type,Height,envt,farm,envt_step,kelploc)
%        Growth only calculated for subsurface and canopy fronds
%
% Output: 
%   Growth, [h-1]
%   gQ, quota-limited growth 
%       from Wheeler and North 1980 Fig. 2   
%   gT, temperature-limited growth
%       piecewise approach taken from Broch and Slagstad 2012 (for sugar
%       kelp) and optimized for Macrocystis pyrifera
%       calls on gT_vX
%   gE, light-limited growth
%       from Dean and Jacobsen 1984
%       calls on gE_vX
%   gH, Sigmoidal rate decrease starting at 5% below max height
%       essential a mathematical function to have smoothed, slowed growth
%       rather than abrupt growth changes at max_height


global param

        
%% gQ

    % Fronds, either subsurface or canopy, don't calc for senescing
    frond = Type == 1 | Type == 2;

    gQ_fr = NaN(length(Q),1); 
    gQ_fr(frond) = (Q(frond) - param.Qmin) ./ (param.Qmax - param.Qmin); 
    gQ_fr(gQ_fr > 1) = 1;
    gQ_fr(gQ_fr < 0) = 0;
    
    
%% gT    

    gT = gT_v1(envt.T(1:farm.z_cult,envt_step));          
    gT_fr = repmat(gT,1,length(Q))';
    clear gT
   
    
%% gE
% light varies across the farm, so extract correct light field with kelploc
% Bertalanffy Growth Equation (Dean and Jacobsen 1984)
% Input: PAR [W/m2]

    PAR = squeeze(envt.PAR_field(kelploc(1),kelploc(2),1:farm.z_cult));
    
    % k sets curvature of relationship. qualitatively fit to match Dean
    % and Jacobsen 1984; 50% growth at PAR ~2.5 and near 100% growth at
    % PAR ~7+
        
        gE = 1-exp(param.kPAR*(PAR-param.PARc));
        
        % If values < 0 replace with zero. We are explicitely modeling
        % mortality and so growth shouldn't be negative.
        gE(gE < 0) = 0;
        
    gE_fr = repmat(gE,1,length(Q))';
    clear gE PAR
       
        
%% gH            

    gH_fr = 0.5 + 0.5 .* tanh(-(Height_tot - (param.Hmax-0.05*param.Hmax)));
    

%% Growth

    Growth = param.umax .* gQ_fr .* gT_fr .* gE_fr .* gH_fr;


end