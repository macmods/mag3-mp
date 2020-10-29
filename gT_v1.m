function gT=gT_v1(Temp)
% Temperature-dependent growth; CONDITIONAL FORMULA
%
% Piece-wise formulation; conditional
% linear increase below gT1
% optimal between gT1 and gT2
% linear decrease above gT2 to gT3
% zero growth above gT3
% Adapted from Broch and Slagstad 2012; gT from Schiel and Foster 2015

global param

        gT = NaN(size(Temp));
        
        gT(Temp < param.Tmin) = 1/param.Tmin * Temp(Temp < param.Tmin);
        gT(Temp >= param.Tmin & Temp < param.Tmax) = 1;
        
            % Solve systems of equations where Tmax = 1; Tlim = 0;
            % Linear decrease from Tmax to Tlim with intercept b and slope m
            
            b = 1 / (-param.Tmax / param.Tlim + 1);
            m = 1 / (param.Tmax - param.Tlim);
            
        gT(Temp >= param.Tmax & Temp <= param.Tlim) = m .* Temp(Temp >= param.Tmax & Temp <= param.Tlim) + b;
        gT(Temp > param.Tlim) = 0;

        
end
