function [Nf_new, Ns_new, Nf_capacity_new, Age_new, ID] = Si_v1(Q,ID,Nf_new,Ns_new,Nf_capacity_new,Age_new,farm)
% Making a new frond, finds the last frond and initiates at the next row
% A new frond is equivalent Nf of 1 m and enough Ns for Q to match average
% Q of rest of fronds

global param

%% What frond are we on?
next = min(find(isnan(ID)));

%% New frond
% Nf
Nf_new(next,farm.z_cult) = param.Nf_capacity_subsurface;

% Ns equivalent to average Q of fronds in area ... 
Ns_new(next,farm.z_cult) = ((nanmean(Q)-param.Qmin)*param.Nf_capacity_subsurface)/param.Qmin;

% Frond characteristics
Nf_capacity_new(next) = param.Nf_capacity_subsurface;
Age_new(next) = 0;
ID(next) = next;
           
end