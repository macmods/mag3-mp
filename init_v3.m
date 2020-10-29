function [kelp_fr, kelp_ar] = init_v3(farm,time)
% Kelp, pre-allocation of MAG variables
% INPUT: farm, time
% OUTPUT:kelp_fr.()
%   xXX_yXX; data stored as a structure with single structure for each x
%   and y location on farm where kelp is outplanted
%     Ns; stored nitrogen (mg N/m frond)
%     Nf; fixed nitrogen (mg N/m frond)
%     Nf_capacity (mg N/m3)
%     Age [h]
%     ID; frond number, useful for tracking last frond
%     lastFrond; time that last frond was initiated, used by frond
%       initiation function
%     kelploc; [X,Y], used to translate from structure to farm grid
% OUTPUT: kelp_ar.()
%     variables tracked per area are the sum of all fronds within given
%     location
          

%% Initialize kelp state variables for all plant spacing
    
        for xx = 1:length(farm.gline_grid)
        for yy = 1:length(farm.bline_grid)

            for gg = 0:farm.gline_length-1
            
                kelpid = sprintf('x%d_y%d',[farm.gline_grid(xx) farm.bline_grid(yy)+gg]);
            
                kelp_fr.(kelpid) = configkelp_v4(farm,time);
                kelp_fr.(kelpid).kelploc = [farm.gline_grid(xx) farm.bline_grid(yy)+gg];
                
            end
            
        end                 
        end
        
        clear xx yy gg kelpid
      
%% Initialize farm area
% Nf, rDON, rPON tracked per area

kelp_ar.Nf = NaN(size(farm.gridx));
kelp_ar.rDON = zeros(size(farm.gridx));
kelp_ar.rPON = zeros(size(farm.gridx));


end