function [] = create_netcdf_v1(kelp_fr,time,filename)


%% GENERATE NETCDF

    ncid = netcdf.create(filename,'NC_WRITE');
    
    % Define and Set Dimensions
    numloc = length(fieldnames(kelp_fr));
    numdim = 2; % x and y 
    numt = length(time.timevec_Gr);

    dim_numloc = netcdf.defDim(ncid,'Location',numloc);
    dim_numdim = netcdf.defDim(ncid,'X-by-Y',numdim);
    dim_numt = netcdf.defDim(ncid,'Time',numt);
    
    
%% Variables to be Saved
% Only need to save Biomass total and canopy
    
    % Identifiers: Location, Time, Depth
    netcdf.defVar(ncid,'Frond Location','NC_DOUBLE',[dim_numloc dim_numdim]); % Create NetCDF variable
    netcdf.defVar(ncid,'Time','NC_DOUBLE',dim_numt); % Create NetCDF variable
    
    % KELP
    netcdf.defVar(ncid,'Btot','NC_DOUBLE',[dim_numloc dim_numt]); % Create NetCDF variable
    netcdf.defVar(ncid,'Bcan','NC_DOUBLE',[dim_numloc dim_numt]); % Create NetCDF variable
    
    
    netcdf.endDef(ncid) % End NetCDF file define mode
    netcdf.close(ncid);
    
    
%% Location of FRONDS

    kelplist = fieldnames(kelp_fr);
    for iField = 1:numel(fieldnames(kelp_fr))

        kelpid = char(kelplist(iField));
        kelploc(iField,:) = kelp_fr.(kelpid).kelploc;

    end
    clear iField kelplist kelpid    


ncwrite(filename,'Frond Location',kelploc)
  

end
    
    