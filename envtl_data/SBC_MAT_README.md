Environmental data extracted from ROMS 1-km grid by D. Dauhajre

lat = 34.36833333
lon= -119.73333333

Each variable has its own mat-file for 1997-2005.  

The files each contain the following variable array structures defined by the following keys:
1) the variable name (e.g., ‘PAR’), 2D array of (time, depth) 
2) ‘ocean_time’ —> 1D array (time) of time in seconds since model initialization (=Jan 1, 1994)
3) ‘depth’ —> 1D array (depth) of the z-grid of depths (1 m —> 60 m in 1 meter bins).

For reference all of the relevant data is on tethys at: /data/project3/sharedfiles/SBC_Farm/
The directory structure in the above directory is also follows:
/Samp_Out_Project3/ —> zoomed in netcdf files containing the ROMS output in the SBC region of the 1-km simulation (here the data is on the ROMS sigma levels)
/zslice_nc/  —> zoomed in netcdf files containing the ROMS output in the SBC region interpolated to a vertical grid of depths (1—>60 m)
/MAT_FILES/ —> the individual files I am attaching here

Wave Data: from NDBC buoy station 46053
34.252 N 119.853 W
https://www.ndbc.noaa.gov/station_page.php?station=46053
downloaded text files were converted to a single mat file to be imported by MAG

Variables:
Wave period (Tw)
Sig. wave height (Hs)
