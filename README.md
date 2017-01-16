[![Travis-CI Build Status](https://travis-ci.org/pmarchand1/waver.svg?branch=master)](https://travis-ci.org/pmarchand1/waver)

## Purpose ##

The `waver` package was developed to provide R code for the calculation of fetch length - the open water distance over which wind can blow along a given direction - and predict wave energy based on the calculated fetch and user-supplied wind and wave monitoring data.

The  `fetch_len` and `fetch_len_multi` functions serve to calculate fetch at single or multiple locations, respectively. These functions accept vector spatial data types (SpatialPoints, SpatialLines and SpatialPolygons) in either projected or unprojected (long/lat) coordinates. Users can specify the bearings for which fetch should be calculated, the number of radii to draw at each bearing (and their angular separation), as well as the maximum fetch length to consider.

These functions are based on existing algorithms published by the USGS[1] and Natural Capital Project[2]. These tools are implemented in ArcGIS with Python code and use raster land layers, whereas our fetch calculation uses vector spatial data.

The `wave_energy` function calculates the wave energy flux from either 
(a) the significant wave height and peak period, as provided by wave models
such as NOAA's [WAVEWATCH III](http://polar.ncep.noaa.gov/waves/index2.shtml); 
or (b) the wind speed, fetch and water depth. The latter calculation is based
on the method presented in the US Army Corps of Engineers Coastal Engineering 
Manual[3].


## References ##

[1] J. Rohweder et al. (2008) "Application of Wind Fetch and Wave Models for Habitat Rehabilitation and Enhancement Projects." USGS report. http://www.umesc.usgs.gov/management/dss/wind_fetch_wave_models.html

[2] InVEST Coastal Vulnerability Model. http://data.naturalcapitalproject.org/invest-releases/documentation/2_2_0/coastal_vulnerability.html

[3] Resio, D.T., Bratos, S.M., and Thompson, E.F. (2003). Meteorology
 and Wave Climate, Chapter II-2. Coastal Engineering Manual.
 US Army Corps of Engineers, Washington DC, 72pp.


## Acknowledgements ##
 
This work was supported by the National Socio-Environmental Synthesis Center (SESYNC) under funding received from the National Science Foundation DBI-1052875.
