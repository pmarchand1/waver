waver 0.3.0
==================

* Code updated to use spatial geometry functions from the sf package to replace
imports from the retired rgdal and rgeos packages.

* The package now internally uses sf objects instead of the sp (SpatialPoints/Lines/Polygons)
object format, but sp objects are still accepted as inputs for backwards compatibility.

waver 0.2.1
==================

* Fixed a bug in the default arguments to `fetch_len_multi`


waver 0.2.0
==================

* First release on CRAN
