# SPIRES

Spectral unmixing with SCAGD placeholder name. Does two sensors right now. Tested w/ MATLAB 2019b.

1) Clone this repository and add the base directory and all of its sub-directories to your MATLAB path.

2) Download example data: ftp://ftp.snow.ucsb.edu/pub/org/snow/users/nbair/SCAGDexamples/
Some of the files, like the p42r34 Landsat 8 topography file are very large, 32GB.

3) For MODIS, the m file to run is smooth_and_run_modis.m
See the comments and ftp://ftp.snow.ucsb.edu/pub/org/snow/users/nbair/SCAGDexamples/modis/modis_april_2011_example.mat for all of the included inputs. Adjust input paths as necessary.

4) For LandSat 8 OLI, the m file to run is run_scagd_landsat.m
See the comments and ftp://ftp.snow.ucsb.edu/pub/org/snow/users/nbair/SCAGDexamples/landsat8/ for all of the included inputs. Adjust inputs path as necessary.

Papers on methods and validation are forthcoming

NB 2020-1-31
