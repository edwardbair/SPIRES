# Snow Property Inversion from Remote Sensing (SPIReS)

## Landsat 8 OLI

1. To run an example scence from 20160426 for p42r34, download the zip and m file (2.2 GB) from:
https://snow.ucsb.edu/products/SPIRES/Landsat8/example/

2. Checkout the code, https://github.com/edwardbair/SPIRES.git

3. Add the code directory and all its subdirectories to your MATLAB path, "addpath(genpath([location where you checked out the code to]))"

4. run L8_spires_example.m from MATLAB. Tested using R2022B. This latest version of MATLAB will produce some warnings about pixcenters that can be ignored.

## MODIS
Requires MATLAB + Parallel Computing Toolbox + mapping & some other toolboxes

1. Download R0, cc, watermask, fice, and Z for a given tile or set of tiles from:
https://snow.ucsb.edu/products/SPIRES/MODIS/Inputs/MODIS/

For example, for h09v04, the R0 file is https://snow.ucsb.edu/products/SPIRES/MODIS/Inputs/MODIS/R0/h09v04R0.mat

3. Checkout the code, https://github.com/edwardbair/SPIRES.git

4. Adjust "RunScripts/run_batch_spires_example.sh"
   l9 - codedir - where you checout code to
   l10 - mccmfile - where the .net file in lives, i.e., https://github.com/edwardbair/SPIRES/blob/master/MccM/net.mat
   l11 - Ffile - lookup table, i.e., https://snow.ucsb.edu/products/SPIRES/MODIS/Sierra/ExampleData/lut_modis_b1to7_3um_dust.mat
   l12 - HDF MOD09GA reflectance files, i.e. https://snow.ucsb.edu/products/SPIRES/MODIS/Sierra/ExampleData/mod09ga/
   l15-20 - paths to inputs from 1.
   l 24 - # of cores for parpool
   l 27-41 - keep as is

   l 44 adjust WYs if needed
   l 47 adjust tiles if needed

9. Execute run_batch_spires_example.sh at the terminal or run as a SLURM job https://github.com/edwardbair/SPIRES/blob/master/RunScripts/SPIRES_tlun_slurm.sh

Notes:

Running the full year will take a long time and a lot of RAM, depending on the number of cores used. Using 60 AMD EPYC cores, plan on about 6 hours and about 400GB RAM, or about 7 GB RAM/core. For testing, you can run with a single core.

The minimum amount of time that'll work for smoothSPIREScube is 1 calendar month, i.e dom 1 through 28 to 31, but the smoothing needs a full water year. You can increase "tolval" to speed up computations and decrease quality. You're smart, you'll figure it out.

Reference:

Bair, E.H., Stillinger, T., and Dozier, J. (2021) Snow Property Inversion from Remote Sensing (SPIReS), IEEE Transactions on Remote Sensing and Geoscience, doi: 10.1109/TGRS.2020.3040328

NB 2024-01-05
