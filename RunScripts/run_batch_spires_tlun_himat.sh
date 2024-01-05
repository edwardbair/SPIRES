#! /bin/bash
# set up desired release env
. /home/matlab/env/matlab-R2021b.bashrc
#set high priority for killing when out of memory
echo 1000 > /proc/$$/oom_score_adj

#machine specific parameters
#file an directories, note closing /
codedir='/home/nbair/code/SPIRES/'
mccmfile='/home/snowhydro/nbair/MccM/net.mat'
Ffile='/home/snowhydro/nbair/datasets/LUT/lut_modis_b1to7_3um_dust.mat'
hdfbasedir='/home/snowhydro/sandbox-snowhydro/MODIS/mod09ga/6.1/'


base='/home/snowhydro/nbair/datasets/SPIRESinputs/MODIS/'
R0dir=$base'R0/'
maskdir=$base'watermask/'
topodir=$base'Z/'
ccdir=$base'cc/'
ficedir=$base'fice/'

outbase='/home/snowhydro/nbair/datasets/SPIRES/'

cores=64

#parameters that shouldn't change, at least across machines
b_R=2.7
dust_rg_thresh=300
dust_thresh=0.90
el_cutoff=500
fsca_thresh=0.10
grain_thresh=0.30
maxdust=950
maxgrainradius=1190
mindust=0
mingrainradius=40
shade=0
tolval=0.05
windowSize=40 #45 2021 IEEE value
windowThresh=20 #13 2021 IEEE value
Nd=7

#for each water year
for WY in {2001..2021}
#for each tile
do
for tile in h22v04 h22v05 h23v04 h23v05 h23v06 h24v04 h24v05 h24v06 h25v04 h25v05 h25v06 h26v05 h26v06
do
R0file=$R0dir$tile'R0.mat'
maskfile=$maskdir$tile'watermask.mat'
topofile=$topodir$tile'Topography.h5'
outloc=$outbase$tile
ccfile=$ccdir'cc_21yr_medianyr_'$tile.'mat'
ficefile=$ficedir$tile'.mat'

matlab -batch "runSPIRESTileyear("\'$codedir\'","\'$R0file\'","\'$maskfile\'","\'$topofile\'",\
"\'$outloc\'","\'$Ffile\'","\'$hdfbasedir\'","\'$tile\'","\'$ccfile\'","\'$ficefile\'",\
"\'$mccmfile\'",$cores,$WY,$b_R,$dust_rg_thresh,$dust_thresh,$el_cutoff,$fsca_thresh,$grain_thresh,\
$maxdust,$maxgrainradius,$mindust,$mingrainradius,$shade,$tolval,$windowSize,$windowThresh,$Nd)"

done
done
