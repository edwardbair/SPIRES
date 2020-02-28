source = '/raid/sandbox/snowhydro/scratch/ERA5/adaptor.mars.internal-1582298516.9798064-7131-35-5c472213-9953-4661-a5a5-b58870190dc3.nc';

finfo = ncinfo(source);

idx = strcmp({finfo.Variables.Name},'q');
sizeData = finfo.Variables(idx).Size;

%get lat lon fields
lon = ncread(source,'longitude');
lat = ncread(source,'latitude');

%grab data over CA.
latBB_CA = [32.5 42.5];
lonBB_CA = [-124.25 -114.5]+360;

latIdx = find(lat>=latBB_CA(1) & lat<latBB_CA(2));
lonIdx = find(lon>=lonBB_CA(1) & lon<lonBB_CA(2));
startIdx = [min(lonIdx) min(latIdx) 1 1];
countIdx = [max(lonIdx)-min(lonIdx) max(latIdx)-min(latIdx) inf inf];

varname ='q';
q = ncread(source,'q',startIdx,countIdx); 
u = ncread(source,'u',startIdx,countIdx); 
v = ncread(source,'v',startIdx,countIdx); 

