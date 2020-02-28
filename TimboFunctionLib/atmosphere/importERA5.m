%'C:\raid\scratch\ERA5\Dec2019_IVT_components\';
gribFN = 'C:\raid\scratch\ERA5\Dec2019_IVT_components\adaptor.mars.internal-1582226738.744073-24801-29-d029a0b8-e43e-4840-a0af-26cb3787c6a0.grib';

irec=-1;
grib_struct=read_grib(gribFN,irec,'HeaderFlag',true,'DataFlag',false','ScreenDiag',false);%,irec,varargin)

%grbfile = ncgeodataset(gribFN);