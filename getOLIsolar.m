function [solarZ,phi0]=getOLIsolar(ldir)
%get solarZ and phi data for L8 OLI

d=dir(fullfile(ldir,'*MTL.txt'));
fid=fopen(fullfile(d.folder,d.name));
c=textscan(fid,'%q %q','Delimiter','=');
fclose(fid);
n=strcmp('SUN_AZIMUTH ',c{1});
phi0=str2double(c{2}{n});
n=strcmp('SUN_ELEVATION ',c{1});
solarZ=90-str2double(c{2}{n});