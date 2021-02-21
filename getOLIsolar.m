function [solarZ,phi0,acquisitionDate]=getOLIsolar(ldir)
%get solarZ and phi data for L8 OLI

d=dir(fullfile(ldir,'*.xml'));
xDoc=xmlread(fullfile(d.folder,d.name));
globalmeta = xDoc.getElementsByTagName('global_metadata');
list=globalmeta.item(0);
solarZ=NaN;
phi0=NaN;
for k=0:list.getLength-1
    thisListItem=list.item(k);
    if strcmp(thisListItem.getNodeName,'solar_angles')
        Alist=thisListItem.getAttributes;
        solarZA=Alist.getNamedItem("zenith");
        solarZ=str2double(solarZA.getFirstChild.getData.toString);
        phiA=Alist.getNamedItem("azimuth");
        phi0=str2double(phiA.getFirstChild.getData.toString);
        
    elseif strcmp(thisListItem.getNodeName,'acquisition_date')
        acquisitionDate=thisListItem.getFirstChild.getData.toString;
        acquisitionDate=datenum(char(acquisitionDate),'yyyy-mm-dd');
    end
    
end
% fid=fopen(fullfile(d.folder,d.name));
% c=textscan(fid,'%q %q','Delimiter','=');
% fclose(fid);
% n=strcmp('SUN_AZIMUTH ',c{1});
% phi0=str2double(c{2}{n});
% n=strcmp('SUN_ELEVATION ',c{1});
% solarZ=90-str2double(c{2}{n});

    
    
    