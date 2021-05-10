function [solarZ,phi0,acquisitionDate]=getOLIsolar(ldir)
%get solarZ and phi data for L8 OLI

d=dir(fullfile(ldir,'*.xml'));
xDoc=xmlread(fullfile(d.folder,d.name));
globalmeta = xDoc.getElementsByTagName('global_metadata');
list=globalmeta.item(0);
c2flag=false;
if isempty(list) %
    c2flag=true; %collection 2
    image_attributes = xDoc.getElementsByTagName('IMAGE_ATTRIBUTES');
    list=image_attributes.item(0);
end
solarZ=NaN;
phi0=NaN;

if c2flag
    for k=0:list.getLength-1
        thisListItem=list.item(k);
        if strcmp(thisListItem.getNodeName,'SUN_ELEVATION')
            solarZ=90-str2double(thisListItem.getFirstChild.getData.toString);
        elseif strcmp(thisListItem.getNodeName,'SUN_AZIMUTH')
            phi0=str2double(thisListItem.getFirstChild.getData.toString);            
        elseif strcmp(thisListItem.getNodeName,'DATE_ACQUIRED')
            acquisitionDate=thisListItem.getFirstChild.getData.toString;
            acquisitionDate=datenum(char(acquisitionDate),'yyyy-mm-dd');
        end
    end
else
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
end




