function [solarZ,phi0,acquisitionDate]=getOLIsolar(ldir)
%get solarZ and phi data for L8 OLI
%works for collection 1,2, and HLS(S30)
%input:ldir - directory for reflectance files

d=dir(fullfile(ldir,'*.xml'));
xDoc=xmlread(fullfile(d.folder,d.name));
globalmeta = xDoc.getElementsByTagName('global_metadata');
list=globalmeta.item(0);
product='collection1';
if isempty(list) %collection2
    product='collection2';
    image_attributes = xDoc.getElementsByTagName('IMAGE_ATTRIBUTES');
    list=image_attributes.item(0);
end
if isempty(list) %HLS
   product='HLS';
   list = xDoc.getElementsByTagName('AdditionalAttribute');
end

solarZ=NaN;
phi0=NaN;

switch product
    case 'collection2'
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
    case 'collection1'
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
    case 'HLS'
        for k=0:list.getLength-1
            thisListItem=list.item(k);
            c1=thisListItem.item(1); %name field
            fieldname=char(c1.getFirstChild.getData.toString);
            c2=thisListItem.item(3); 
            c3=c2.getChildNodes;%value field
            x=c3.item(1).getFirstChild.getData.toString;
            switch fieldname
                case 'MEAN_SUN_ZENITH_ANGLE'
                    solarZ=str2double(x);
                case 'MEAN_SUN_AZIMUTH_ANGLE'
                    phi0=str2double(x);            
                case'SENSING_TIME'
                    acquisitionDate=datenum(char(x),'yyyy-mm-dd');
            end
        end
    end
end




