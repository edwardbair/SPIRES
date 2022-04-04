function R=getS2sr(sdir,res,varargin)
%retrieve Sentinel2 surface refl
%input: sdir - directory of SR jp2 files, IMG_DATA, e.g.
% ../S2A MSIL2A 20201028T184511 N0214 R070 T11SLB 20201028T210424/
% S2A_MSIL2A_20201028T184511_N0214_R070_T11SLB_20201028T210424.SAFE/
% GRANULE/L2A_T11SLB_A027952_20201028T184911/IMG_DATA/

% res - string for resolution needed, eg '20' for 20m
% subset pixels, [row 1 row2;col1 col2], optional
%output R - struct with fields bands & RR, a map cell ref object

d=dir(fullfile(sdir,['R' res 'm'],'*_B*.jp2'));

if isempty(d)
    error('cannot find surface refl files');
end
cropflag=false;
if nargin==3 %crop
    if ~isempty(varargin{1})
        c=varargin{1};
        cropflag=true;
    end
end
for i=1:length(d)
    fname=fullfile(d(i).folder,d(i).name);
    X=single(imread(fname));
    X(X==0)=NaN;
    X=X*1e-4; %rescale
    if i==1
        %read metadata file 2 directory up
        s=strsplit(d(i).folder,filesep);
        s=s(1:end-2);
        xmlfname=fullfile(strjoin(s,filesep),'MTD_TL.xml');
        T=readtable(xmlfname);
        %xml2struct fields stopped being compatible w/ below starting in
        %R2021b
%         S=xml2struct(xmlfname);
        % geolocation
%         gp=S.n1:Level_dash_2A_Tile_ID.n1_colon_Geometric_Info.Tile_Geocoding.Geoposition;
%         for ii=1:length(gp)
%             if strcmp(gp{ii}.Attributes.resolution,res)
%                 idx=ii;
%             end
%         end
%         ulx=str2double(gp{idx}.ULX.Text);
%         uly=str2double(gp{idx}.ULY.Text);
%         dx=str2double(gp{idx}.XDIM.Text);
%         dy=str2double(gp{idx}.YDIM.Text);
        idx=T.resolutionAttribute==str2double(res);
        ulx=T.ULX(idx);
        uly=T.ULY(idx);
        dx=T.XDIM(idx);
        dy=T.YDIM(idx);

%         epsg=S.n1_colon_Level_dash_2A_Tile_ID.n1_colon_Geometric_Info.Tile_Geocoding.HORIZONTAL_CS_CODE.Text;
%         epsgfid = fopen('temp.txt','rt');
fid=fopen(xmlfname);
CC = textscan(fid,'%s');
CC=CC{1};
cc=regexp(CC,'(?<=\<HORIZONTAL_CS_CODE>EPSG:).[0-9]+','match');
tt=cellfun(@isempty,cc);
epsg=str2double(cc{~tt});
cc=regexp(CC,'<Mean_Sun_Angle>','match');
tt=cellfun(@isempty,cc);
idtt=find(~tt);
R.solarzenith=str2double(regexp(CC{idtt+2},...
    '(?<=unit="deg">).*(?=<\/ZENITH_ANGLE)','match'));
R.solarzenith=str2double(regexp(CC{idtt+4},...
    '(?<=unit="deg">).*(?=<\/AZIMUTH_ANGLE)','match'));
%solar data
%         SA=S.n1_colon_Level_dash_2A_Tile_ID.n1_colon_Geometric_Info.Tile_Angles.Mean_Sun_Angle;
%         R.solarzenith=str2double(SA.ZENITH_ANGLE.Text);
%         R.solarazimuth=str2double(SA.AZIMUTH_ANGLE.Text);
fclose(fid);

        sz=size(X);
            xlims=[ulx ulx+sz(2)*dx];
            ylims=[uly+sz(1)*dy uly];
        RR=maprefcells(xlims,ylims,size(X),...
            'ColumnsStartFrom','north');
        if cropflag
            [x,y]=worldGrid(RR,'fullgrid');
            x=x(c(1,1):c(1,2),c(2,1):c(2,2));
            y=y(c(1,1):c(1,2),c(2,1):c(2,2));
            xlims=[x(1,1)-dx/2 x(1,end)+dx/2];
            ylims=[y(end,1)-abs(dy/2) y(1,1)+abs(dy/2)]; %dy is negative
            RR=maprefcells(xlims,ylims,size(x),...
            'ColumnsStartFrom','north');
        end
        RR.ProjectedCRS=projcrs(epsg);
        R.RR=RR;
        
    end
    if cropflag
        X=X(c(1,1):c(1,2),c(2,1):c(2,2));
    end
    R.bands(:,:,i)=X;
end
end