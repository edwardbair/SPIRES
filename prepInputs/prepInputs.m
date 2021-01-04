function out=prepInputs(pshapefile,tiles,demdir,ccdir,glimsfile,waterdir)
%create inputs cropped to basinshape
%inputs:
% shapefile - location of basin shapefile, lat,lon, & geographic projection
% tiles - cell mimimicking tiling layout
% demdir - location where h5 file Z values that can be retrieved using
% GetElevation
% ccdir - location of MOD44B files
% glimsfile - GLIMS shapefile location
% waterdir - location of MOD44W files

%read polyshape, P
load(pshapefile)

vars={'Z','cc','mask'};
out_tile_R=zeros(3,2,length(vars));
allocated=false;

for i=1:length(vars)
    for j=1:length(tiles)
        if i==1 %Z, elevation
            d=dir(fullfile(demdir,['*' tiles{j} '*.h5']));
            [out_tile,out_tile_hdr]=GetTopography(fullfile(demdir,d.name),...
                'elevation');
            out_tile_R(:,:,j)=out_tile_hdr.RefMatrix;
        elseif i==2 %cc
            d=dir(fullfile(ccdir,['*' tiles{j} '*.hdf']));
            out_tile=hdfread(fullfile(ccdir,d(end).name),'Percent_Tree_Cover');
            out_tile=single(out_tile);
            out_tile(out_tile>100)=0;
            out_tile=out_tile./100;
            out_tile=imresize(out_tile,0.5); %resize to 463m
        elseif i==3 %mask
            d=dir(fullfile(waterdir,['*' tiles{j} '*.hdf']));
            out_tile=hdfread(fullfile(waterdir,d(end).name),'water_mask');
            out_tile(out_tile<0 | out_tile>1)=0;
            out_tile=logical(out_tile);
            out_tile=imresize(out_tile,0.5,'nearest'); %resize to 463 m
        end
       
        %allocate mosaic
        if ~allocated
            %full size mosaic
            outfull=zeros(size(out_tile,1)*size(tiles,1),...
                size(out_tile,2)*size(tiles,2));
            outfull_RefMatrix=out_tile_R(:,:,1);%upper left
            
            [lon,lat]=P.boundary;

             [x,y]=projfwd(out_tile_hdr.ProjectionStructure,...
                 lat,lon);

            [r,c]=map2pix(outfull_RefMatrix,x,y);
            
            maskfull=poly2mask(c,r,size(outfull,1),size(outfull,2));
            
            [x,y]=pixcenters(outfull_RefMatrix,size(maskfull));
            [r,~]=find(maskfull);
            sub_r=min(r):max(r);
            [c,~]=find(maskfull');
            sub_c=min(c):max(c);
            basinmask=maskfull(sub_r,sub_c);
            x11=x(min(c));
            y11=y(min(r));
            %assume uniform spacing across mosaic
            dx=outfull_RefMatrix(2,1);
            dy=outfull_RefMatrix(1,2);
            out.hdr.RefMatrix=makerefmat(x11,y11,dx,dy);
            %assume uniform projection structure
            out.hdr.ProjectionStructure=out_tile_hdr.ProjectionStructure;
            out.hdr.RasterReference=refmatToMapRasterReference(...
                out.hdr.RefMatrix,size(basinmask));
            allocated=true;
        end
        %fill mosaic
        %get map coordinates for tile
        [x,y]=pixcenters(out_tile_R(:,:,j),size(out_tile));
        %convert to full size mosaic coordinates
        [r,c]=map2pix(outfull_RefMatrix,x,y);
        r=round(r);
        c=round(c);
        %fill full size mosaic
        outfull(r,c)=out_tile;
    end
    %subset full mosaic to basin
    out.(vars{i})=outfull(sub_r,sub_c);
    if i==3 %masked areas are not in basinmask or water
        out.(vars{i})= ~basinmask | out.(vars{i});
    end
end

% treat fice last
S=shaperead(glimsfile);
t1=tic;
fprintf('making ice mask\n')
out.fice=makeIceMask(S,out.hdr,5);
t2=toc(t1);
fprintf('made ice mask in %g min\n',t2/60);

