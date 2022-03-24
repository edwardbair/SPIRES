function [X,mdates,rrb]=reprojectSPIREScube(fname,endmember,x11,y11,x22,...
    y22,px_sz,target_p)
%reproject SPIRES endmember
%fname - h5 filename, e.g. 2021.h5';
%endmember- endmember name, 'snow','grain', or 'dust'
%x11,y11,x22,y22 - bounding box to crop w/ 11 as top (nw) left corner
%px_sz - target px size in target_p units
%target_p - proj CRS target

[x,mdates,hdr]=GetEndmember(fname,endmember);

% %target proj, US Albers
% target_p=projcrs(102003,'authority','esri');
%bounding box
% x11=-2442000;
% x22=130000;
% y11=1774000;
% y22=-558000;
%reproject cube

minX=min(x(x(:)>0));

[X,rrb]=rasterReprojection(x,hdr.RasterReference,...
    'InProj',hdr.ProjectionStructure,'OutProj',...
    target_p,'pixelsize',[px_sz px_sz],....
     'XLimit',[x11 x22],'YLimit',[y22 y11]);
X(X<minX)=0;

end

%plot
% i=100; % day of water year 100
% st=datestr(mdates(i)); %corresponding date
% imagesc(x2(:,:,i));axis image;colorbar; %plot
% title(st);