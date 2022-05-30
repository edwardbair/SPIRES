function [X,mdates,rrb]=reprojectSPIREScube(fname,endmember,x11,y11,x22,...
    y22,px_sz,target_p)
%reproject SPIRES endmember
%fname - h5 filename, e.g. 2021.h5';
%endmember- endmember name, 'snow','grain', or 'dust'
%x11,y11,x22,y22 - bounding box to crop w/ 11 as top (nw) left corner
%px_sz - target px size in target_p units
%target_p - proj CRS target

[x,mdates,hdr]=GetEndmember(fname,endmember);

minX=min(x(x(:)>0));

[X,rrb]=rasterReprojection(x,hdr.RasterReference,...
    'InProj',hdr.ProjectionStructure,'OutProj',...
    target_p,'pixelsize',[px_sz px_sz],....
     'XLimit',[x11 x22],'YLimit',[y22 y11]);
X(X<minX)=0;

end