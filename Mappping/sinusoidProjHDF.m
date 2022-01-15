function [mstruct, RefMat500, RefMat1000, RasRef500, RasRef1000] =...
    sinusoidProjHDF(HDFfilename)
% [mstruct, RefMat500, RefMat1000, RasRef500, RasRef1000] = sinusoidProjHDF(HDFfilename)
%INPUT
% filename - EOS HDF 4 filename
%
%OUTPUT
% mstruct - projection structure
% RefMat500, RefMat1000 - referencing matrices for 500m and 1km data
% RasRef500, RasRef1000 - map raster reference structures
%
% sinsusoidProjHDF still works, but for MODIS files has been replaced by sinusoidProjMODtile

% projection structure
mstruct = defaultm('sinusoid');

% info from file
S = hdfinfo(HDFfilename,'eos');

% referencing matrices - set to null in case of lack of data in file
RefMat1000 = [];
RefMat500 = [];
for k=1:length(S.Grid)
    assert(strcmp(S.Grid(k).OriginCode,'ul'),...
        'grid %s does not start in upper left, need to modify the code',...
        HDFfilename)
    assert(strcmpi(S.Grid(k).Projection.ProjCode,'snsoid'),...
        'grid %s is not in a sinusoidal projection',HDFfilename);
    dx=(S.Grid(k).LowerRight(1)-S.Grid(k).UpperLeft(1))/(S.Grid(k).Columns);
    dy=(S.Grid(k).LowerRight(2)-S.Grid(k).UpperLeft(2))/(S.Grid(k).Rows);
    x11 = S.Grid(k).UpperLeft(1)+dx/2;
    y11 = S.Grid(k).UpperLeft(2)+dy/2;
    if strfind(S.Grid(k).Name,'1km')
        RefMat1000 = makerefmat(x11,y11,dx,dy);
        RasRef1000 = refmatToMapRasterReference(RefMat1000,...
            [S.Grid(k).Rows S.Grid(k).Columns]);
    elseif strfind(S.Grid(k).Name,'500m')
        RefMat500 = makerefmat(x11,y11,dx,dy);
        RasRef500 = refmatToMapRasterReference(RefMat500,...
            [S.Grid(k).Rows S.Grid(k).Columns]);
    else
        warning('file %s does not have the right Grid information',...
            HDFfilename)
        disp(S)
    end
end

mstruct.geoid=S.Grid(1).Projection.ProjParam(1:2);
mstruct = defaultm(mstruct);

end