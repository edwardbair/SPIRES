function [mstruct, RefMat500, RefMat1000, RasRef500, RasRef1000] =...
    sinusoidProjH5(H5filename)
% [mstruct, RefMat500, RefMat1000, RasRef500, RasRef1000] = sinusoidProjH5(H5filename)
%INPUT
% filename - HDF 5 filename
%
%OUTPUT
% mstruct - projection structure
% RefMat500, RefMat1000 - referencing matrices for 500m and 1km data
%
% sinusoidProjH5 still works, but has been replaced by sinusoidProjMODtile

% projection structure
mstruct = defaultm('sinusoid');

% info from file
S = hdfinfo(H5filename,'eos');

% referencing matrices - set to null in case of lack of data in file
RefMat1000 = [];
RefMat500 = [];
for k=1:2
    assert(strcmp(S.Grid(k).OriginCode,'ul'),...
        'grid %s does not start in upper left, need to modify the code',...
        H5filename)
    assert(strcmpi(S.Grid(k).Projection.ProjCode,'snsoid'),...
        'grid %s is not in a sinusoidal projection',H5filename);
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
            H5filename)
        disp(S)
    end
end

mstruct.geoid=S.Grid(1).Projection.ProjParam(1:2);
mstruct = defaultm(mstruct);

end