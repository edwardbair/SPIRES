function [x11,y11,dx,dy,inProj] = cornerCoords(R)
% corner coordinates from raster reference
if contains(class(R),'MapCellsReference') ||...
        contains(class(R),'MapPostingsReference')
    inProj = true;
    [x11,y11] = intrinsicToWorld(R,1,1);
    if contains(class(R),'Cells')
        if strcmp(R.ColumnsStartFrom,'north')
            dy = -R.CellExtentInWorldY;
        else
            dy = R.CellExtentInWorldY;
        end
        if strcmp(R.RowsStartFrom,'west')
            dx = R.CellExtentInWorldX;
        else
            dx = -R.CellExtentInWorldX;
        end
    else % Postings instead of Cells
        if strcmp(R.ColumnsStartFrom,'north')
            dy = -R.SampleSpacingInWorldY;
        else
            dy = R.SampleSpacingInWorldY;
        end
        if strcmp(R.RowsStartFrom,'west')
            dx = R.SampleSpacingInWorldX;
        else
            dx = -R.SampleSpacingInWorldX;
        end
    end
elseif contains(class(R),'GeographicCellsReference') ||...
        contains(class(R),'GeographicPostingsReference')
    inProj = false;
    [y11,x11] = intrinsicToGeographic(R,1,1);
    if contains(class(R),'Cells')
        if strcmp(R.ColumnsStartFrom,'north')
            dy = -R.CellExtentInLatitude;
        else
            dy = R.CellExtentInLatitude;
        end
        if strcmp(R.RowsStartFrom,'west')
            dx = R.CellExtentInLongitude;
        else
            dx = -R.CellExtentInLongitude;
        end
    else % Postings
        if strcmp(R.ColumnsStartFrom,'north')
            dy = -R.SampleSpacingInLatitude;
        else
            dy = R.SampleSpacingInLatitude;
        end
        if strcmp(R.RowsStartFrom,'west')
            dx = R.SampleSpacingInLongitude;
        else
            dx = -R.SampleSpacingInLongitude;
        end
    end
else
    error('raster reference class %s unrecognized',class(R));
end
end