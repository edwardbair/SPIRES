function [ variableList ] = MODSCAGvariables( filename )
%returns list of variables in a MODSCAG file
% [ variableList ] = MODSCAGvariables( filename )
%

if iscell(filename)
    filename = char(filename);
end

if ~isempty(strfind(filename,'.mat'))
    matObj = matfile(filename);
    variableList = who(matObj);
elseif ~isempty(strfind(filename,'.h5'))
    group = findMODISh5group(filename,'500m');
    info = h5info(filename,group);
    variableList = cell(1,length(info.Datasets));
    for k=1:length(info.Datasets)
        variableList{k} = info.Datasets(k).Name;
    end
else
    error('filename %s must be .mat or .h5',filename)
end

end