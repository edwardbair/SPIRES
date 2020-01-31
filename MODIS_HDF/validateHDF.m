function I = validateHDF( file )
% I = validateHDF( file )
%validate whether HDF file is okay
% return information structure if okay, [] if not

try
    I = hdfinfo(file,'eos');
catch
    warning('file %s not valid HDF file',file)
    I = [];
end

end