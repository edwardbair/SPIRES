function h5writeProjection(filename,location,mstruct)
% h5writeProjection(filename,location,mstruct)
% write attributes of a projection structure into HDF 5 file
%
% filename - HDF 5 file
% location - e.g. Group appropriate to this projection structure
% mstruct - MATLAB projection structure

% just the relevant fields of the input structure
fn = fieldnames(mstruct);
t = strfind(fn,'frame');
nk = length(t);
for k=1:length(t)
    if ~isempty(t{k})
        nk = k-1;
        break;
    end
end

% create new structure with just the relevant fields
S = struct;
for k=1:nk
    value = mstruct.(fn{k});
    if ~isempty(value)
        S.(fn{k}) = value;
    end
end

% write these to the HDF 5 file
fn = fieldnames(S);
for k=1:length(fn)
    value = S.(fn{k});
    % any vectors must be row vectors
    if iscolumn(value) && length(value)>1
        value = value';
    end
    h5writeatt(filename,location,fn{k},value)
end
end