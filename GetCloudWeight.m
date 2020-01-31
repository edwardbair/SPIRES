function [ cloudMask, weight, dateval, varargout ] = GetCloudWeight( fileList)
% retrieve MOD09 cloud mask and confidence weights from set of cloud files
% [ cloudMask, weight, dateval [, hdr] ] = GetCloudWeight( h5files )
%
% Input
%   fileList - HDF 5 filename(s)
%
% Output
%   cloudMask - logical cloud mask of dimension rows x cols x
%       length(dateval)
%   weight - floating point (single) weights for each pixel, same size as
%       cloudMask
%   dateval - scalar or column vector of date values
%       in MATLAB datenum format
% Optional output
%   hdr - structure with coordinate information about the dataset

numOut = 3;
if ischar(fileList)
    fileList = {fileList};
end

% assemble
dateval = zeros(length(fileList),1);
for k=1:length(fileList)
    % check to make sure inputs are valid .h5 files
    h5file = fileList{k};
%     try
        % information from HDF file
        group = findMODISh5group( h5file,'500m' );
        info = h5info(h5file,group);
        siz = info.Datasets(1).ChunkSize;
        dataset = {[group '/cloudMask'], [group '/confidenceWeights']};
        
        % date of this file
        dateval(k) = h5readatt(h5file,'/','MATLABdate');
        
        % cloud mask
        C = h5read(h5file,dataset{1},[1 1],siz)>0;
        if k==1
            cloudMask = C;
        else
            cloudMask = cat(3,cloudMask,C);
        end
        
        % weights, converted to single
        W = h5read(h5file,dataset{2},[1 1],siz);
        divisor = h5readatt(h5file,dataset{2},'divisor');
        if k==1
            weight = single(W)/divisor;
        else
            weight = cat(3,weight,single(W)/divisor);
        end
        
        % other output argument?
        if nargout>numOut && k==1
            hdr = GetCoordinateInfo(h5file,group,size(weight));
            varargout{1} = hdr;
        end
        
%     catch
%         warning('%s does not appear to be a valid HDF 5 cloud file',h5file)
%         cloudMask = [];
%         weight = [];
%         dateval = [];
%         if nargout>numOut
%             varargout{1} = [];
%         end
%         return
%     end
end