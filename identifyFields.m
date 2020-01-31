function [ endmember,attribute,varargout ] = identifyFields( S )
% identify the fields, and optionally the weights, cloud mask, and image mask
% in the structure

field = fieldnames(S);
e = 1;
a = 1;
% revised to just include the following
endmember = {'snow_fraction'};
attribute = {'grain_size','deltavis'};
for k=1:length(field)
    if ~isempty(strfind(field{k},'fraction')) && isempty(strfind(field{k},'snow_'))
        rmMember {e} = field{k}; %#ok<AGROW>
        e = e+1;
    elseif ~isempty(strfind(field{k},'drfs'))
        rmAttribute{a} = field{k}; %#ok<AGROW>
        a = a+1;
    end
end
rmThisField = cat(2,rmMember,rmAttribute);

if nargout>2
    % 1st optional argument are the weights
    if isfield(S,'weight')
        weight = S.weight;
    else
        weight = ones(size(S.(endmember{1})));
    end
    varargout{1} = weight;
    
    % 2nd optional output argument is the cloud mask
    if nargout>3
        if isfield(S,'cloudMask')
            cloudMask = S.cloudMask;
        else
            cloudMask = false(size(S.(endmember{1})));
        end
        varargout{2} = cloudMask;
    end
    
    % 3rd optional output argument is the image mask (only a portion of the
    % image is analyzed)
    if nargout>4
        if isfield(S,'imgMask')
            imgMask = S.imgMask;
        else
            imgMask = true(size(S.cloudMask,1),size(S.cloudMask,2));
        end
        varargout{3} = imgMask;
    end
    
    % rth optional argument is the cell list of fields to remove
    if nargout>5
        varargout{4} = rmThisField;
    end
end
end