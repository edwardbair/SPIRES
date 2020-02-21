function [ X, varargout ] = GetEndmember( filename, whichVariable, varargin )
% X = GetEndmember( filename, whichVariable [, date] )
% [X,dates] = GetEndmember( filename, whichVariable [, date] )
% [ X, dates, hdr ] = GetEndmember( filename, whichVariable [, date] )
%
% retrieve MODSCAG endmember fraction or attribute (fractions for snow,
% vegetation, rock, shade, other; grain size in mm, drfs_grnsz in mm,
% deltavis as fraction)
%
% Input
%   filename - HDF 5 or .mat  filename
%   whichVariable - not all are in every file, but choices could be 'snow',
%   'rock' or 'soil', 'vegetation', 'grain_size', 'drfs_grnsz', 'deltavis',
%   'shade','dust'
%   (only first 4 letters are considered, and are case-insensitive)
% Optional input
%   date - followed by single date, vector of length 2 specifying range, or
%       vector of any length specifying individual dates to retrieve
%       (format for dates can be MATLAB datenum,  7-digit ISO date YYYYDDD,
%       or 8-digit ISOdate YYYYMMDD)
%       trick to retrieve just 2 separated dates: vector in reverse chronological
%       order
%
% Output
%   X - output matrix or cube of requested variable, converted to single-
%       precision floating point unless variable is a mask
% Optional output, in order
%   dates - vector of MATLAB dates retrieved
%   hdr - structure with coordinate information about the dataset

% parse input
p = inputParser;
defaultDaterange = [];
validationFcn = @(x) iscell(x) || ischar(x);
addRequired(p,'filename',validationFcn);
addRequired(p,'whichVariable',@ischar);
addOptional(p,'dateval',defaultDaterange,@isnumeric);
parse(p,filename,whichVariable,varargin{:});
daterange = p.Results.dateval;

% check to make sure input is either a .h5 or .mat file
if iscell(filename)
    filename = char(filename);
end
assert(strcmpi('.h5',filename(end-2:end)) ||...
    strcmpi('.mat',filename(end-3:end)),...
    'input must be HDF 5 or .mat file')
h5file = strcmpi('.h5',filename(end-2:end));

% other output arguments?
nout = max(nargout,1) - 1;
% set required and optional output to empty in case we bail
X = [];
for k=1:nout
    switch k
        case 1
            varargout{k} = []; %#ok<AGROW>
        case 2
            varargout{k} = struct([]); %#ok<AGROW>
    end
end

% check if variable supported, case insensitive, only the first four
% characters are considered
% need the variable's full name in order to retrieve
if length(whichVariable)>4
    Var = lower(whichVariable(1:4));
else
    Var = lower(whichVariable);
end
switch Var
    case 'snow'
        whichVariable = 'snow_fraction';
    case 'rock'
        whichVariable = 'rock_fraction';
    case 'soil'
        whichVariable = 'rock_fraction';
    case 'vege'
        whichVariable = 'vegetation_fraction';
    case 'shad'
        whichVariable = 'shade_fraction';
    case 'othe'
        whichVariable = 'other_fraction';
    case 'grai'
        whichVariable = 'grain_size';
    case 'drfs'
        whichVariable = 'drfs_grnsz';
    case 'delt'
        whichVariable = 'deltavis';
    case 'dust'
        whichVariable = 'dust';
    case 'rms'
        whichVariable = 'rms';
    otherwise
        warning('whichVariable = ''%s'' not recognized',whichVariable)
        return
end

% parse the date argument, if specified
if ~isempty(daterange)
    % MATLAB datenums or ISO dates, based on number of digits
    nDigits = ceil(log10(daterange(1)));
    assert(nDigits>=6 && nDigits<=8,...
        'input dates must be MATLAB datenums or 7- or 8-digit ISO dates')
    if nDigits==6
        daterange = round(daterange); % whole dates only
    else
        daterange = iso2datenum(daterange);
    end
end

% find the divisor
if h5file
    % check the group name for backward compatibility. Current name is
    % /Grid/MODIS_Grid_500m, but some old files may just be /Grid/500m.
    % This function will find either.
    group = findMODISh5group(filename,'500m');
    thisDataSet = [group '/' whichVariable];
    % try-catch sequence because could be 'divisor' or 'Divisor' or
    % none
    try
        divisor = h5readatt(filename,thisDataSet,'divisor');
    catch
        try
            divisor = h5readatt(filename,thisDataSet,'Divisor');
        catch
            divisor = 1;
        end
    end
else % .mat file
    m = matfile(filename);
    variableList = who(m);
    list = strcmpi('divisor',variableList);
    if nnz(list)>0
        if nnz(list)>1
            warning('more than 1 ''divisor'' field in file %s, using 1st one',filename)
            k = find(list>0,1,'first');
            divisor = m.(variableList{k});
        else
            divisor = m.(variableList{list});
        end
    else
        divisor = 1;
    end
end
% single date or multi-date
if isempty(daterange)
    alldates = true;
    if nout>=1
        varargout{1} = getDates(filename);
    end
else
    dval = getDates(filename);
    if isempty(dval)
        return
    elseif length(dval)>1
        alldates = false;
    else
        alldates = true;
        if nout>=1
            varargout{1} = dval;
        end
    end
end

% read the data
if alldates
    if h5file
        try
            data = h5read(filename,thisDataSet);
        catch
            warning('dataset %s not found in file %s',...
                thisDataSet, filename)
            return
        end
        info = h5info(filename,thisDataSet);
        fillvalue = info.FillValue;
    else
        vlq = who(m,whichVariable);
        if isempty(vlq)
            warning('dataset %s not found in file %s',...
                whichVariable, filename)
            return
        else
            data = m.(whichVariable);
        end
        fillvalue = [];
    end
else % retrieval by specific date(s)
    try
        info = h5info(filename,thisDataSet);
    catch
        warning('dataset %s not found in file %s',...
            thisDataSet, filename)
        return
    end
    siz = info.Dataspace.Size;
    fillvalue = info.FillValue;
    [idx,option,dates] = dateIndices(dval,daterange);
    if nout>=1
        varargout{1} = dates;
    end
    if isempty(idx)
        warning(['date(s) ' num2str(daterange) ' not found in file ' filename])
        return
    end
    switch option
        case 1
            siz(3) = 1;
            data = h5read(filename,thisDataSet,[1 1 idx],siz);
        case 2
            siz(3) = 1+idx(2)-idx(1);
            data = h5read(filename,thisDataSet,[1 1 idx(1)],siz);
        case 3
            siz(3) = 1;
            for k=1:length(idx)
                D = h5read(filename,thisDataSet,[1 1 idx(k)],siz);
                if k==1
                    data = zeros(siz(1),siz(2),length(idx),'like',D);
                end
                data(:,:,k) = D;
            end
    end
end

X = single(data)/divisor;
if fillvalue~=0
    X(data==fillvalue) = NaN;
end

if nout>=2
    if h5file
        hdr = GetCoordinateInfo(filename,group,size(X));
        varargout{2} = hdr;
    else
        matObj = matfile(filename);
        vlq = who(matObj,'hdr');
        if ~isempty(vlq)
            h = matObj.hdr;
            vl = fieldnames(h);
            % some raster reference fields empty in some files, so fix all
            list = strcmpi('rasterreference',vl);
            if nnz(list)>0
                h = rmfield(h,vl{list});
            end
            if strcmp('projected',h.gridtype)
                h.RasterReference =...
                    refmatToMapRasterReference(h.RefMatrix,[size(X,1) size(X,2)]);
            elseif strcmp('geographic',h.gridtype)
                h.RasterReference =...
                    refmatToGeoRasterReference(h.RefMatrix,[size(X,1) size(X,2)]);
            end
            varargout{2} = h;
        end
    end
end
end

function [idx,option,dates] = dateIndices(dval,daterange)
% returns indices of dates to retrieve
% option = 1 - single date, single value of idx
% option = 2 - date range, idx vector of length 2
% option = 3 - specific dates, idx vector of length >= 2
if length(daterange)==1
    idx = find(dval==daterange);
    option = 1;
    dates = daterange;
elseif length(daterange)==2 && daterange(1)<daterange(2)
    idx1 = find(dval==daterange(1));
    if isempty(idx1)
        idx1 = 1;
    end
    idx2 = find(dval==daterange(2));
    if isempty(idx2)
        idx2 = length(dval);
    end
    idx = [idx1 idx2];
    option = 2;
    dates = dval(idx1:idx2);
else % vector of dates
    daterange = sort(daterange);
    m = 1;
    option = 3;
    for k=1:length(daterange)
        id = find(dval==daterange(k));
        if ~isempty(id)
            idx(m) = id; %#ok<AGROW>
            dates(m) = dval(id); %#ok<AGROW>
            m = m+1;
        end
    end
end
if isrow(dates)
    dates = dates';
end
end

function dateval=getDates(filename)
% get single or multiple MATLAB dates from file
h5file = strcmpi('.h5',filename(end-2:end));
dateval = []; % in case of error
if h5file
    try
        dateval = h5readatt(filename,'/','MATLABdate');
    catch
        try
            dateval = h5readatt(filename,'/','MATLABdates');
        catch
            warning('apparently no date(s) in file %s',filename)
            return
        end
    end
else % .mat file
    m = matfile(filename);
    varList = who(m);
    list = strcmpi('MATLABdate',varList) | ...
        strcmpi('MATLABdates',varList) | ...
        strcmpi('datenum',varList) | ...
        strcmpi('matdate',varList);
    if nnz(list)==1
        dateval = m.(varList{list});
    elseif nnz(list)>1
        warning('more than 1 field matching a date variable in file %s',...
            filename)
        k = find(list>0,1,'first');
        dateval = m.(varList{k});
    else
        warning('apparently no date(s) in file %s',filename)
        return
    end
end
dateval = round(dateval); % whole days only
if isrow(dateval)
    dateval = dateval';
end
end