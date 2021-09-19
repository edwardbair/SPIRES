function [broadbandAlbedo, nIRAlbedo, visAlbedo ] = AlbedoLookup(radius,...
    cosineSolarZ,elevation,varargin)
%
% Estimates of snow albedos for a mid-latitude winter atmosphere, v1.1
% This version updates the original code:
%
% EH Bair. Snow albedo lookup tables. Zenodo.
% http://doi.org/10.5281/zenodo.3228428
%
% referenced in Appendix A of
%   Bair, E.H., Rittger, K., Skiles, S.M., and Dozier, J. (2019), An
%   examination of snow albedo estimates from MODIS and their impact on
%   snow water equivalent reconstruction, Water Resources Research
%
% Note: Jeff Dozier wrote the original function which is mostly repeated
% here with a few changes
%
%Input values can be of any size and dimension. All must be the same size
%and dimension, except any can be a scalar.
%   radius - effective optical radius of the snow grains, micrometers
%   cosineSolarZ - cosine of the solar zenith angle
%   elevation - in kilometers
%Optional inputs, Name value pairs - assumed clean if skipped
%   LAPname - either 'dust' or 'soot'
%   LAPconc - concentration of light-absorbing particles, as a mass ratio
%       so therefore dimensionless
%
%Output
%   broadbandAlbedo (0.28 to 4 um)
%   nIRAlbedo (0.8 4 um)
%   visAlbedo (0.380 0.750 um)
%   All of same size and dimension as inputs
%
%   Range of lookup table includes radii from 30 to 1800 micrometer 
%   grain radius,
%   cosineSolarZ from 0.05 to 1.0 (87 deg to 0 deg),
%   elevations from 0 to 8 km, dust concentrations from 0 to .001 (1000 ppm),
%   and soot concentrations from 0 to 5e-5 (50 ppm).
%   Albedo for values outside a range are set to that of the boundary
%
% example usage: [broadbandAlbedo,nIRAlbedo,visAlbedo] = AlbedoLookup(500,...
%     0.8,3,'dust',200e-6)

persistent already S

narginchk(3,7)
nargoutchk(0,3)
assert(nargin~=4,'both LAPname and LAPconc must be specified, or neither')

%load lookup table
if isempty(already) || ~already
    already = true;
    S = load('albedoLUT.mat');
end
    F = S.F;

%parse input values
p = inputParser;
addRequired(p,'radius',@(x) isnumeric(x) && all(x(:)>0))
addRequired(p,'cosineSolarZ',@(x) isempty(x) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<=1)))
addRequired(p,'elevation',@isnumeric)
addOptional(p,'LAPname','',@(x) isempty(x) ||...
    (ischar(x) && (strcmpi(x,'soot') || strcmpi(x,'dust'))))
addOptional(p,'LAPconc',[],@(x) isempty(x) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<1)))
parse(p,radius,cosineSolarZ,elevation,varargin{:})

%dust (1st dim=1) or soot(1st dim=2), or clean
if isempty(p.Results.LAPname)
    dim1=1; %dust
    LAPconc = 0;
else
    switch lower(p.Results.LAPname)
        case 'dust'
            dim1=1;
        case 'soot'
            dim1=2;
        otherwise % shouldn't reach as should be caught by input parser
            error('LAPname %s not recognized',p.Results.LAPname)
    end
    LAPconc = p.Results.LAPconc;
end

%check input solar and illumination angles
assert(~(isempty(p.Results.cosineSolarZ)),...
    'cosineSolarZ cannot be empty')

%check input sizes
[radius,cosZ,elevation,LAPconc] =...
    checkSizes(p.Results.radius,cosineSolarZ,p.Results.elevation,LAPconc);

%lookup albedo
%grid vectors are:
%laptype(1=dust;2=soot),LAPConc,elevation,cosZ,bandPasses(1=broadband;2=nIR;3=vis),grainSize;

for i=1:3
    A= F(dim1,double(LAPconc),...
        double(elevation),double(cosZ),i,double(radius));
    switch i
        case 1
            broadbandAlbedo=A;
        case 2
            nIRAlbedo=A;
        case 3
            visAlbedo=A;
     end
end
end

function [ varargout ] = checkSizes( varargin )
% [ varargout ] = checkSizes( varargin )
%general function to check sizes of inputs
%all inputs that are not scalars must be same size
%all inputs that are scalars are expanded to same size as non-scalar input(s)

% check number of input and output arguments
narginchk(1,Inf)
nargoutchk(nargin,nargin)

% which ones not scalar (set to scalar if just replicated scalar)
scalar = false(size(varargin));
for k=1:length(varargin)
    scalar(k) = isscalar(varargin{k});
    if ~scalar(k)
        if isscalar(unique(varargin{k}(:)))
            scalar(k) = true;
            varargin{k} = unique(varargin{k}(:));
        end
    end
end

% trivial case, all scalar
if all(scalar)
    for k=1:length(varargin)
        varargout{k} = varargin{k}; %#ok<*AGROW>
    end
    
    % just one non-scalar, set others to same size
elseif sum(~scalar)==1
    n = find(~scalar);
    for k=1:length(varargin)
        if k==n
            varargout{k} = varargin{k};
        else
            varargout{k} = repmat(varargin{k},size(varargin{n}));
        end
    end
    
    % general case
else
    % vectors must be row or column, but not both
    for k=1:length(varargin)
        if isrow(varargin{k})
            varargin{k} = varargin{k}';
        end
    end
    n = find(~scalar);
    N = size(varargin{n(1)});
    for k=2:length(n)
        assert(isequal(N,size(varargin{n(k)})),...
            'variable %d has different size than variable %d',n(k),n(1))
    end
    for k=1:length(varargin)
        if scalar(k)
            varargout{k} = repmat(varargin{k},N);
        else
            varargout{k} = varargin{k};
        end
    end
end
end