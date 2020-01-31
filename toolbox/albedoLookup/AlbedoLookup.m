function broadbandAlbedo = AlbedoLookup(radius,cosineSolarZ,cosineIllum,elevation,varargin)
%broadbandAlbedo = AlbedoLookup(radius,cosineSolarZ,cosineIllum,elevation [,LAPname,LAPconc])
%Estimates broadband albedo of snow for mid-latitude winter atmosphere,
%corresponding to Bair et al. (2019):
%   Bair, E.H., Rittger, K., Skiles, S.M., and Dozier, J. (2019), An
%   examination of snow albedo estimates from MODIS and their impact on
%   snow water equivalent reconstruction, Water Resources Research
%
%Input values can be of any size and dimension. All must be the same size
%and dimension, except any can be a scalar. Either cosineSolarZ or
%cosineIllum can be empty, but not both, to use same values.
%   radius - effective optical radius of the snow grains, micrometers
%   cosineSolarZ - cosine of the solar zenith angle ([] for same as
%       cosineIllum)
%   cosineIllum - cosine of the local illumination angle ([] for same as
%       cosineSolarZ)
%   elevation - in kilometers
%Optional inputs
%   LAPname - either 'dust' or 'soot'
%   LAPconc - concentration of light-absorbing particles, as a mass ratio
%       so therefore dimensionless
%
%Output
%   broadbandAlbedo of same size and dimension as inputs
%   Range of lookup table includes radii from 30 to 2000 micrometers,
%   cosineSolarZ and cosineIllum from 0.05 to 1.0 (87 deg to 0 deg),
%   elevations from 1.5 to 4.5 km,dust concentrations from 0 to .001 (1000 ppm),
%   and soot concentrations from 0 to 5e-5 (50 ppm).
%   Albedo for values outside a range are set to that of the boundary.

persistent already S

narginchk(4,8)
nargoutchk(0,1)
assert(nargin~=5,'both LAPname and LAPconc must be specified, or neither')

%load lookup table
if isempty(already) || ~already
    already = true;
    S = load('albedoLUT.mat');
end

%parse input values
p = inputParser;
addRequired(p,'radius',@(x) isnumeric(x) && all(x(:)>0))
addRequired(p,'cosineSolarZ',@(x) isempty(x) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<=1)))
addRequired(p,'cosineIllum',@(x) isempty(x) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<=1)))
addRequired(p,'elevation',@isnumeric)
addOptional(p,'LAPname','',@(x) isempty(x) ||...
    (ischar(x) && (strcmpi(x,'soot') || strcmpi(x,'dust'))))
addOptional(p,'LAPconc',[],@(x) isempty(x) ||...
    (isnumeric(x) && all(x(:)>=0) && all(x(:)<1)))
parse(p,radius,cosineSolarZ,cosineIllum,elevation,varargin{:})

%dust or soot, or clean
if isempty(p.Results.LAPname)
    F = S.Fdust;
    LAPconc = 0;
else
    switch lower(p.Results.LAPname)
        case 'dust'
            F = S.Fdust;
        case 'soot'
            F = S.Fsoot;
        otherwise % shouldn't reach as should be caught by input parser
            error('LAPname %s not recognized',p.Results.LAPname)
    end
    LAPconc = p.Results.LAPconc;
end

%check input solar and illumination angles
assert(~(isempty(p.Results.cosineSolarZ) && isempty(p.Results.cosineIllum)),...
    'cosineSolarZ and cosineIllum cannot both be empty')
if isempty(p.Results.cosineSolarZ)
    cosineSolarZ = p.Results.cosineIllum;
    cosineIllum = p.Results.cosineIllum;
elseif isempty(p.Results.cosineIllum)
    cosineIllum = p.Results.cosineSolarZ;
    cosineSolarZ = p.Results.cosineSolarZ;
else
    cosineSolarZ = p.Results.cosineSolarZ;
    cosineIllum = p.Results.cosineIllum;
end

%check input sizes
[radius,cosZ,cosZp,elevation,LAPconc] =...
    checkSizes(p.Results.radius,cosineSolarZ,cosineIllum,p.Results.elevation,LAPconc);

%lookup albedo

broadbandAlbedo = F(radius,LAPconc,elevation,cosZ,cosZp);
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