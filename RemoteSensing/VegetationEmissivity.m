function [ emiss ] = VegetationEmissivity(veg,varargin)
% emiss = VegetationEmissivity(veg,lambda [,units])
% or VegetationEmissivity('help') to get list of vegetation types
%emissivity of specified vegetation
%
%Input
%   veg - type of vegetation or 'help' to get list
%   lambda - wavelength, scalar or vector/matrix (can omit if veg is 'help')
%Optional input
%   units - any length unit from 'nm' to 'm', 'mum' is default
%
%Output
%   emissivity corresponding to wavelengths

persistent lastVeg eInterp

warning('deprecated function, use ECOSTRESS_SpecRefl instead and subtract reflectance from 1')

p = inputParser;
positiveFcn = @(x) isnumeric(x) && all(x>0);
addRequired(p,'veg',@ischar)
addOptional(p,'lambda',[],positiveFcn)
addOptional(p,'units','mum',@ischar)
parse(p,veg,varargin{:})

% all vegetation types
allVegType = {'pine'};
emissFile = {'PineEmissivity.mat'};

if ~strcmpi(veg,lastVeg)
    switch p.Results.veg
        case 'help'
            emiss = [];
            disp(['function ' mfilename ' vegetation types'])
            disp(allVegType)
            return
        otherwise
            found = false;
            for k=1:length(allVegType)
                if strcmpi(p.Results.veg,allVegType{k})
                    load(emissFile{k});
                    found = true;
                    break;
                end
            end
            if ~found
                error('veg type ''%s'' not recognized',veg)
            end
            lastVeg = veg;
            % built interpolating function
            eInterp = griddedInterpolant(LeafEmissivity.wavelength,...
                LeafEmissivity.emissivity,'pchip','none');
    end
end

% wavelengths to micrometers
if isempty(p.Results.lambda)
    error('unless ''help'' is specified as veg, lambda must be specified too')
end
lambda = convertLengthUnits(p.Results.lambda,p.Results.units,'mum');

emiss = eInterp(lambda);

end

