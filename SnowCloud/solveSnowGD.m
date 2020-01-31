function [grainSize,stats,varargout] = solveSnowGD(cosZ,reflectance,lambda,varargin )
% [grainsize,stats [, contamConc] ] = solveSnowGD(cosZ,reflectance,lambda,varargin )
%solves for grain size and, optionally, contaminant concentration given
%snow reflectance in one or more band passes
%
%Input
%   cosZ - cosine of solar illumination angle (scalar or vector of length
%       of reflectance)
%   reflectance - scalar or vector of reflectance values at angle arccos(cosZ)
%   lambda - if reflectance is scalar, vector of size [1 2] of wavelength
%       range, or if reflectance is a vector of length N, then bandpass is
%       a matrix of size [N 2]
%
%Optional input, name-value pairs, case insensitive
%   'GSunits' - units for grain size, default um
%   'lambdaUnits' - units for wavelength, default um
%   'contam' - contaminant, either 'dust' or 'soot'
%   'particulateSize' - effective radius of particulates, same 'units' as for
%       grain size (default 2.5 um for dust, 0.3 um for soot)
%
%Output
%   grainSize - in units specified, 'um' default
%   stats - statistics on solution
%Optional output
%   particulateConc - concentration of dust or soot

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

p = inputParser;
validZ = @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1;
validR = @(x) isnumeric(x) && all(x(:)>0) && all(x(:)<1);
validB = @(x) isnumeric(x) && size(x,2)==2 && all(x(:)>0);
addRequired(p,'cosZ',validZ)
addRequired(p,'reflectance',validR)
addRequired(p,'lambda',validB)
defaultContam = 'neither';
S = SnowCloudLimits();
defaultUnits = S.unitsSize;
addParameter(p,'gsunits',defaultUnits,@ischar)
addParameter(p,'lambdaunits',defaultUnits,@ischar)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'contamradius',[],@isnumeric)
parse(p,cosZ,reflectance,lambda,varargin{:})

% variables to pass to solver
cosineZ = p.Results.cosZ;
GSunits = p.Results.gsunits;
lambdaUnits = p.Results.lambdaunits;
Reflectance = p.Results.reflectance(:);
BandPass = sort(p.Results.lambda,2);
Contaminant = p.Results.contam;
ParticulateSize = p.Results.contamradius;
% make sure wavelenghts long enough to get grain size
assert(convertUnits(max(BandPass(:)),lambdaUnits,'nm')>=1060,...
    'maximum wavelength must be >= %f %s to retrieve grain size',...
    convertUnits(1060,'nm',lambdaUnits),lambdaUnits)
% could be same time of day or different times
assert(isscalar(cosineZ) || isequal(size(cosineZ),size(Reflectance)),...
    'if not scalar, cosZ must be same size as reflectance vector')

% default values for bounds
defaultGSbounds = convertUnits(S.snowRadius,defaultUnits,GSunits);
% initial guess use mean of grain size range
% reflectance close to linear w sqrt(grainSize) so use that
initialGS = sqrt(mean(defaultGSbounds));

% solving method depends on input
if strcmpi(p.Results.contam,defaultContam) %clean snow, just one unknown
    if nargout>2
        warning('if ''contam'' not specified, then just 2 output variables')
        varargout{1} = [];
    end
    
    % one unknown, one measurement, so exact solution
    if numel(Reflectance)==1 && size(BandPass,1)==1
        [x,fval,exitflag,output] = fzero(@cleanSnow,sqrt(initialGS));
        grainSize = x^2;
        stats.fval = fval;
        stats.exitflag = exitflag;
        stats.output = output;
    elseif numel(Reflectance)==size(BandPass,1) % least-squares solution
        lb = sqrt(min(defaultGSbounds));
        ub = sqrt(max(defaultGSbounds));
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@cleanSnow,sqrt(initialGS),lb,ub);
        grainSize = x^2;
        stats.resnorm = resnorm;
        stats.residual = residual;
        stats.exitflag = exitflag;
        stats.output = output;
        stats.lambda = lam;
        stats.jacobian = jacobian;
    else
        error('size(lambda,1) must == numel(reflectance)')
    end
else % dirty snow
    assert(strcmpi(Contaminant,'dust') || strcmpi(Contaminant,'soot'),...
        '''contam'' must be ''dust'' or ''soot''')
    if isempty(ParticulateSize)
        if strcmpi(Contaminant,'dust')
            ParticulateSize = S.dustRadius;
        else
            ParticulateSize = S.sootRadius;
        end
    end
    assert(nargout==3,'you need an output variable for the dust/soot concentration')
    % make sure wavelenghts short enough to get particulates
    assert(convertUnits(min(BandPass(:)),lambdaUnits,'nm')<=700,...
        'minimum wavelength must be <= %f %s to retrieve particulate concentration',...
        convertUnits(700,'nm',lambdaUnits),lambdaUnits)
    %initial guess
    if strcmpi(p.Results.contam,'dust')
        x0 = [sqrt(initialGS) mean(S.dust)];
        lb = [sqrt(min(defaultGSbounds)) min(S.dust)];
        ub = [sqrt(max(defaultGSbounds)) max(S.dust)];
    else
        x0 = [sqrt(initialGS) mean(S.soot)];
        lb = [sqrt(min(defaultGSbounds)) min(S.soot)];
        ub = [sqrt(max(defaultGSbounds)) max(S.soot)];
    end
    
    if numel(Reflectance)==size(BandPass,1) % least-squares solution
        [x,resnorm,residual,exitflag,output,lam,jacobian] =...
            lsqnonlin(@dirtySnow,x0,lb,ub);
        grainSize = x(1)^2;
        varargout{1} = x(2);
        stats.resnorm = resnorm;
        stats.residual = residual;
        stats.exitflag = exitflag;
        stats.output = output;
        stats.lambda = lam;
        stats.jacobian = jacobian;
    else
        error('size(lambda,1) must == numel(reflectance)')
    end
end

    function diffR = cleanSnow(sqrtGS)
        if numel(Reflectance)==1
            refl = SnowCloudBandPassReflectance(cosineZ,sqrtGS^2,GSunits,...
                BandPass, lambdaUnits);
            diffR = refl-Reflectance;
        else
            refl = zeros(numel(Reflectance),1);
            for k=1:numel(Reflectance)
                if isscalar(cosineZ)
                    mu0 = cosineZ;
                else
                    mu0 = cosineZ(k);
                end
                refl(k) = SnowCloudBandPassReflectance(mu0,sqrtGS^2,GSunits,...
                    BandPass(k,:), lambdaUnits);
            end
            diffR = refl-Reflectance;
        end
    end

    function diffR = dirtySnow(x)
        GS = x(1)^2;
        conc = x(2);
        refl = zeros(numel(Reflectance),1);
        for k=1:numel(Reflectance)
            if isscalar(cosineZ)
                mu0 = cosineZ;
            else
                mu0 = cosineZ(k);
            end
            refl(k) = SnowCloudBandPassReflectance(mu0,GS,GSunits,...
                BandPass(k,:),lambdaUnits,'contam',Contaminant,...
                'contamradius',ParticulateSize,'contamConc',conc);
        end
        diffR = refl-Reflectance;
    end
end