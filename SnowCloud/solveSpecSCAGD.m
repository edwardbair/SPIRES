function outS = solveSpecSCAGD(cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin)
% S = solveSpecSCAGD(cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin)
%solves for grainSize and fSCA of clean snow
%
%Input
%   cosZ - cosine illumination angle
%   lambda - vector of wavelengths
%   lambdaUnits - typically 'um', 'mum', or 'nm'
%   pixelReflectance - reflectance of the pixel corresponding to the bands
%   backgroundReflectance - reflectance of the non-snow endmember
%Optional input, name-value pairs
%   'radiusUnits' followed by any length unit for the output grainSize,
%       typically 'mm' or 'um', default 'um'
%   'contam', followed by 'soot' or 'dust' (program deals with either, but
%       not both in combination, default is 'neither', in which case snow is
%       assumed clean)
%   'contamRadius', followed by effective radius of contaminant, same 'radiusUnits'
%       as for grain size (default 1 um for dust, 0.1 um for soot)
%   'scale', true or false, to scale reflectances by their sum, default
%       false (option to be deprecated)
%   'muTimesS', product mu.*cosd(S), if a scalar it's the same for both pixel
%       reflectance and background reflectance, it it's a vector the first
%       element is for pixel reflectance, the second for background, or if
%       it's a 2x2 matrix then the first column are the bounding values for
%       pixel reflectance and the second column are for background
%       reflectance (default is [0.1 0.1; 1 1])
%   'startguess', followed by vector of 2 or 3 for initial guesses--2 if
%       clean snow, 3 otherwise--as [fSCA grainSize contamConc]
%   'method', solution method, either 'lsqnonlin' (default), 'fmincon', or
%       'fminsearch'
%   'lookup', use lookup table for Mie and radiative transfer calculations
%       (default true)
%
%Output, structure with fields
%   fSCA - fractional snow cover
%   grainSize - in units specified
%   contamConc - concentration of contaminant as a mass fraction
%   solvedReflectance - pixel reflectance at the solution
%   L - optional output from lsqnonlin or fmincon, whichever is used

persistent passWeight

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

p = inputParser;
defaultUnits = 'mum';
defaultDustRadius = 1;
defaultSootRadius = 0.1;
defaultContam = 'neither';
defaultMethod = 'lsqnonlin';
defaultX0 = [];
default_muTimesS = [.1 .1; 1 1];
rangeFcn = @(x) isnumeric(x) && all(x(:)<=1) && all(x(:)>=0);
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'lambda',positiveFcn)
addRequired(p,'lambdaUnits',@ischar)
addRequired(p,'pixelReflectance',positiveFcn)
addRequired(p,'backgroundReflectance',positiveFcn)
addParameter(p,'radiusUnits',defaultUnits,@ischar)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'contamradius',0,positiveFcn)
addParameter(p,'scale',false,@islogical)
addParameter(p,'startguess',defaultX0,@isnumeric)
addParameter(p,'method',defaultMethod,@ischar)
addParameter(p,'muTimesS',default_muTimesS,@isnumeric)
addParameter(p,'lookup',true,@islogical)

parse(p,cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin{:})
assert(isscalar(cosZ),'''cosZ'' must be scalar')
cosZ = double(cosZ);

% check inputs
[lambda,pixelReflectance,backgroundReflectance] =...
    checkSizes(lambda,pixelReflectance,backgroundReflectance);
assert(isrow(lambda) || iscolumn(lambda),...
    'lambda must be row or column vector (maybe fixed in later version)')

passDoScale = p.Results.scale;
if passDoScale
    warning('''scale'' option will be deprecated, replaced by ''muTimesS''')
end

contaminant = p.Results.contam;
switch contaminant
    case defaultContam
        cleanSnow = true;
        contamRadius = p.Results.contamradius;
    case 'dust'
        cleanSnow = false;
        if p.Results.contamradius==0
            contamRadius = defaultDustRadius;
        else
            contamRadius = p.Results.contamradius;
        end
    case 'soot'
        cleanSnow = false;
        if p.Results.contamradius==0
            contamRadius = defaultSootRadius;
        else
            contamRadius = p.Results.contamradius;
        end
    otherwise
        error('''contam'' ''%s'' not recognized',p.Results.contam)
end
% make all the vectors row vectors
pixR = double(pixelReflectance(:)');
backR = double(backgroundReflectance(:)');
passCosine = cosZ;
passContam = p.Results.contam;
passContamRadius = contamRadius;
passLambda = lambda(:)';
passLambdaUnits = lambdaUnits;
passLookup = p.Results.lookup;

% lower and upper bounds
S = SnowCloudLimits();
if cleanSnow
    lb = [0 sqrt(S.snowRadius(1))];
    ub = [1 sqrt(S.snowRadius(2))];
else
    switch p.Results.contam
        case 'dust'
            lb = [0 sqrt(S.snowRadius(1)) S.dust(1)];
            ub = [1 sqrt(S.snowRadius(2)) S.dust(2)];
        case 'soot'
            lb = [0 sqrt(S.snowRadius(1)) S.soot(1)];
            ub = [1 sqrt(S.snowRadius(2)) S.soot(2)];
    end
end

% if muTimesS is a matrix, then it's among the unknowns
assert(max(p.Results.muTimesS(:))<=max(default_muTimesS(:)) &&...
    min(p.Results.muTimesS(:))>=min(default_muTimesS(:)),...
    '''muTimesS'' must be >=%f and <=%f',...
    min(default_muTimesS(:)),max(default_muTimesS(:)))
if isequal(size(p.Results.muTimesS),[2 2])
    lb = cat(2,lb,p.Results.muTimesS(1,:));
    ub = cat(2,ub,p.Results.muTimesS(2,:));
elseif isscalar(p.Results.muTimesS)
    passMuTimesS = p.Results.muTimesS*ones(1,2);
elseif numel(p.Results.muTimesS)==2
    passMuTimesS = [p.Results.muTimesS(1) p.Results.muTimesS(2)];
else
    error('''muTimesS'', if specified, must be a scalar, vector of length 2, or 2x2 matrix')
end

% weights based on differences across grain size and contaminant amounts for snow
% and between snow reflectance and the background
passWeight = spectralWeight(passCosine,passLambda,passLambdaUnits,backR,'dust');

% starting guess: x(1)=fSCA, x(2)=sqrt(grainSize), x(3)=contamConc
% (sqrt transform because reflectance is closer to linear in sqrt(r))
if ~isempty(p.Results.startguess) % x0 provided
    x0 = p.Results.startguess(:)';
    % convert units if necessary
    if ~strcmp(p.Results.radiusUnits,defaultUnits)
        x0(2) = convertLengthUnits(x0(2),p.Results.units,defaultUnits);
    end
    if cleanSnow
        assert(length(x0)==2,...
            'for clean snow, starting guesses are fSCA and grain size')
    else
        assert(length(x0)==3,...
            'for dirty snow, starting guesses are fSCA, grain size, and contaminant concentration')
    end
    % check validity of starting guess
    x0(2) = sqrt(x0(2));
    assert(x0(1)<=ub(1) && x0(1)>=lb(1) && x0(2)<=ub(2) && x0(2)>=lb(2),...
        ['some starting values out of range, x0=' num2str(x0)])
    if ~cleanSnow
        assert(x0(3)>=lb(3) && x0(3)<=ub(3),...
            ['some starting values out of range, x0=' num2str(x0)])
    end
else % use default starting guesses
    if cleanSnow
        x0 = [.5 sqrt(mean(S.snowRadius))];
    else
        switch p.Results.contam
            case 'dust'
                x0 = [.5 sqrt(mean(S.snowRadius)) mean(S.dust)];
            case 'soot'
                x0 = [.5 sqrt(mean(S.snowRadius)) mean(S.soot)];
            otherwise
                error('''contam'' must be ''dust'' or ''soot''')
        end
    end
end
% add initialization if solving for muTimesS
if length(lb)>length(x0)
    x0 = cat(2,x0,(lb(end-1)+ub(end-1))/2,(lb(end)+ub(end))/2);
end

% solution
switch p.Results.method
    case 'lsqnonlin'
        options = optimoptions('lsqnonlin','FiniteDifferenceType','forward',...
            'Display','off');
        [x,resnorm,res,eflag,output] = lsqnonlin(@compareSoln,x0,lb,ub,options);
        L.resnorm = resnorm;
        L.res = res;
        L.eflag = eflag;
        L.output = output;
    case 'fmincon'
        options = optimoptions('fmincon','FiniteDifferenceType','forward',...
            'Display','iter');
        [x,fval,exitflag,output,flambda,grad,hessian] =...
            fmincon(@sumsq,x0,[],[],[],[],lb,ub,[],options);
        L.fval = fval;
        L.exitflag = exitflag;
        L.output = output;
        L.lambda = flambda;
        L.grad = grad;
        L.hessian = hessian;
    case 'fminsearch'
        options = optimset('Display','iter');
        [x,fval,exitflag,output] = fminsearch(@sumsq,x0,options);
        L.fval = fval;
        L.exitflag = exitflag;
        L.output = output;
    otherwise
        error('''method'' is %s; choices are ''lsqnonlin'', ''fmincon''',...
            p.Results.method)
end
effectiveRadius = x(2)^2;
fSCA = x(1);
if ~cleanSnow
    contamConc = x(3);
else
    contamConc = 0;
end

% reflectance values at solution
solvedReflectance = mixedR(x);
if ~strcmp(p.Results.radiusUnits,defaultUnits)
    effectiveRadius = convertLengthUnits(effectiveRadius,defaultUnits,p.Results.units);
end

% results into structure
outS.fSCA = fSCA;
outS.effectiveRadius = effectiveRadius;
outS.contamConc = contamConc;
if length(x)>3
    outS.muTimesS = [x(end-1) x(end)];
else
    outS.muTimeS = passMuTimesS;
end
outS.pixelReflectanceScaled = pixR(:)*passCosine/outS.muTimesS(1);
outS.solvedReflectance = solvedReflectance(:);
snowRefl = snowR(x);
outS.snowReflectance = snowRefl(:);
outS.backgroundReflectanceScaled = backR(:)*passCosine/outS.muTimesS(2);
outS.method = p.Results.method;
outS.L = L;

    function R = snowR(x)
        % snow reflectance in all wavelengths for given grain radius
        r = x(2)^2;
        if cleanSnow
            R = SnowCloudSpectralReflectance(passCosine,r,defaultUnits,...
                passLambda,passLambdaUnits,'lookup',passLookup);
        else
            conc = x(3);
            R = SnowCloudSpectralReflectance(passCosine,r,defaultUnits,...
                passLambda,passLambdaUnits,'contam',passContam,'contamRadius',passContamRadius,...
                'contamConc',conc,'lookup',passLookup);
        end
    end

    function R = mixedR(x)
        % reflectance of mixed pixel, x(1)=fSCA, x(2)=grainSize, x(3)=contamConc
        % corrected for topography (eq 13b Dozier 1989, with no diffuse)
        % Reff = R0*cos(S)*mu/mu0, so R0 = mu0*Reff/(cos(S)*mu))
        sca = x(1);
        refl = snowR(x);
        % adjust for angles
        if length(x)>3
            passMuTimesS = [x(end-1) x(end)];
        end
        % convert backR (i.e., Reff) to R0 of background
        R = sca*refl+(1-sca)*backR*passCosine/passMuTimesS(2);
    end

    function rdiff = compareSoln(x)
        % adjust for angles
        if length(x)>3
            passMuTimesS = [x(end-1) x(end)];
        end
        % convert simulated reflectance (i.e. R0) to Reff
        mix = mixedR(x)*passMuTimesS(1)/passCosine;
        if passDoScale
            rdiff = ((pixR-mix)./(pixR+mix)).*passWeight;
        else
            rdiff = (pixR-mix).*passWeight;
        end
        rdiff = double(rdiff);
    end

    function s = sumsq(x)
        s = sum(compareSoln(x).^2);
    end

end