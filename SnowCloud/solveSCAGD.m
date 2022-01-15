function S = solveSCAGD(cosZ,sensor,bands,pixelReflectance,backgroundReflectance,varargin)
% S = solveSCAGD(cosZ,sensor,bands,pixelReflectance,backgroundReflectance,varargin)
%solves for fSCA, fBack, fShade, grain size, and contaminant concentration of snow
%
%Input
%   cosZ - cosine illumination angle
%   sensor - e.g., 'modis', 'LandsatOLI', or any of the sensors in
%       SensorTable.m
%   bands - cell or numeric vector of bands to use (numbers are character
%       strings in some sensors that also use letters)
%   pixelReflectance - reflectance of the pixel corresponding to the bands
%   backgroundReflectance - reflectance of the non-snow endmember
%Optional input, name-value pairs
%   'fShade' 0 to 1, though unlikely to be 1, specify if shade fraction
%       known, for example from topographic shading
%   'units' followed by any length unit for the output grainSize,
%       typically 'mm' or 'um', default 'um'
%   'scale', true or false, if true scales reflectance to the mean of the
%       bands, if false (the default) uses absolute reflectances
%       (if false, reflectances need to be adjusted for local illumination angle)
%   'contam', followed by 'soot' or 'dust' (program deals with either, but
%       not both in combination, default is 'neither', in which case snow is
%       assumed clean)
%   'contamRadius', followed by effective radius of contaminant, same 'units'
%       as for grain size (default 1 um for dust, 0.1 um for soot)
%   'startguess', followed by vector of 2 or 3 for initial guesses--2 if
%       clean snow, 3 otherwise--as [fSCA grainSize contamConc]
%   'method', solution method, either 'lsqnonlin' (default), 'fmincon', or
%       'fminsearch'
%
%Output, a structure
%   fSCA - fractional snow cover
%   fBack - fractional background
%   fShade - (1-fSCA-fBack)
%   fSCAnorm - fSCA normalized for shade
%   fBackNorm - fBack normalized for shade
%   grainSize - in units specified
%   contamConc - concentration of contaminant as a mass fraction
%   solvedReflectance - pixel reflectance at the solution
%   L - additional output from lsqnonlin or fmincon, whichever is used

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

p = inputParser;
defaultUnits = 'um';
defaultScale = false;
defaultDustRadius = 1.e-3;
defaultSootRadius = 1.e-4;
defaultContam = 'neither';
defaultMethod = 'lsqnonlin';
defaultShade = [];
defaultX0 = [];
validationFcn = @(x) iscell(x) || isnumeric(x);
rangeFcn = @(x) isnumeric(x) && all(x<=1) && all(x>=0);
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'sensor',@ischar)
addRequired(p,'bands',validationFcn)
addRequired(p,'pixelReflectance',@isnumeric)
addRequired(p,'backgroundReflectance',@isnumeric)
addParameter(p,'units',defaultUnits,@ischar)
addParameter(p,'scale',defaultScale,@islogical)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'contamradius',0,@isnumeric)
addParameter(p,'startguess',defaultX0,@isnumeric)
addParameter(p,'method',defaultMethod,@ischar)
addParameter(p,'fShade',defaultShade,rangeFcn)

parse(p,cosZ,sensor,bands,pixelReflectance,backgroundReflectance,varargin{:})

% fShade known or unknown
if isempty(p.Results.fShade)
    passShadeKnown = false;
else
    passShadeKnown = true;
    passShadeValue = p.Results.fShade;
end

% bands must be a string function
if iscell(p.Results.bands)
    bands = p.Results.bands;
else
    bands = cell(size(p.Results.bands));
    for k=1:length(p.Results.bands)
        bands{k} = num2str(p.Results.bands(k));
    end
end

% check inputs
assert(isequal(length(bands),length(pixelReflectance),...
    length(backgroundReflectance)),...
    'number of bands must be same as number of pixel and background reflectances')
doScale = p.Results.scale;
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
passB = bands(:)';
pixR = pixelReflectance(:)';
backR = backgroundReflectance(:)';
passSensor = sensor;
passCosine = cosZ;
passContam = p.Results.contam;
passContamRadius = contamRadius;

% lower and upper bounds
S = SnowCloudLimits();
if cleanSnow
    lb = [0 0 sqrt(S.snowRadius(1))];
    ub = [1 1 sqrt(S.snowRadius(2))];
else
    switch p.Results.contam
        case 'dust'
            lb = [0 0 sqrt(S.snowRadius(1)) S.dust(1)];
            ub = [1 1 sqrt(S.snowRadius(2)) S.dust(2)];
        case 'soot'
            lb = [0 0 sqrt(S.snowRadius(1)) S.soot(1)];
            ub = [1 1 sqrt(S.snowRadius(2)) S.soot(2)];
    end
end

% starting guess: x(1)=fSCA, x(2)=fBack x(3)=sqrt(grainSize), x(4)=contamConc
% (sqrt transform because reflectance is closer to linear in sqrt(r))
if ~isempty(p.Results.startguess) % x0 provided
    x0 = p.Results.startguess(:)';
    % convert units if necessary
    if ~strcmp(p.Results.units,defaultUnits)
        x0(3) = convertLengthUnits(x0(3),p.Results.units,defaultUnits);
    end
    x0(3) = sqrt(x0(3));
    if cleanSnow
        assert(length(x0)==3,...
            'for clean snow, starting guesses are fSCA, fBack, and grain size')
    else
        assert(length(x0)==4,...
            'for dirty snow, starting guesses are fSCA, fBac, grain size, and contaminant concentration')
    end
    % check validity of starting guess
    for k=1:length(x0)
        assert(x0(k)<=ub(k) && x0(k)>=lb(k),...
        ['some starting values out of range, x0=' num2str(x0)])
    end
else % use default starting guesses
    if cleanSnow
        x0 = [.5 .4 sqrt(mean(S.snowRadius))];
    else
        switch p.Results.contam
            case 'dust'
                x0 = [.5 .5 sqrt(mean(S.snowRadius)) mean(S.dust)];
            case 'soot'
                x0 = [.5 .5 sqrt(mean(S.snowRadius)) mean(S.soot)];
            otherwise
                error('''contam'' must be ''dust'' or ''soot''')
        end
    end
end

% solution
switch p.Results.method
    case 'lsqnonlin'
        options = optimoptions('lsqnonlin','FiniteDifferenceType','forward',...
            'Display','iter');
        [x,resnorm,res,eflag,output] = lsqnonlin(@compareSoln,x0,lb,ub,options);
        L.resnorm = resnorm;
        L.res = res;
        L.eflag = eflag;
        L.output = output;
    case 'fmincon'
        options = optimoptions('fmincon','FiniteDifferenceType','forward',...
            'Display','iter');
        [x,fval,exitflag,output,lambda,grad,hessian] =...
            fmincon(@sumsq,x0,[],[],[],[],lb,ub,[],options);
        L.fval = fval;
        L.exitflag = exitflag;
        L.output = output;
        L.lambda = lambda;
        L.grad = grad;
        L.hessian = hessian;
    case 'fminsearch'
        options = optimset('Display','iter');
        [x,fval,exitflag,output] = fminsearch(@sumsq,x0,options);
        L.fval = fval;
        L.exitflag = exitflag;
        L.output = output;
    otherwise
        error('''method'' is %s; choices are ''lsqnonlin'', ''fmincon'' or ''fminsearch''',...
            p.Results.method)
end
grainSize = x(3)^2;
fSCA = x(1);
fBack = x(2);
fShade = 1-fSCA-fBack;
fSCAnorm = fSCA/(1-fShade);
fBackNorm = fBack/(1-fShade);
if ~cleanSnow
    contamConc = x(4);
else
    contamConc = 0;
end

% reflectance values at solution
doScale = false;
solvedReflectance = mixedR(x);
if iscolumn(solvedReflectance)
    solvedReflectance = solvedReflectance';
end
if ~strcmp(p.Results.units,defaultUnits)
    grainSize = convertLengthUnits(grainSize,defaultUnits,p.Results.units);
end

%fill structure
S.fSCA = fSCA;
S.fBack = fBack;
S.fShade = fShade;
S.fSCAnorm = fSCAnorm;
S.fBackNorm = fBackNorm;
S.grainSize = grainSize;
S.contamConc = contamConc;
S.solvedReflectance = solvedReflectance;
S.L = L;

    function R = snowR(x)
        % snow reflectance in all bands for given grain radius
        r = x(3)^2;
        if cleanSnow
            R = SnowCloudSensorReflectance(passCosine,r,defaultUnits,...
                passSensor,passB,'lookup',true);
        else
            conc = x(4);
            R = SnowCloudSensorReflectance(passCosine,r,defaultUnits,...
                passSensor,passB,'contam',passContam,'contamRadius',passContamRadius,...
                'contamConc',conc,'lookup',true);
        end
    end

    function R = mixedR(x)
        % reflectance of mixed pixel, x(1)=fSCA, x(2)=grainSize,
        % x(3)=contamConc, x(4)=fBack
        sca = x(1);
        back = x(2);
        refl = snowR(x);
        if passShadeKnown
            shade = passShadeValue;
            back = 1-sca-shade;
        else
            shade = 1-sca-back; %#ok<NASGU>
        end
        R = sca*refl+back*backR;
    end

    function rdiff = compareSoln(x)
        if doScale
            M = mixedR(x);
            rdiff = (pixR-M)/(pixR+M);
        else
            rdiff = pixR-mixedR(x);
        end
    end

    function s = sumsq(x)
        s = sum(compareSoln(x).^2);
    end

end