function outS = solveSpecAbsorption(cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin)
% outS = solveSpecAbsorption(cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin)
%solves for grainSize and fSCA of snow based on absorption regions
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
%   'startguess', followed by vector of 2 or 3 for initial guesses--2 if
%       clean snow, 3 otherwise--as [fSCA grainSize contamConc]
%   'method', solution method, either 'lsqnonlin' (default), 'fmincon', or
%       'fminsearch'
%   'lookup', use lookup table for Mie calculations
%       (default true)
%
%Output, structure with fields
%   fSCA - fractional snow cover
%   effectiveRadius - in units specified
%   contamRadius - effective radius of particulate contaminant
%   contamConc - concentration of contaminant as a mass fraction
%   solvedReflectance - only in the wavelengths of ScaledSnowAbsorption
%   reflDiff - only in the wavelengths of ScaledSnowAbsorption
%   L - optional output from lsqnonlin or fmincon, whichever is used

warning('function %s deprecated, use the toolbox/SnowCloudReflectance instead', mfilename)

p = inputParser;
defaultUnits = 'mum';
defaultContam = 'neither';
defaultMethod = 'lsqnonlin';
rangeFcn = @(x) isnumeric(x) && all(x(:)<=1) && all(x(:)>=0);
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
addRequired(p,'cosZ',rangeFcn)
addRequired(p,'lambda',positiveFcn)
addRequired(p,'lambdaUnits',@ischar)
addRequired(p,'pixelReflectance',positiveFcn)
addRequired(p,'backgroundReflectance',positiveFcn)
addParameter(p,'radiusUnits',defaultUnits,@ischar)
addParameter(p,'contam',defaultContam,@ischar)
addParameter(p,'method',defaultMethod,@ischar)
addParameter(p,'lookup',true,@islogical)

parse(p,cosZ,lambda,lambdaUnits,pixelReflectance,backgroundReflectance,varargin{:})
assert(isscalar(cosZ),'''cosZ'' must be scalar')
cosZ = double(cosZ);

% check inputs
[lambda,pixelReflectance,backgroundReflectance] =...
    checkSizes(lambda,pixelReflectance,backgroundReflectance);
assert(isrow(lambda) || iscolumn(lambda),...
    'lambda must be row or column vector (maybe fixed in later version)')

switch p.Results.contam
    case defaultContam
        cleanSnow = true;
    case 'dust'
        cleanSnow = false;
    case 'soot'
        cleanSnow = false;
    otherwise
        error('''contam'' ''%s'' not recognized',p.Results.contam)
end
% make all the vectors row vectors
pixR = double(pixelReflectance(:)');
backR = double(backgroundReflectance(:)');
passCosine = cosZ;
passContam = p.Results.contam;
passLambda = lambda(:)';
passLambdaUnits = lambdaUnits;
passLookup = p.Results.lookup;

% lower and upper bounds
S = SnowCloudLimits();
passRadiusUnits = S.unitsSize;
if cleanSnow
    lb = [0 sqrt(S.snowRadius(1))];
    ub = [1 sqrt(S.snowRadius(2))];
else
    switch p.Results.contam
        case 'dust'
            lb = [0 sqrt(S.snowRadius(1)) S.dust(1) sqrt(S.dustRadius(1))];
            ub = [1 sqrt(S.snowRadius(2)) S.dust(2) sqrt(S.dustRadius(2))];
        case 'soot'
            lb = [0 sqrt(S.snowRadius(1)) S.soot(1) sqrt(S.sootRadius(1))];
            ub = [1 sqrt(S.snowRadius(2)) S.soot(2) sqrt(S.sootRadius(2))];
    end
end

% starting guess: x(1)=fSCA, x(2)=sqrt(grainSize), x(3)=contamConc, x(4)=sqrt(contamRadius)
% (sqrt transform because reflectance is closer to linear in sqrt(r))

% use default starting guesses
if cleanSnow
    x0 = [.5 sqrt(mean(S.snowRadius))];
else
    switch p.Results.contam
        case 'dust'
            x0 = [.5 sqrt(mean(S.snowRadius)) mean(S.dust) sqrt(mean(S.dustRadius))];
        case 'soot'
            x0 = [.5 sqrt(mean(S.snowRadius)) mean(S.soot) sqrt(mean(S.sootRadius))];
        otherwise
            error('''contam'' must be ''dust'' or ''soot''')
    end
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
    contamRadius = x(4)^2;
else
    contamConc = 0;
    contamRadius = [];
end

% results into structure
outS.fSCA = fSCA;
outS.effectiveRadius = convertLengthUnits(effectiveRadius,passRadiusUnits,p.Results.radiusUnits);
outS.radiusUnits = p.Results.radiusUnits;
outS.wavelength = lambda;
outS.waveUnits = p.Results.lambdaUnits;
if ~cleanSnow
    outS.contamConc = contamConc;
    outS.contamRadius = convertLengthUnits(contamRadius,defaultUnits,p.Results.radiusUnits);
    outS.snowReflectance = SnowCloudSpectralReflectance(cosZ,outS.effectiveRadius,...
        p.Results.radiusUnits,lambda,p.Results.lambdaUnits,...
        'contam',p.Results.contam,'contamRadius',outS.contamRadius,...
        'contamConc', contamConc);
else
    outS.snowReflectance = SnowCloudSpectralReflectance(cosZ,outS.effectiveRadius,...
        p.Results.radiusUnits,lambda,p.Results.lambdaUnits);
end
outS.method = p.Results.method;
outS.L = L;

    function R = snowR(x)
        % snow reflectancefor given
        % grain radius and contaminant properties
        r = x(2)^2;
        if cleanSnow
            R = SnowCloudSpectralReflectance(passCosine,r,passRadiusUnits,...
                passLambda,passLambdaUnits);
        else
            R = SnowCloudSpectralReflectance(passCosine,r,passRadiusUnits,...
                passLambda,passLambdaUnits,'contam',passContam,...
                'contamConc',x(3),'contamRadius',x(4)^2);
        end
    end

    function mixV = mixedRdiff(x)
        % scaled reflectance difference of mixed pixel, x(1)=fSCA, x(2)=sqrt(grainSize),
        % x(3)=contamConc, x(4) sqrt(contamRadius)
        sca = x(1);
        modelRefl = sca*snowR(x)+(1-sca)*backR;
        mixedS = ScaledSnowAbsorptionActual(passLambda,passLambdaUnits,...
            modelRefl,'snowRadius',x(2)^2,'cosZ',passCosine,'lookup',passLookup);
        fn = fieldnames(mixedS);
        for k=1:length(fn)
            thisS = mixedS.(fn{k});
            if k==1
                mixV = [thisS.wavelength thisS.reflDiff];
            else
                mixV = [mixV; thisS.wavelength thisS.reflDiff]; %#ok<AGROW>
            end
        end
        mixV = sortrows(mixV);
    end

    function rdiff = compareSoln(x)
        % convert simulated reflectance difference to actual
        mixV = mixedRdiff(x);
        thisWave = mixV(:,1);
        thisR = pchip(passLambda,pixR,thisWave);
        compS = ScaledSnowAbsorptionActual(thisWave,passLambdaUnits,...
            thisR,'snowRadius',x(2)^2,'cosZ',passCosine,'lookup',passLookup);
        fn = fieldnames(compS);
        for k=1:length(fn)
            thisS = compS.(fn{k});
            if k==1
                dV = [thisS.wavelength thisS.reflDiff];
            else
                dV = [dV; thisS.wavelength thisS.reflDiff]; %#ok<AGROW>
            end
        end
        dV = sortrows(dV);
        rdiff = dV(:,2)-mixV(:,2);
    end

    function s = sumsq(x)
        s = sum(compareSoln(x).^2);
    end

end