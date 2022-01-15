function [oStruct,varargout] = SPIReS_inv(Rfun,substance,unknowns,P0,varargin )
% [oStruct] = SPIReS_inv(Rfun,unknowns,prescription,prop/val )
% [oStruct] = SPIReS_inv(Rfun,substance,unknowns,prescription,prop/val )
% [oStruct,stats] = SPIReS_inv(Rfun,unknowns,prop/val)
% [oStruct,stats,P] = SPIReS_inv(Rfun,unknowns,prop/val)
%
%solves for snow or cloud properties based on spectral reflectance
%   (i.e. the inverse of  SPIReS_fwd)
%
%Input
% Rfun - interpolant object to return measured reflectance from wavelength
% substance - 'snow','iceCloud','waterCloud','mixedCloud'
% unknowns - cell vector of snow or cloud properties to solve for, with
%       possibilities depending on the input substance and the number of
%       reflectance values available (reasonable abbreviations work)
%       for 'snow' - 'radius', 'fSCA', 'wetness', 'waterEquivalent',
%       'lapfraction', 'lapradius', 'muS', (to-be-added 'imagRI')
%       for 'iceCloud' - 'radius', 'waterEquivalent'
%       for 'waterCloud' - 'radius', 'waterEquivalent'
%       for 'mixedCloud' - 'radius', 'waterRadius', 'wetness', 'waterEquivalent'
% prescription - prescription that covers the information unrelated to the
%   snow or cloud properties along with the initial estimates of those properties
%
%Optional input
% 'atmosphere' - followed by a function that getAtmosProp uses to calculate
%   atmospheric transmission as functions of wavelength, elevation, and
%   solar zenith (if omitted, reflectance values are not weighted)
% 'method' - solution method to solve the inversion, either 'lsq' (default),
%   'normResiduals', or 'spectralAngle' (4-letter truncations okay)
% 'useParallel' - default false, use parallel processing in the optimization
%   (generally leave as default unless lots of calculations)
%
%Output
%   oStruct - snow or cloud properties, depending on inputs
%Optional output
%   stats - statistics about solution
%   P - prescription at last call to forward function

narginchk(4,10)
nargoutchk(0,3)

p = inputParser;
addRequired(p,'Rfun',@(x) contains(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
addRequired(p,'substance',@(x) ischar(x) || isstring(x) || iscategorical(x))
addRequired(p,'unknowns',@iscell)
addRequired(p,'P0',@isstruct)
addParameter(p,validatestring('atmosphere',{'atm','atmos','atmosphere'}),...
    [],@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
addParameter(p,'method','lsq',@ischar)
addParameter(p,'useparallel',false,@(x) isnumeric(x) || islogical(x))
parse(p,Rfun,substance,unknowns,P0,varargin{:});

% substances must match
if ischar(substance)
    substance = categorical({substance});
elseif isstring(substance)
    substance = categorical(substance);
end
assert(P0.substance==substance,...
    'input substance ''%s'' different than prescription substance ''%s''',...
    char(substance),char(P0.substance))

% parallel processing
useParallel = logical(p.Results.useparallel);

% parse unknowns
Spectrum = P0.Spectrum;
assert(length(unknowns)<=length(Spectrum.wavelength),...
    'length of reflectance vector must be >= number of unknowns')
method = p.Results.method;
if strcmpi(method,'lsq')
    method = 'leastSquares';
elseif contains(method,'norm','IgnoreCase',true)
    method = 'normResiduals';
elseif contains(method,'spec','IgnoreCase',true)
    method = 'spectralAngle';
end

[solveFor,prescription] = setUnknowns(substance,unknowns,P0);
Spectrum = prescription.Spectrum; % probably didn't change, but just in case

% stupid warning - leastSquares is better
% if strcmpi(method,'leastSquares') && any(contains(solveFor,'fsca','IgnoreCase',true))
%     warning('if ''fSCA'' is among the unknowns, ''method'' ''normResiduals'' is preferred')
% end

% intial values and limits
[x0,lb,ub,solveFor] = setBounds(solveFor,prescription,method);

% spectral weights
if ~isempty(p.Results.atmosphere)
    passWeight = getAtmosProp(p.Results.atmosphere,...
        prescription.Spectrum.wavelength,prescription.Illumination.elevation,...
        prescription.Illumination.cosZ);
else
    passWeight = ones(size(Spectrum.wavelength));
end

measuredReflectance = Rfun(Spectrum.wavelength);

switch method
    case 'leastSquares'
        functionToMinimize = @SnowCloudDiff;
    case 'normResiduals'
        functionToMinimize = @SnowCloudNormResiduals;
    case 'spectralAngle'
        functionToMinimize = @SnowCloudSpectralAngle;
    otherwise
        error('solution ''method'' %s not recognized',method)
end
passSolver = method;
mR = []; % initialize so it's a passed variable
% for all cases, use central differences
options = optimset('Display','off','FinDiffType','central',...
    'UseParallel',useParallel);
if strcmpi(method,'leastSquares')
    [x,fval,~,exitflag,output] = lsqnonlin(functionToMinimize,x0,lb,ub,options);
else
    t = contains(solveFor,'fSCA','IgnoreCase',true);
    if nnz(t)
        Aeq = zeros(1,length(x0));
        Aeq(find(t):end) = 1;
        beq = 1;
    else
        Aeq = [];
        beq = [];
    end
    A = [];
    b = [];
    [x,fval,exitflag,output] =...
        fmincon(functionToMinimize,x0,A,b,Aeq,beq,lb,ub,[],options);
end
stats.x = x;
switch method
    case 'leastSquares'
        stats.normResiduals = sqrt(fval);
        oStruct.solver = 'lsqnonlin';
    case 'normResiduals'
        stats.normResiduals = fval;
        oStruct.solver = 'fmincon';
    case 'spectralAngle'
        stats.sineAngle = fval;
        oStruct.solver = 'fmincon';
end
% calculate regression between measured and modeled
warnID = 'curvefit:fit:iterationLimitReached';
warnStruct = warning('off',warnID);
[F,G,~] = fit(mR, measuredReflectance,'poly1','robust','bisquare','weights',passWeight);
warning(warnStruct);
ci = confint(F);
oStruct.method = method;
stats.exitflag = exitflag;
stats.slope = F.p1;
stats.confidenceInterval = [ci(1,1) ci(2,1)]';
stats.rmse = G.rmse;
stats.goodness = G.adjrsquare;
stats.output = output;
stats.options = options;
% not retrieving or saving lambda, grad, or hessian
if exitflag<=0
    warning('%s, %s, %s',oStruct.solver,oStruct.method,output.message)
end

% put solution into output structure
for k=1:length(solveFor)
    if contains(solveFor{k},'1') || contains(solveFor{k},'2')
        sn = solveFor{k};
        nv = str2double(sn(end));
        fn = sn(1:end-1);
        oStruct.(fn)(nv) = x(k);
    else
        oStruct.(solveFor{k}) = x(k);
    end
end

if nargout>1
    varargout{1} = stats;
    if nargout>2
        varargout{2} = setargc(prescription,x,solveFor,passSolver);
    end
end

    function D = SnowCloudDiff(x)
        % vector of differences between measured and modeled snow or cloud reflectance
        newScrip = setargc(prescription,x,solveFor,passSolver);
        mR =  SPIReS_fwd(newScrip);
        % vector of signed difference between model and measure
        D = passWeight.*(mR-measuredReflectance);
        prescription = newScrip;
    end

    function N = SnowCloudNormResiduals(x)
        % norm of residuals between measured and modeled snow or cloud reflectance
        D = SnowCloudDiff(x);
        % norm of difference between model and measure
        N = norm(D);
    end

    function S = SnowCloudSpectralAngle(x)
        % norm of residuals between measured and modeled snow or cloud reflectance
        newScrip = setargc(prescription,x,solveFor,passSolver);
        mR =  SPIReS_fwd(newScrip);
        % account for angle between model and measure
        cosA = dot(passWeight.*mR,passWeight.*measuredReflectance)/...
            (norm(passWeight.*mR)*norm(passWeight.*measuredReflectance));
        % want the sine of the angle (because 0 means close to similarity)
        S = sqrt((1-cosA)*(1+cosA));
        prescription = newScrip;
    end
end

function newScrip=setargc(fscript,z,solveName,solver)
% new prescription must combine multiple endmembers into
% fractional coverage vector, so length may be shorter than solveName
subStruct = char(fscript.substance);
for m=1:length(solveName)
    if strcmpi(solveName{m},'fSCA')
        if strcmpi(solver,'leastSquares')
            fscript.snow.fSCA = [z(m) 1-z(m)];
        else
            fscript.snow.fSCA = z(m:end);
        end
    elseif contains(solveName{m},'muS','IgnoreCase',true)
        fscript.Illumination.muS = z(m);
    elseif contains(solveName{m},'1') || contains(solveName{m},'2')
        sn = solveName{m};
        nv = str2double(sn(end));
        fn = sn(1:end-1);
        fscript.snow.(fn)(nv) = z(m);
    else
        fscript.(subStruct).(solveName{m}) = z(m);
    end
end
newScrip = fscript;
end