function [oStruct,varargout] = invertSnowCloudIntgRefl(solarTbl,...
    reflectance,unknowns,varargin )
% [oStruct] = invertSnowCloudIntgRefl(solarTbl,reflectance,unknowns,prescription,prop/val )
% [oStruct] = invertSnowCloudIntgRefl(__________, substance,prescription,prop/val )
% [oStruct,stats] = invertSnowCloudIntgRefl(__________)
% [ostruct,stats,P] = invertSnowCloudIntgRefl(__________)
%
%solves for snow or cloud properties based on band-integrated reflectance
%   (i.e. the inverse of SnowCloudIntgRefl
%
%Required inputs, for the spectrum of the incoming irradiance
%   solarTbl - columns wavelength and irradiance, irradiance can have 2
%   columns, direct and diffuse
%   Make sure solarTbl.Properties.VariableUnits is set, e.g. = {'nm','W/m2/nm'}
%Required inputs, for the inversion
%   reflectance - measured reflectances at band passes, as determined by
%       'bandPass' (as specified below) or 'sensor'
%   unknowns - cell vector of snow or cloud properties to solve for, with
%       possibilities depending on the input substance and the number of
%       reflectance values available (reasonable abbreviations work)
%       for 'snow' - 'radius' or 'ssa', 'fSCA', 'wetness', 'waterEquivalent',
%           'dust', 'dustRadius', 'soot', 'sootRadius',
%           'corrFactor' (cosS * muS)
%       for 'iceCloud' - 'radius', 'waterEquivalent', 'dust',
%           'dustRadius', 'soot', 'sootRadius'
%       for 'waterCloud' - 'radius', 'waterEquivalent', 'dust',
%           'dustRadius', 'soot', 'sootRadius'
%       for 'mixedCloud' - 'radius', 'waterRadius', 'wetness', 'waterEquivalent',
%           'dust','dustRadius', 'soot', 'sootRadius'
%
%optional arguments all name-value pairs, except substance and prescription
%arguments, if present, must come before the name-value pairs if either or
%both are specified
%   'substance', either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%       (but any unambiguous abbreviation beginning with first letter works)
%   prescription, to use and possibly modify an existing prescription
%Other inputs, name-value pairs (either 'sensor' or 'wavelength' must be specified)
%   'sensor' - any multispectral sensor known to SensorTable or
%   'bandPass' - Nx2 matrix of wavelengths, where equal to length(reflectance)
%   'waterEquivalent' - if not an unknown to solve for
%   'R0' - reflectance of underlying substrate if 'waterEquivalent' is specified
%       or of other endmember if 'fractionalCoverage' or 'fSCA' is
%       specified, must be scalar, length of 'wavelength', or most commonly
%       is a table with columns wavelength and reflectance
%   'cosZ' - cosine of illumination angle on flat surface, scalar (or
%       vector if using measurements at more than one time of day)
%   'muS' - cosine of illumination angle on slope, same size as cosZ
%   'cosS' - cosine of slope, same size as cosA
%   'sizeUnit' - units for sizes of scatterer and impurities, default 'mum'
%   'waveUnit' - units for wavelength, default 'nm'
%   'radius' - effective optical radius of scatterer, scalar
%   'ssa' - surface-specific area, kg/m^2, as an alternative to 'radius'
%   'solutionMethod' - to solve the inversion, either 'lsqnonlin' (default)
%       or 'spectralAngle'
%
%Output
%   oStruct - snow or cloud properties, depending on inputs
%Optional output
%   stats - statistics about solution
%   P - prescription at last call to snow reflectance function

narginchk(3,Inf)
nargoutchk(0,3)
passP = struct([]);
passTbl = solarTbl;

% parse inputs (pass to SetSnowCloud)
fscript = SetSnowCloud(mfilename,varargin{:});
reflectance = double(reflectance(:));
assert(length(reflectance)==size(fscript.bandPass,1),...
    'length of reflectance vector must equal length of bandPass matrix')

% parse unknowns
assert(length(unknowns)<=length(reflectance),...
    'length of reflectance vector must be >= number of unknowns')
validStrings = {'fSCA','wetness','dust','dustRadius',...
    'soot','sootRadius','corrFactor','waterRadius','radius',...
    'waterEquivalent','fractionalCoverage','ssa'};
equivStrings = {'fSCA','fractionalCoverage'};

solveFor = cell(length(unknowns),1);
for k=1:length(unknowns)
    matchedstr = validatestring(unknowns{k},validStrings);
    if contains(matchedstr,equivStrings(:),'IgnoreCase',true)
        if contains(matchedstr,equivStrings(:),'IgnoreCase',true)
            solveFor{k} = equivStrings{1};
        end
    else
        solveFor{k} = matchedstr;
    end
end

% make sure wavelenghts long enough to get size of scatterer, if asked for
if any(contains(solveFor,'iceRadius','IgnoreCase',true)) ||...
        any(contains(solveFor,'ssa','IgnoreCase',true)) ||...
        any(contains(solveFor,'waterRadius','IgnoreCase',true))
    assert(convertLengthUnits(max(fscript.bandPass(:)),fscript.waveUnit,'nm')>=1060,...
        'maximum bandPass must be >= %f %s to retrieve size of snow or cloud scatterer',...
        convertLengthUnits(1060,'nm',fscript.waveUnit),fscript.waveUnit)
end
% make sure wavelengths short enough to get dust or soot characterization
if any(contains(solveFor,'dust','IgnoreCase',true)) ||...
        any(contains(solveFor,'soot','IgnoreCase',true))
    assert(convertLengthUnits(min(fscript.bandPass(:)),fscript.waveUnit,'nm')<=700,...
        'minimum bandPass must be <= %f %s to retrieve dust or soot properties',...
        convertLengthUnits(700,'nm',fscript.waveUnit),fscript.waveUnit)
end

% intial values and limits
[x0,lb,ub,extraEndMember] = setBounds(solveFor,fscript);

% solving method depends on input
passWeight = 1;
useParallel=true;
switch fscript.solutionMethod
    % inversion method lsqnonlin uses the signed differences between measurement and model
    case 'lsqnonlin'
        options = optimset('Display','off','FinDiffType','central','UseParallel',useParallel);
        [x,resnorm,residual,exitflag,output,lambda,jacobian] =...
            lsqnonlin(@SnowCloudDiff,x0,lb,ub,options); %#ok<ASGLU>
        stats.resnorm = resnorm;
        stats.residual = residual;
        stats.exitflag = exitflag;
        stats.output = output;
        %         stats.lambda = lambda;
        %         stats.jacobian = jacobian;
        if exitflag<=0
            warning('lsqnonlin: %s',output.message)
        end
    case 'spectralAngle'
        options = optimset('Display','off','FinDiffType','central','UseParallel',useParallel);
        [x,fval,exitflag,output,lambda,grad,hessian] =...
            fmincon(@SnowCloudSpectralAngle,x0,...
            [],[],[],[],lb,ub,[],options); %#ok<ASGLU>
        stats.spectralAngle = fval;
        stats.exitflag = exitflag;
        stats.output = output;
        %         stats.lambda = lambda;
        %         stats.grad = grad;
        %         stats.hessian = hessian;
        if exitflag<=0
            warning('fmincon: %s',output.message);
        end
    otherwise
        error('''solutionMethod'' ''%s'' not recognized',fscript.solutionMethod)
end

% put solution into output structure
for k=1:length(solveFor)
    oStruct.(solveFor{k}) = x(k);
    if strcmpi(solveFor{k},'fSCA')
        if size(fscript.R0.reflectance,2)>1
            oStruct.otherEndMem = [x(end) 1-x(end)-x(k)];
        else
            oStruct.otherEndMem = 1-x(k);
        end
    end
end

if nargout>1
    varargout{1} = stats;
    if nargout>2
        varargout{2} = passP;
    end
end

    function diffR = SnowCloudDiff(x)
        % difference between measured and model snow or cloud reflectance
        argc = cell(2*length(solveFor),1);
        v = 1;
        for m=1:length(solveFor)
            if strcmpi(solveFor{m},'fSCA')
                argc{v} = 'fractionalCoverage';
                if extraEndMember
                    argc{v+1} = [x(m) x(end) 1-x(m)-x(end)];
                else
                    argc{v+1} = [x(m) 1-x(m)];
                end
            else
                argc{v} = solveFor{m};
                argc{v+1} = x(m);
            end
            v = v+2;
        end
        % correct for topography either by recomputing effective modeled
        % reflectance (i.e. when viewed from above assuming flat surface)
        % or by recomputing reflectance on slope
        if any(contains(solveFor,'corrFactor','IgnoreCase',true))
            v = contains(solveFor,'corrFactor','IgnoreCase',true);
            m = find(v);
            cf = argc{2*m};
            doCorrFactor = true;
        elseif ~isempty(fscript.corrFactor)
            doCorrFactor = true;
            cf = fscript.corrFactor;
        else
            doCorrFactor = false;
        end
            [Rtbl,passP] = SnowCloudIntgRefl(passTbl,fscript,argc{:});
            modelRefl = Rtbl.reflectance;
        
        % correct for topography by recomputing measured reflectance, i.e.
        % if viewed normal to slope
        if doCorrFactor
            diffR = passWeight.*(reflectance.*fscript.cosZ./cf-modelRefl);
        else
            try
                diffR = passWeight.*(reflectance-modelRefl);
            catch
                pause
            end
        end
    end

    function ang = SnowCloudSpectralAngle(x)
        % spectral angle between model and actual reflectance
        
        argc = cell(2*length(solveFor),1);
        v = 1;
        for m=1:length(solveFor)
            if strcmpi(solveFor{m},'fSCA')
                argc{v} = 'fractionalCoverage';
                if extraEndMember
                    argc{v+1} = [x(m) x(end) 1-x(m)-x(end)];
                else
                    argc{v+1} = [x(m) 1-x(m)];
                end
            else
                argc{v} = solveFor{m};
                argc{v+1} = x(m);
            end
            v = v+2;
        end
        [Rtbl,passP] = SnowCloudIntgRefl(passTbl,fscript,argc{:});
        modelRefl = Rtbl.reflectance;
        
        %spectral angle
        cosAng = dot(passWeight.*modelRefl,passWeight.*reflectance)/...
            (norm(passWeight.*modelRefl)*norm(passWeight.*reflectance));
        ang = acosd(cosAng);
    end
end