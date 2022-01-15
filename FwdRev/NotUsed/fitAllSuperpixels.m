function [RS,LAP] = fitAllSuperpixels(S,Ftrans,Fdf,method,allUnknowns)
% [RS,LAP] = fitAllSuperpixels(S,Ftrans,Fdf,method,allUnknowns)
%run fitSuperpixel on all
narginchk(5,5)
nargoutchk(0,3)
% allUnknowns ={'radius','fSCA','muS','dust','dustRadius'};
possibleVariables = {'radius','fSCA','muS','LAPfraction','LAPradius'};
% nSig = 3; % # of significant digits in error term
if contains(method,'spec')
    errorTerm = 'sineAngle';
else
    errorTerm = 'normResiduals';
end

% arrays to hold all results
nVals = size(S.spectralSP,1);
idx = zeros(nVals,1);
regressionSlope = zeros(size(idx));
confidenceInterval = zeros(length(idx),2);
goodness = zeros(size(regressionSlope));
% fSCA
if any(contains(allUnknowns,'fSCA','IgnoreCase',true))
    % extra column if shade
    if any(contains(allUnknowns,'shade','IgnoreCase',true))
        fSCA = zeros(nVals,3);
    else
        fSCA = zeros(nVals,2);
    end
    fSCA(:,1) = 1;
end
if any(contains(allUnknowns,'dust','IgnoreCase',true))
    cleanSnow = false;
    LAP = 'dust';
elseif any(contains(allUnknowns,'soot','IgnoreCase',true))
    cleanSnow = false;
    LAP = 'soot';
else
    cleanSnow = true;
    LAP = '';
end
if ~cleanSnow
    LAPfraction = zeros(nVals,1);
    if any(contains(allUnknowns,'dustRadius','IgnoreCase',true)) ||...
            any(contains(allUnknowns,'sootRadius','IgnoreCase',true))
        LAPradius = zeros(nVals,1);
    end
end
if any(contains(allUnknowns,'radius'))
    radius = nan(nVals,1);
end
if any(contains(allUnknowns,'mus','IgnoreCase',true))
    muS = nan(nVals,1);
end
errorValue = nan(nVals,1);
funcCount = zeros(nVals,1);
passVals = {allUnknowns,S.wavelength,S.backReflectance,...
    Ftrans,Fdf,method,errorTerm,possibleVariables};
W = parallel.pool.Constant(passVals);
spectralSP = S.spectralSP;
topography = S.topography(:,1:4);
parfor k=1:size(spectralSP,1)
    pv = W.Value;
    % pv{1-6}: unknowns, wavelength, backRefl, Ftrans, Fdf, method
    [ostruct,stats,LP] = fitSuperpixel(pv{1},...
        pv{2},spectralSP(k,:),topography(k,:),pv{3},...
        'atmos',pv{4},'diffuse',pv{5},'method',pv{6});
    goodExit = stats.exitflag>0;
    if goodExit
        idx(k) = k;
        errorValue(k) = stats.(pv{7}); % errorTerm
        regressionSlope(k) = stats.slope;
        confidenceInterval(k,:) = stats.confidenceInterval;
        goodness(k) = stats.goodness;
        fn = fieldnames(ostruct);
        t = contains(fn,'solver') | contains(fn,'method');
        fn = fn(~t);
        possV = pv{8};
        for m=1:length(fn)
            if any(contains(possV,fn{m}))
                switch fn{m}
                    case 'radius'
                        radius(k) = LP.snow.radius;
                    case 'muS'
                        muS(k) = LP.Illumination.muS;
                    case 'fSCA'
                        fSCA(k,:) = LP.snow.fSCA;
                    case 'LAPfraction'
                        LAPfraction(k) = LP.snow.LAPfraction;
                    case 'LAPradius'
                        LAPradius(k) = LP.snow.LAPradius;
                end
            else
                warning('output %s not among possible variables',fn{m})
            end
        end
    end
end
% consolidate the results
RS = struct;
RS.idx = idx;
if exist('radius','var')
    RS.radius = radius;
end
if exist('muS','var')
    RS.muS = muS;
end
if exist('fSCA','var')
    RS.fSCA = fSCA;
end
if strcmp(errorTerm,'normResiduals')
    RS.normResiduals = errorValue;
elseif strcmp(errorTerm,'sineAngle')
    RS.sineAngle = errorValue;
end
if exist('regressionSlope','var')
    RS.regressionSlope = regressionSlope;
end
if exist('confidenceInterval','var')
    RS.confidenceInterval = confidenceInterval;
end
if exist('goodness','var')
    RS.goodness = goodness;
end
if ~cleanSnow
    if exist('LAPfraction','var')
        RS.LAPfraction = LAPfraction;
    end
    if exist('LAPradius','var')
        RS.LAPradius = LAPradius;
    end
end
RS.funcCount = funcCount;
RS.topography = topography;
end