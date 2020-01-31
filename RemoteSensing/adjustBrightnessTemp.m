function [ Tbfix ] = adjustBrightnessTemp(sensor,band,Tb,emiss,solar)
% [ Tbfix ] = adjustBrightnessTemp(sensor,band,Tb,emiss,solar)
%adjust brightness temperature to account for reflected solar radiation,
%needed to interpret daytime thermal signal in the 4 um region
%
%Input
%   sensor - e.g. 'modis'
%   band - numeric or char
%   (following arguments are scalars or vectors/matrices for each pixel)
%       Tb - measured brightness temperatures in band
%       emiss - emissivity, scalar averaged over each pixel for the band
%       solar - total (all wavelengths) solar irradiance on pixel(s)
%
%Output
%   Tbfix - estimated brightness temperatures

p = inputParser;
validationFcn = @(x) isnumeric(x) || ischar(x);
rangeFcn = @(x) isnumeric(x) && all(x(:)>=0 & x(:)<=1);
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
nonnegativeFcn = @(x) isnumeric(x) && all(x(:)>=0);
addRequired(p,'sensor',@ischar)
addRequired(p,'band',validationFcn)
addRequired(p,'Tb',positiveFcn)
addRequired(p,'emiss',rangeFcn)
addRequired(p,'solar',nonnegativeFcn)
parse(p,sensor,band,Tb,emiss,solar)

[Tb,emiss,solar] = checkSizes(p.Results.Tb,p.Results.emiss,p.Results.solar);

% band must be a cell vector (for SensorTable)
if isnumeric(p.Results.band)
    bandc = num2str(p.Results.band);
else
    bandc = {p.Results.band};
end

% sensor table entry for the band
X = SensorTable(p.Results.sensor);
t = strcmpi(X.Band,bandc);
X = X(t,:);

% correct each pixel
origSize = size(Tb);
Tb = Tb(:);
solar = solar(:);
emiss = emiss(:);
Tbfix = zeros(size(Tb));
lambda = [X.LowerWavelength X.UpperWavelength];
% wavelength in meters for calling Planck function
lambdaM = convertUnits(lambda,'um','m');
% reflected solar radiation, W/m^2, for each pixel
solarRad = solar.*(1-emiss)*...
    integral(@SolarScale,lambda(1),lambda(2));
% loop through the pixels
for k=1:numel(Tb)
    if solarRad(k)==0
        Tbfix(k) = Tb(k);
    else
        planckRad = PlanckIntegral(lambdaM(1),lambdaM(2),Tb(k));
        passRadiance = (planckRad-solarRad(k)/pi);
        passLambda = lambdaM;
        Tbfix(k) = fzero(@radDiff,Tb(k));
    end
end
Tbfix = reshape(Tbfix,origSize);

    function radD = radDiff(T)
        Rbright = PlanckIntegral(passLambda(1),passLambda(2),T);
        radD = Rbright-passRadiance;
    end

end