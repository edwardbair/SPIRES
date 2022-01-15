function [wt,varargout]=spectralSnowWeight(Prescription,Ftrans)
% [wt [,snowDiff,backDiff,solarTrans]]=spectralSnowWeight(Prescription)
% wt=spectralWeight(Prescription)
% weights based on differences across grain size and contaminant amounts for snow
% and between snow reflectance and the background
%
%Input
% Prescription - from fwdPrescription (but this routine will set some values)
% Ftrans - function to calculate atmospheric transmission, squared
%
%Output
% wt - spectral weights for use in inversion
%Optional output - components of the weight
% snowDiff - normalized difference between fine, clean and coarse, dirty snow
% backDiff - normalized difference between snow and background
% solarTrans - solar radiation transmittance


narginchk(2,2)
nargoutchk(0,4)

p = inputParser;
addRequired(p,'Prescription',@isstruct);
addRequired(p,'Ftrans',@(x) strcmpi(class(x),'cfit') ||...
    contains(class(x),'interpolant','IgnoreCase',true))
parse(p,Prescription,Ftrans)

Prescription.Illumination.muS = Prescription.Illumination.cosZ;
Prescription.snow.fSCA = 1;
Prescription.snow.fOther = 0;
S = SnowCloudLimits;
fineSnowRadius = 100;
coarseSnowRadius = 1000;
fractionMaxLAP = 0.5;
if strcmpi(Prescription.snow.LAP,'dust')
    Prescription.snow.LAPfraction = max(S.dust)*fractionMaxLAP;
else
    Prescription.snow.LAPfraction = max(S.soot)*fractionMaxLAP;
end
atmosWeight = getAtmosProp(Ftrans,Prescription.Spectrum.wavelength,...
    Prescription.Illumination.elevation,Prescription.Illumination.cosZ);

%clean snow, fine and coarse
Prescription.snow.cleanSnow = true;
Prescription.snow.radius = fineSnowRadius;
R1 =  SPIReS_fwd(Prescription);
Prescription.snow.radius = coarseSnowRadius;
R2 =  SPIReS_fwd(Prescription);

%dirty snow, coarse
Prescription.snow.cleanSnow = false;
R3 =  SPIReS_fwd(Prescription);

% snow vs background
F = Prescription.Substrate.fR0{1};
backRefl = F(Prescription.Spectrum.wavelength);
backDiff = rescale((R1-backRefl)+(R3-backRefl));

% weights for the spectral reflectance
snowDiff = max([rescale(R1-R2) rescale(R2-R3)],[],2);

% weight is max of backDiff, snowDiff times atmos trans
wt = rescale(max([snowDiff(:) backDiff(:)],[],2).*atmosWeight(:),.01,.99);

%optional output
if nargout>1
    for k=2:nargout
        n = k-1;
        switch n
            case 1
                varargout{n} = snowDiff; %#ok<*AGROW>
            case 2
                varargout{n} = backDiff;
            case 3
                varargout{n} = atmosWeight;
        end
    end
end
end