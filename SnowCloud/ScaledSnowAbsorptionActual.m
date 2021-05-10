function outS = ScaledSnowAbsorptionActual(w,waveUnits,R,varargin)
%outS = ScaledSnowAbsorptionActual(w,waveUnits,R,[,...])
%ScaledSnowAbsorptionActual - relative difference and integral across
%absorption bands
%Input
%   w - wavelength
%   waveUnits - units of wavelength
%   R - reflectance for values of w
%Optional input, name-value pairs, needed to cover vis/nir dust/soot
%absorption
%   'snowRadius' - effective radius of snow
%   'radiusUnits' - default 'um'
%   'cosZ' - cosine of illumination angle
%   'lookup' - use Mie lookup tables, default true
%
%Output structure, 2 or 3 substructures corresponding to dust/soot
%(if snowRadius & cosZ specified), nm1030, and nm1240, each containing
%   name - name of feature
%   wavelength across feature
%   reflDiff - difference continuum minus actual, scaled to continuum
%       values, corresponding to wavelengths
%   continuumR - continuum reflectance, corresponding to wavelengths
%   actualR - actual reflectance, corresponding to wavelengths
%   integralDiff - integral of reflDiff

p = inputParser;
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
nonnegativeFcn = @(x) isnumeric(x) && all(x(:)>=0);
rangeFcn = @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1;
addRequired(p,'w',positiveFcn);
addRequired(p,'waveUnits',@ischar);
addRequired(p,'R',nonnegativeFcn);
addParameter(p,'snowradius',[],positiveFcn)
addParameter(p,'cosz',[],rangeFcn)
addParameter(p,'radiusunits','um',@ischar)
addParameter(p,'lookup',true,@islogical)
parse(p,w,waveUnits,R,varargin{:});

dirtySnow = ~isempty(p.Results.snowradius) && ~isempty(p.Results.cosz);
if dirtySnow
    N = 3;
else
    N = 2;
end
for k=1:N
    switch k
        case 1
            nm1030.name = '1030 nm absorption';
            nm1030.wavelength = (930:5:1120)';
            tx{1} = nm1030.wavelength<940 | nm1030.wavelength>1100;
        case 2
            nm1240.name = '1240 nm absorption';
            nm1240.wavelength = (1140:5:1370)';
            tx{2} = nm1240.wavelength<1150 | nm1240.wavelength>1350;
        case 3
            dustsoot.name = '400-800 nm dust/soot absorption';
            dustsoot.wavelength = (400:5:800)';
            %             tx(:,3) = false(length(dustsoot.wavelength),1);
    end
end
% wavelength range over the absorption bands
allWave = [400 1370];
w = convertUnits(w,p.Results.waveUnits,'nm');

%check that wavelengths cover the absorption bands
assert(min(w)<=min(allWave) && max(w)>=max(allWave),...
    'input wavelengths must cover range 400 to 1370 nm')

for k=1:N
    switch k
        case 1
            tempS = nm1030;
        case 2
            tempS = nm1240;
        case 3
            tempS = dustsoot;
    end
    % continuum reflectance
    if k==1 || k==2
        Rband = pchip(w,R,tempS.wavelength);
        rc = pchip(tempS.wavelength(tx{k}),Rband(tx{k}),tempS.wavelength);
        % scaled difference
        tempS.reflDiff = (rc-Rband)./rc;
        tempS.continuumR = rc;
        tempS.actualR = Rband;
    else % dirty snow so 'continuum' is clean snow
        % scale R to match reflectance in the 1030 and 1240 absorption regions
        wBand = [nm1030.wavelength; nm1240.wavelength];
        Rclean = SnowCloudSpectralReflectance(p.Results.cosz,p.Results.snowradius,...
            p.Results.radiusunits,w,'nm','lookup',p.Results.lookup);
        modelR = pchip(w,Rclean,wBand);
        measR = pchip(w,R,wBand);
        FB = fit(measR,modelR,'poly1','Robust','bisquare');
        fixR = FB(pchip(w,R,dustsoot.wavelength));
        % clean snow in dust/soot wavelengths
        rc = pchip(w,Rclean,dustsoot.wavelength);
        tempS.reflDiff = (rc-fixR)./rc;
        tempS.continuumR = rc;
        tempS.actualR = fixR;
    end
    tempS.integralDiff = trapz(tempS.wavelength,tempS.reflDiff);
    
    %results into structure
    switch k
        case 1
            nm1030 = tempS;
        case 2
            nm1240 = tempS;
        case 3
            dustsoot = tempS;
    end
end

outS.nm1030 = nm1030;
outS.nm1240 = nm1240;
if dirtySnow
    outS.dustsoot = dustsoot;
end

end