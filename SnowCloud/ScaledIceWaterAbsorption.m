function outS = ScaledIceWaterAbsorption(w,waveUnits,R)
%outS = ScaledSnowAbsorptionActual(w,waveUnits,R,[,...])
%ScaledIceWaterAbsorption - absolute and relative differences and integrals
%between values and continuum across near-infrared absorption bands for ice
%and water
%
%Input
%   w - wavelength
%   waveUnits - units of wavelength
%   R - reflectance for values of w
%
%Output structure, 4 subs corresponding to 2 ice and 2 water absorption bands,
%identified by their central wavelengths, 1030 & 1270 for ice, 980 & 1200
%for water
%   name - name of feature
%   wavelength - across feature
%   reflDiff - difference continuum minus actual, scaled to continuum
%       values, corresponding to wavelengths
%   continuumR - continuum reflectance, corresponding to wavelengths
%   actualR - actual reflectance, corresponding to wavelengths
%   integralDiff - integral of reflDiff

p = inputParser;
positiveFcn = @(x) isnumeric(x) && all(x(:)>0);
nonnegativeFcn = @(x) isnumeric(x) && all(x(:)>=0);
addRequired(p,'w',positiveFcn);
addRequired(p,'waveUnits',@ischar);
addRequired(p,'R',nonnegativeFcn);
parse(p,w,waveUnits,R);

N = 4;
iceMin = [950 1115 1150 1330];
waterMin = [900 1070 1260];

% wavelength range over the absorption bands
allWave = [min(waterMin) max(iceMin)];
w = convertLengthUnits(w,p.Results.waveUnits,'nm');

%check that wavelengths cover the absorption bands
assert(min(w)<=min(allWave) && max(w)>=max(allWave),...
    'input wavelengths must cover range %d to %d nm',...
    min(waterMin),max(iceMin))

for k=1:N
    switch k
        case 1
            tempS.name = '1030 nm ice absorption';
            tempS.wavelength = (iceMin(1):5:iceMin(2))';
        case 2
            tempS.name = '1270 nm ice absorption';
            tempS.wavelength = (iceMin(3):5:iceMin(4))';
        case 3
            tempS.name = '980 nm water absorption';
            tempS.wavelength = (waterMin(1):5:waterMin(2))';
        case 4
            tempS.name = '1200 nm water absorption';
            tempS.wavelength = (waterMin(2):5:waterMin(3))';
    end
    % continuum reflectance
    Fband = fit(w,R,'pchipinterp');
    Rband = Fband(tempS.wavelength);
    Fc = fit([min(tempS.wavelength) max(tempS.wavelength)]',...
        Fband([min(tempS.wavelength) max(tempS.wavelength)]'),'poly1');
    rc = Fc(tempS.wavelength);
    % scaled difference
    tempS.reflDiff = (rc-Rband)./rc;
    tempS.continuumR = rc;
    tempS.actualR = Rband;
    tempS.integralDiff = trapz(tempS.wavelength,tempS.reflDiff);
    
    %results into structure
    switch k
        case 1
            outS.ice1030 = tempS;
        case 2
            outS.ice1270 = tempS;
        case 3
            outS.water980 = tempS;
        case 4
            outS.water1200 = tempS;
    end
end

end