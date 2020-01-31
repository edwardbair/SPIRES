function [outStruct] = invertCUES(cosZ,spectrum,specDirect,specDiff,waveBands,albedo)
% [outStruct] = invertCUES(cosZ,spectrum,specDirect,specDiff,waveBands,albedo)
%quick & dirty function to estimate snow properties from measurements at
%CUES, will be modified to be more general
%
%Input
%   cosZ - vector of length N of cosine zenith angles at which measurements
%       were made
%   spectrum - #L wavelengths, nm, of spectral distribution of radiation, for
%       example from SMARTS295
%   specDirect - direct irradiances, e.g. from SMARTS295, matrix of size LxN,
%       with irradiances along the columns, with one column for each cosZ
%   specDiff - diffuse irradiances, again of size LxN
%   waveBands - matrix of size Mx2, whereby we have M radiometers with
%       wavelength ranges with the lower wavelength in column 1 and the
%       upper wavelengths in column 2, in nm
%   albedo - matrix of size NxM with albedo measurements, one row for each
%       cosZ, albedo for waveband 1 in column 1, for waveband 2 in column
%       2, and so forth
%
%Output
%   structure with best estimates of snow properties
%       if 2 measurements are specified, solve for grain radius and dust
%       content
%       if 3 or more measurements are specified, solve for the same 2
%       variables but also solve for the dust radius

S = SnowCloudLimits;
% unknowns for solving for grain size and dust content
xTwo0 = [S.defaultSnowRadius mean(S.dust)];
lbTwo = [S.snowRadius(1) 0];
ubTwo = [S.snowRadius(2) S.dust(2)];
% unknowns for solving for grain size, dust content, dust radius
xThree0 = cat(2,xTwo0,S.defaultDustRadius);
lbThree = cat(2,lbTwo,S.dustRadius(1));
ubThree = cat(2,ubTwo,S.dustRadius(2));

% variables to pass to functions that calculate difference between model
% and measurement
passCos = cosZ;
passWave = spectrum;
passDirect = specDirect;
passDiffuse = specDiff;
passBands = waveBands;
passMeas = albedo;

% solve for grain size and dust
[x,resnorm,residual,exitflag,output] = lsqnonlin(@reflDiff,...
    xTwo0,lbTwo,ubTwo);
X2.radius = x(1);
X2.dust = x(2);
X2.resnorm = resnorm;
X2.residual = residual;
X2.exitflag = exitflag;
X2.output = output;

%solve for grain size, dust, and dust radius
if numel(albedo)>2
    [x,resnorm,residual,exitflag,output] = lsqnonlin(@reflDiff,...
        xThree0,lbThree,ubThree);
    X3.radius = x(1);
    X3.dust = x(2);
    X3.dustRadius = x(3);
    X3.resnorm = resnorm;
    X3.residual = residual;
    X3.exitflag = exitflag;
    X3.output = output;
end

% output structure
outStruct.RadiusDust = X2;
if exist('X3','var')
    outStruct.RadiusDustDustR = X3;
end

    function Diff = reflDiff(x)
        dval = zeros(length(passCos),size(passBands,1));
        whichFun = length(x);
        for k=1:length(passCos)
            switch whichFun
                case 2
                    R = SnowCloudIntgRefl(passWave,'nm',...
                        [passDirect(:,k) passDiffuse(:,k)],...
                        'snow','radius',x(1),'dust',x(2),...
                        'cosZ',passCos(k),'wavelength',passBands,'waveU','nm');
                case 3
                    R = SnowCloudIntgRefl(passWave,'nm',...
                        [passDirect(:,k) passDiffuse(:,k)],...
                        'snow','radius',x(1),'dust',x(2),'dustR',x(3),...
                        'cosZ',passCos(k),'wavelength',passBands,'waveU','nm');
            end
            for b=1:size(passBands,1)
                dval(k,b) = passMeas(k,b)-R.reflectance(b);
            end
        end
        Diff = dval(:);
    end
end