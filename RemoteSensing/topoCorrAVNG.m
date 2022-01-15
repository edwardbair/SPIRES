function [outputfile] = topoCorrAVNG(avngTemplate,baseFolder)
% [outputfile] = topoCorrAVNG(avngTemplate,baseFolder)
%corrects AVIRIS-NG reflectance file for topography

reflFile = fullfile(baseFolder,'reflectance',...
    [avngTemplate '_rfl_v2m3'],[avngTemplate '_corr_v2m3_img.mat']);
obsFile = fullfile(baseFolder,'radiance',...
    [avngTemplate '_rdn_v2m3'],[avngTemplate '_rdn_v2m3_obs_ort.mat']);
R = load(reflFile);
O = load(obsFile);

% geometry from obs file
slope = O.obs(:,:,7);
aspect = 180-O.obs(:,:,8);
mu0 = cosd(O.obs(:,:,5));
phi0 = 180-O.obs(:,:,4);

% solar angle on slope
mu = sunslope(mu0,phi0,slope,aspect);
corrFactor = mu0./mu; % inverse of mu./mu0;
corrFactor(isinf(corrFactor)) = NaN;

% corrected values
refl = bip2bsq(R.reflectance);
for b=1:size(refl,3)
    X = squeeze(refl(:,:,b)).*corrFactor;
    refl(:,:,b) = X;
end
R.reflectance = bsq2bip(refl);
R.description = 'AVIRIS-NG atmospherically and topographically corrected reflectance';
R.corrFactor = corrFactor;
outputfile = fullfile(baseFolder,'reflectance',...
    [avngTemplate '_rfl_v2m3'],[avngTemplate '_topocorr_v2m3_img.mat']);
save(outputfile,'-struct','R','-v7.3')
end