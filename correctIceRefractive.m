function [waveIce, N] = correctIceRefractive(file,wvice,iceOrig)
% [waveIce, N] = correctIceRefractive(file,wvice,iceOrig)
%correct Warren-Brandt ice refractive index (imag part) with later
%measurements by Picard et al.
%
%Input
% file - with the Picard data, but in um
% wvice - log(wavelength), um
% iceOrig - complex index (real col 1, imag col 2)
%
%Output
% waveIce - log(consolidated wavelengths), um
% N - consolidated new refractive index (col 1 real, col 2 imaginary)

M = matfile(file); % file with the Picard data, convered to um

% insert the Picard et al. wavelengths
wvPicard = log(M.wvPicard);
t = wvice>=wvPicard(1) & wvice<=wvPicard(end);
lowend = ~t & wvice<wvPicard(1);
highend = ~t & wvice>wvPicard(end);
newWave = unique([wvice(~t); wvPicard]);

% real part
F = fit(wvice,iceOrig(:,1),'pchip');
realpart = F(newWave);

% imag part
newImag = cat(1,iceOrig(lowend,2),M.kicePicard,iceOrig(highend,2));
N = [realpart newImag];

waveIce = newWave;

end

