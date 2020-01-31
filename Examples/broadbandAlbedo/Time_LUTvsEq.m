function [timeStruct] = Time_LUTvsEq(totalNum)
% [timeStruct] = Time_LUTvsEq(totalNum)
%time the LUT approach vs that using power function and rational equation
%
%Input
%   totalNum - number of combined values of cosZ, radius, dust, altit
%       (or for power function, just the cosZ and radius are used)

%values for the LUT
n = power(totalNum,.2);
nDust = ceil(n);
nAlt = floor(n);
remain = totalNum/(nDust*nAlt);
n = sqrt(remain);
nRad = ceil(n);
nCos = floor(n);

%values for power + rational
n = sqrt(totalNum);
pRad = ceil(n);
pCos = floor(nRad*nCos*nDust*nAlt/pRad);

%generate random values for input
% gaussian mixture for radii, linear for others
tic
gm = gmdistribution(250,100);
[rad,dust,alt,cosZ] = ndgrid(random(gm,nRad),1e-3*rand(nDust,1),...
    (4.5-1.5)*rand(nAlt,1)+1.5,rand(nCos,1));
[radEq,mu0] = ndgrid(random(gm,pRad),rand(pCos,1));
elapsedTime = toc;
fprintf('for LUT, arrays of size %g created\nfor Eq, arrays of size %g created\ntime = %f sec\n',...
    numel(rad),numel(radEq),elapsedTime);
timeStruct.LUTsize = numel(rad);
timeStruct.EqSize = numel(radEq);

% coefficients for equation (from Bair et al Table 2)
P = [-9.025001 -6.853901 -6.360441;...
    0.05785986 0.273218 0.1890732;...
    0.07632736 1.017243 0.4149719];
Q = [1 92.35081 27.87415;...
    1 1.28665 1.53981;...
    0 1 0.3373872];

% time for the LUT solution
tic
refl = AlbedoLookup(rad,cosZ,[],alt,'dust',dust);
elapsedTime = toc;
fprintf('time to run the lookup table = %f seconds\n',elapsedTime)
timeStruct.LUT_time = elapsedTime;
timeStruct.reflLUT = refl;

% time for the equation
tic
mu0sq = mu0.^2;
a = (P(1,1)*mu0sq+P(1,2)*mu0+P(1,3))./(Q(1,1)*mu0sq+Q(1,2)*mu0+Q(1,3));
b = (P(2,1)*mu0sq+P(2,2)*mu0+P(2,3))./(Q(2,1)*mu0sq+Q(2,2)*mu0+Q(2,3));
c = (P(3,1)*mu0sq+P(3,2)*mu0+P(3,3))./(Q(3,1)*mu0sq+Q(3,2)*mu0+Q(3,3));
refl = a.*radEq.^b+c;
tz = mu0>=0.05;
mumin = min(mu0(tz));
nfill = find(mu0==mumin,1,'first');
tsmall = mu0<mumin;
refl(tsmall) = refl(nfill);
elapsedTime = toc;
fprintf('time to run the equation (power+rational) = %f seconds\n',elapsedTime)
timeStruct.EqTime = elapsedTime;
timeStruct.reflEq = refl;

% time for the LUT solution with clean snow, one elevation
tic
refl = AlbedoLookup(rad,cosZ,[],alt);
elapsedTime = toc;
fprintf('time to run the LUT with clean snow = %f seconds\n',elapsedTime)
timeStruct.LUTclean_time = elapsedTime;
timeStruct.reflLUTclean = refl;

end