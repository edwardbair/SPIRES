function [ W ] = findHimawariW( file )
% [ W ] = himawariW( file )
%find FWHM for Himawari-8 or -9 AHI sensor

bands = (1:16)';
W = zeros(length(bands),2);

for k=1:length(bands)
    Y = xlsread(file,['Band ' num2str(k)]);
    disp(['band ' num2str(k) ' read, size matrix = ' num2str(size(Y))])
    G = griddedInterpolant(Y(:,1),Y(:,3)-max(Y(:,3))/2,'pchip','nearest');
    fwhm = @(x) G(x);
    xlow = fzero(fwhm,min(Y(:,1)));
    xhigh = fzero(fwhm,max(Y(:,1)));
    W(k,:) = [xlow xhigh];
end