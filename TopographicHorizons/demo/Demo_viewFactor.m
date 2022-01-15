function et = Demo_viewFactor(Z,R,useParallel)
%figure showing view factor

tic; % start the timer

nHorz = 16;
if useParallel
    [A,H,~] = horizonAllDirections(double(Z),R,'nHorz',nHorz,'parallel','rotate');
else
    [A,H,~] = horizonAllDirections(double(Z),R,'nHorz',nHorz);
end
vf = viewFactor(A,H,double(Z),R);

fprintf('This code %s reproduces Fig. 5 but based on only %d horizon directions\n',...
    mfilename,nHorz);
ptile = 5;
ppart = prctile(vf(:),ptile);
fprintf('The %d-th percentile view factor is %.2f, so the image colors represent the range [%.2f 1.0]\n',...
    ptile, ppart, ppart);
fprintf('The white parts of the image represent values outside this range\n')

% image, shaded relief, shade by slope, shade by horizon
% shaded relief
figure('Name','Fig. 5 View Factor')
ax = setAxes(R,true); %#ok<NASGU>
vf(vf<ppart) = NaN;
geoshow(vf,R,'DisplayType','Surface')
colorbar('Location','SouthOutside')

et = toc;

end