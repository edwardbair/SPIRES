function et = Demo_profile(m)
%figure to demo the profile across the diagonal of the Himachal Pradesh DEM

fprintf('This code %s reproduces Fig. 2\n',mfilename)
z = m.z;

tic; % start the timer

% plot profile
dkm = m.d/1000;
figure('Name','Fig. 2 Elevation/Horizon profile')
plot(dkm,z,'k','LineWidth',1)
xlabel('distance (km)')
ylabel('elevation (m)')
hold on;

% horizons
[horzAng,horzDis,horzAzm,horzPts] = horizonAlongProfile(m.lat,m.lon,z,m.ellipsoid); %#ok<ASGLU>
hpu = unique(horzPts); % unique horizon points
fprintf('%d horizon pts, angles ranging from %.1f to %.1f degrees\n',...
    length(hpu),min(horzAng),max(horzAng))
disp([num2str(length(hpu)) ' horizon pts']);
% number of points to each horizon pt
npts = zeros(size(hpu));
for k=1:length(npts)
    npts(k) = nnz(horzPts==hpu(k));
end

% find the horizon point at the highest elevation
khigh = find(z==max(z(hpu)));
% find the horizon point beyond khigh that serves the most points
knext = hpu(npts==max(npts(hpu>khigh)));
kplot = [khigh knext];

% those above 97th percentile in terms of number of points served
t97 = npts>prctile(npts,97);
hp97 = hpu(t97);
disp([num2str(length(hp97)) ' horizon pts shown'])
disp(['forward direction azimuth ' num2str(horzAzm)...
    ' counterclockwise from South'])

% plot the paths for the horizons at kplot
for m=1:2
    thisK = kplot(m);
    t = horzPts==thisK;
    hplot = find(t);
    for k=1:length(hplot)
        plot([dkm(hplot(k)) dkm(thisK)],[z(hplot(k)) z(thisK)],...
            'Color',[.8 .8 .8])
    end
end

% plot the 97th percentile points
scatter(dkm(hp97),z(hp97),200,'r.')

pbaspect([1.6 1 1])
hold off;

et = toc;

end
