function MainDemo()
% Live function to reproduce Figures 1 through 5 in the paper Revisiting the
% topographic horizon problem in the era of big data and parallel computing.
% submitted to IEEE Geoscience and Remote Sensing Letters. The demo uses a
% 0.25x0.25 degree subset of the 1x1 degree data, so the input file is
% smaller. The user can also run the demo with a different DEM.
%
% Note that the profile in Figure 3 was derived from the full 1x1 degree
% image and is a separate file.

% Change the next lines if needed to identify the folder and file with the demo data
folder = '..\demodata';
file = 'SmallDEM.tif';
profileFile = 'profileDemo.mat';

% Change the next line to choose an appropriate datetime for the sun angle
% calculations in Fig. 4, which illustrates the effect of slopes and horizons
% on topographic shading. If the time zone is omitted, UTC is assumed. A
% MATLAB datenum can also be used, with the time zone UTC.
dt = datetime('2020-12-21 09:00','TimeZone','Asia/Calcutta');

% read the elevation data
if verLessThan('map','5.0')
        [Z,R] = geotiffread(fullfile(folder,file));
else
        [Z,R] = readgeoraster(fullfile(folder,file));
end
% demo is for a lat-lon dataset, could be revised to use a projected DEM
assert(contains(class(R),'geographic','IgnoreCase',true),...
    'input file must be geographic (code could be revised for projected DEM)')

% print the overall message
fprintf('This Demo reproduces Figures 1 to 5 for the paper Revisiting the\ntopographic horizon problem in the era of big data and parallel\ncomputing. submitted to IEEE Geoscience and Remote Sensing Letters.\n\n')

% data for Demo_profile are in a .mat file
m = matfile(fullfile(folder,profileFile));

elapsed = Demo_imageDEM(Z,R);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

elapsed = Demo_profile(m);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

elapsed = Demo_rotate(Z);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

% if you want to run without parallelism, just set useParallel = false
delete(gcp('nocreate'));
poolobj = parpool;
useParallel = poolobj.NumWorkers>1;
elapsed = Demo_shade(Z,R,dt,useParallel);
if useParallel
    fprintf('Elapsed time = %f seconds using %d cores\n\n',elapsed,poolobj.NumWorkers);
else
    fprintf('Elapsed time = %f seconds\n\n',elapsed);
end

elapsed = Demo_viewFactor(Z,R,useParallel);
if useParallel
    fprintf('Elapsed time = %f seconds using %d cores\n\n',elapsed,poolobj.NumWorkers);
else
    fprintf('Elapsed time = %f seconds\n\n',elapsed);
end
end