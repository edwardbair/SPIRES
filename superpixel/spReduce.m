% Draft function to implement jeff dozier's recomendations for using
% superpixels to reduce image size.
%
%make work with any superpixel function (SLIC or SLICSAM)
%
% How to identify the superpixels? Timbo can start with this. 
%The superpixels function in MATLAB can use a 3-band image.
%I think it’s better than using superpixels3.
%For the 3 bands, I suggest: 
%(i) norm of the reflectances; 
%(ii) spectral angle of the reflectances compared to some other spectrum, could be R0 or nominal snow pixel or . . .; 
%(iii) something else, could be NDSI or illumination angle. 
%I also suggest setting the ‘compactness’ value to a low number 
%and the number of superpixels perhaps 0.001 x image size or 0.0005 x image size.
%
% Characterize the superpixels by taking the median value of each band, and maybe include the NDSI or illumination angle if you want. 
%Make a table, either from a single image or a set of images.
%
% Reduce the number of superpixels as follows: 
%eliminate some of small size depending on a threshold, 
%convert each value to an ‘int8’ using my float2integer function, 
%which would set any NaN or out of range values to intmin(‘int8’). 
%Make a new table, eliminate any rows where any value is intmin(‘int8’), 
%and then run unique on the table to come up with the set of values to use 
%in building the LUT or the ML function.
%
% If desired (depending on how we develop the shortcut) 
%convert the values in the table back to floating point using the divisors 
% and offsets from the float2integer function.
%
% Run DCAS on the values in the table using Ned’s approach to handle dust 
% using the pixels above some fSCA threshold (e.g. 0.95). 
% Set some values below some fSCA threshold (e.g. 0.10) to 0.00.
% Build a LUT or a ML function to characterize the relationships.
% Now run each image for each year.



% imag data (raw Rs or TOA)

%calulate pixel charateristics

%MAKE THIS A FUNCTION

%For the 3 bands, I suggest: 
%(i) norm of the reflectances; 
%(ii) spectral angle of the reflectances compared to some other spectrum, could be R0 or nominal snow pixel or . . .; 
%(iii) something else, could be NDSI or illumination angle. 

%generate superpixels - either from charatersitcs or raw band values.

%%option for SLIC or SLICSAM

%superpixel charateristcs (Label matrics and region props) -extract from
%my code.

spFTable = spFeatureSets(L,pixelCharateristics);

%median band reflectances of each superpixel +
%illuminationangle,NDSI,canopy cover, R0?

%reduce the number of superpixels

%eliminate some of small size depending on a threshold.

%convert each value to an ‘int8’ using my float2integer function, 
%which would set any NaN or out of range values to intmin(‘int8’). 

%Make a new table, eliminate any rows where any value is intmin(‘int8’), 
%and then run unique on the table to come up with the set of values to use 
%in building the LUT or the ML function.


