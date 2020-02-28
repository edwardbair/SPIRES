function [ rescaledValues ] = convertLandsat8DN( Qcal, MTL_info, band)
%[ Rp ] = ConvertLandsat8DN( Qcal, MTL, band) - bands 1-9
%[ BT ] = ConvertLandsat8DN( Qcal, MTL, band) - bands 10,11
%
%DNtoRp Convert 16bit landsat 8 digital number values to top
%of atmosphere planetarty reflectance corrected for sun angle for bands 1-9
%and convert DNtoBT for thermal bands 10 & 11.
%
%Inputs
%   Qcal - Quantized and calibrated Landsat 8 pixel values (16 bit digital number)
%   band - Landsat 8 band the Qcal values are from
%   MTL - Landsat 8 scene metadata file
%
%Formulas from:
%   http://landsat.usgs.gov/Landsat8_Using_Product.php
%
%
% 8/16/2014 - add ability to compute thermal bands BT

b = num2str(band);

if ischar(band)
    band = str2double(band);
end

%convert Qcal to double
if ~isfloat(Qcal)
    Qcal = cast(Qcal,'double');
end

%change 0 values in Qcal to NaN
Qcal(Qcal == 0) = NaN;

if nnz(band == 1:9) == 1
    %multiplicative rescaling factor
    multi_param = 'REFLECTANCE_MULT_BAND_';
    MpB = strcat(multi_param, b);
    M = MTL_info.RADIOMETRIC_RESCALING.(MpB);
    
    %additive rescaling factor
    add_param = 'REFLECTANCE_ADD_BAND_';
    ApB = strcat(add_param, b);
    A = MTL_info.RADIOMETRIC_RESCALING.(ApB);
    
    %sun elevation
    sunElv =  MTL_info.IMAGE_ATTRIBUTES.SUN_ELEVATION;
    
    % Top of atmosphere reflectance w/ correction for sun angle
    Rprime = M .* Qcal + A;
    Rp = Rprime ./ sind(sunElv);
    
    rescaledValues = Rp;
    
elseif nnz(band == 10:11) == 1
    %multiplicative rescaling factor
    multi_param = 'RADIANCE_MULT_BAND_';
    MpB = strcat(multi_param, b);
    M = MTL_info.RADIOMETRIC_RESCALING.(MpB);
    
    %additive rescaling factor
    add_param = 'RADIANCE_ADD_BAND_';
    ApB = strcat(add_param, b);
    A = MTL_info.RADIOMETRIC_RESCALING.(ApB);
    
    %band specific thermal conversion constants
    K1_param = 'K1_CONSTANT_BAND_';
    K1B = strcat(K1_param, b);
    K1 = MTL_info.TIRS_THERMAL_CONSTANTS.(K1B);
    
    K2_param = 'K2_CONSTANT_BAND_';
    K2B = strcat(K2_param, b);
    K2 = MTL_info.TIRS_THERMAL_CONSTANTS.(K2B);
    
    % Top of Atmosphere Radiance (Watts/(m2*srad*um)
    L = M .* Qcal + A;
    
    %At-Satellite Brightness Temperature
    BT = K2 ./ log((K1./L) + 1);
    
    rescaledValues = BT;
    
else
    warning('convertLandsat8DN.m only works with LS8 bands 1-11.')
end
end
