function [S] = SnowCloudLimits()
%limits of snow and cloud properties used in various functions
%
%Output
%   S - structure of 2-element vectors with limits for the following
%       variables
%   snowRadius, mum
%   iceCloudRadius, mum
%   waterCloudRadius, mum
%   iceCloudWE, mm
%   waterCloudWE, mm
%   dust (in snow) mass fraction
%   soot (in snow) mass 
%   dust radius
%   soot radius

S.snowRadius = [30 2000];
S.waterInSnowRadius = [10 100];
S.iceCloudRadius = [5 50];
S.waterCloudRadius = [1 25];
S.dustRadius = [.1 10];
S.sootRadius = [0.01 5];
S.defaultSnowRadius = 250;
S.defaultDustRadius = 3;
S.defaultSootRadius = .5;
S.defaultIceCloudRadius = 10;
S.defaultWaterCloudRadius = 5;
S.iceCloudWE = [.01 3];
S.waterCloudWE = [.01 10];
S.dust = [10*eps 1000]/1e6;
S.soot = [eps 5000]/1e8;
S.wetSnow = [0 .2];
S.unitsSize = 'mum';
S.unitsWE = 'mm';
S.snowSSA = sort(radius2SSA(S.snowRadius,S.unitsSize));
S.iceCloudSSA = sort(radius2SSA(S.iceCloudRadius,S.unitsSize));

end