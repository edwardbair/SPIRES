% SUNPOSITION
% Functions
%   EarthEphemeris - Returns the solar declination, Earth-Sun radius vector
%       (distance in astronomical units), solar longitude, and optionally the
%       equation of time.
%   sunang - Returns cosine of sun angle and its azimuth in degrees, either
%       clockwise from north or counter-clockwise from south, and airmass
%       (path length through atmosphere accounting for Earth curvature).
%       Optionally accounts for refraction.
%   sunRiseSet - Times of sunrise and sunset.
%   sunslope - Cosine of illumination angle on slope
%
%About Orientation for Azimuths for sun angles and slopes
%   Years ago when I started coding for solar geometry in mountains, my go-to
%   text was William Sellers’ Physical Climatology (1965). Based on his conventions,
%   I adopted 0° toward South, positive to the East and negative to the West,
%   putting North at ±180°. This convention is consistent with a right-hand
%   coordinate system, which we conventionally use by making latitude positive
%   to the North and longitude positive to the East.
%   I realize that another common convention puts 0° toward North, with other
%   azimuths clockwise from north, so I have added this option to sunang.
%   Set your preference in the function azimuthPreference. Because this preference
%   is likely persistent, you can set it in the function and not worry about
%   it again.
%   The sunslope.m code will work without modification as long as the slope
%   azimuth uses the same convention as the solar azimuth.
