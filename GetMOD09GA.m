function X = GetMOD09GA( file, whichVariable )
% X = GetMOD09GA( file, whichVariable )
%extract variable from MOD09GA input file (or MYD09GA)
%
% Input
%   file - name of MOD09GA or MYD09GA .hdf file
%   whichVariable - choose among (case insensitive)
%       1 km spacing
%           'num_obs_1km', 'state', 'SensorZenith', 'SolarZenith','SolarAzimuth'
%           (others could be added but these are ones we use)
%       500 m spacing
%           'num_obs_500m', 'band1', 'band2', 'band3', ... 'band7',
%           'QC'
%           or 'allbands' to get bands 1 thru 7
%
% Output
%   variable scaled to appropriate values and units (the angle and band
%   variables are returned as floating point single precision, state and QC
%   are returned as structures, others as native size integers)

if iscell(file)
    file = char(file);
end
angle = false;
refl = false;
switch lower(whichVariable)
    case 'num_obs_1km'
        var = 'num_observations_1km';
    case 'state'
        var = 'state_1km_1';
    case lower('SensorZenith')
        var = 'SensorZenith_1';
        angle = true;
    case lower('SolarZenith')
        var = 'SolarZenith_1';
        angle = true;
    case lower('SolarAzimuth')
        var = 'SolarAzimuth_1';
        angle = true;
    case 'num_obs_500m'
        var = 'num_observations_500m';
    case lower('QC')
        var = 'QC_500m_1';
    case lower('allbands')
        X = GetMOD09GA(file,'band1');
        for b=2:7
            X = cat(3,X,GetMOD09GA(file,['band' num2str(b)]));
        end
        return
    otherwise
        if strncmpi(whichVariable,'band',4)
            bandNo = str2double(whichVariable(5:end));
            var = ['sur_refl_b' num2str(bandNo,'%02d') '_1'];
            refl = true;
        else
            error('Variable %s unrecognized',whichVariable)
        end
end

angleScaleDivisor = 100;
reflScaleDivisor = 10000;

iX = hdfread(file,var);
if angle
    X = single(iX)/angleScaleDivisor;
    X(X<-180) = NaN;
elseif refl
    X = single(iX)/reflScaleDivisor;
    X(X<0) = NaN;
elseif strncmpi(whichVariable,'state',5)
    X = unpackMOD09state(iX);
elseif strncmpi(whichVariable,'qc',2);
    X = unpackMOD09QC(iX);
else
    X = iX;
end

end