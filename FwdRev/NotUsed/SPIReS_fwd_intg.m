function [R,varargout ] = SPIReS_fwd_intg(varargin)
%band-integrated reflectance and transmittance of snow
% usage 1: R=SPIReS_fwd_intg(substance,prop1,val1,prop2,val2,...)
% usage 2: R=SPIReS_fwd_intg(substance,prescription,prop1,val1,prop2,val2,...)
% usage 3: R=SPIReS_fwd_intg(prescription,prop1,val1,prop2,val2,...)
% in all cases, you can call as [R,P] = to retrieve the prescription used,
% as modified
%
%Inputs
%Substance and prescription arguments, if specified, must come before the
%name-value pairs if either or both are specified
%   substance, either 'snow', 'iceCloud', 'waterCloud', or 'mixedCloud'
%       (but any unambiguous abbreviation beginning with first letter works)
%   prescription, if you want to modify an existing prescription
%
%Other inputs, name-value pairs, SEE SetSpectral.m FOR ALL POSSIBILITIES
%for example
% 'cosZ' - cosine of illumination angle, scalar, vector, or matrix
% 'elevation' - in meters, needed to calculate spectral distribution of
%   irradiance
% 'radius' or 'SSA' - effective optical radius, or specific surface area,
%   of scatterer, scalar, vector, or matrix
%   (if not scalars, 'cosZ' and 'radius' must be same size as wavelength)
% 'lookup' - logical, use lookup tables for Mie parameters, default false
% 'bandpass' - Nx2 wavelengths for which to calculate values
%   or
% 'sensor' - any multispectral sensor in the SensorTable.m function
%   and
% 'bands' - default all bands for that sensor
% 'LAP' - followed either by 'dust','soot', or a table with wavelength and
%   complex refractive indices.
% 'lfraction' - LAP fraction by weight, dimensionless
%
%Output
% R - table with snow or cloud reflectance and transmittance, dimensionless,
%   same size dimensions as number of bands or multiples of bands x size of
%   other variables
% P - optional, input prescription, returned for reuse

narginchk(1,Inf)
nargoutchk(0,2)

%% call SolarScale once to get the wavelengths needed
iStruct = SetSpectral(mfilename,varargin{:});
STbl = SolarScale('units',iStruct.waveUnit);
%just the wavelengths we need
wavelength = wavelengthsNeeded(STbl.wavelength,iStruct.bandPass);

R = table;
[elevation,cosZ,muS,radius] = checkSizes(iStruct.elevation,iStruct.cosZ,...
    iStruct.muS,iStruct.iceRadius);
for k=1:length(elevation)
    if k==1 || elevation(k)~=elevation(k-1) || cosZ(k)~=cosZ(k-1)
        STbl = SolarScale('units',iStruct.waveUnit,'elevation',elevation(k),...
            'cosZ',cosZ(k),'wavelength',wavelength,'smooth',false);
    end
    STbl.irradiance = [STbl.Global.*(1-STbl.DiffuseFraction)...
        STbl.Global.*STbl.DiffuseFraction];
    iStruct.wavelength=STbl.wavelength;
    [RTbl,P]=SPIReS_fwd(iStruct);
%     [RTbl,P] =  SPIReS_fwd(iStruct,'wavelength',STbl.wavelength,...
%         'cosZ',muS(k),'muS',muS(k),'radius',radius(k));
    if k==1 || elevation(k)~=elevation(k-1) || cosZ(k)~=cosZ(k-1) ||...
            radius(k)~=radius(k-1)
        Rd =  SPIReS_fwd(P,'cosZ',[],'muS',[]);
    end
    
    %% call bandPassReflectance
    reflTbl = table(RTbl.wavelength,[RTbl.reflectance Rd.reflectance],...
        'VariableNames',{'wavelength','reflectance'});
    reflTbl.Properties.VariableUnits = {P.waveUnit,''};
    thisR = bandPassReflectance(reflTbl,STbl,...
        'bandPass',iStruct.bandPass,'bands',iStruct.bands,...
        'muS',muS(k),'mu0',cosZ(k));
    R = [R;table(elevation(k),radius(k),cosZ(k),muS(k),thisR,'VariableNames',...
        {'elevation','radius','cosZ','muS','reflectance'})]; %#ok<AGROW>
end

%% optional output
if nargout>1
    fn = fieldnames(iStruct);
    for k=1:length(fn)
        if isnumeric(iStruct.(fn{k}))
            if isscalar(unique(iStruct.(fn{k})))
                iStruct.(fn{k}) = unique(iStruct.(fn{k}));
            end
        end
    end
    varargout{1} = iStruct;
end

end