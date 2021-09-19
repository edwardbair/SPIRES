function [lap,fval]=EffectiveLAP(deltavis,rg,mu0,Z,varargin)
%solve for effective lap concentration of dust (ppmw) given deltavis
%inputs
%deltavis, expressed as a fraction, eg 0.05
%rg - grain radius, um
%mu0 - solar zenith 0-1
%Z - altitude in km of snowcover, e.g. 3.5
%optional
%pollutant is assumed to be dust but if
%'soot' is provided as a 6th argument, that will be used instead

lap='dust';
x1=0;
x2=1000; %ppmw
dvp=0.63; %Portion of the spectrum that deltavis encompasses.
%Note that deltavis is not properly named as it includes nIR region
% 0.350 to 0.876 um
%Bair et al 2019 WRR

if nargin==5
    if strcmp(varargin{1},'soot')
        lap='soot';
        x2=50; %ppmw
    end
end

[lap,fval]=fminbnd(@ADiff,x1,x2);

%find the rms difference of clean-dirty snow, varying only dust conc
%return the "lap", i.e. effective dust concentration
    function y=ADiff(x)
        dv=AlbedoLookup(rg,mu0,Z)-AlbedoLookup(rg,mu0,Z,lap,x*1e-6);
        y=sqrt((dv-dvp.*deltavis).^2);
    end

end
