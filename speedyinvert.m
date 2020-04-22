function [out,modelRefl] = speedyinvert(R,R0,solarZ,Ffile,pshade,...
    dustmask,dust,cc)
%stripped down inversion for speed
% input:
%   R - Nx1 band reflectance as vector, center of bandpass
%   R0 - Nx1 band background reflectance
%   solarZ - solar zenith angle for flat surface, deg, scalar
%   Ffile, location of griddedInterpolant with 4 inputs: radius (um),
% dust (ppm), solarZ (deg),
%   and band for a specific sensor, e.g. LandSat 8 OLI or MODIS
%   pshade:  shade spectra (bx1)
%   dustmask - only retrieve dust values where this is true
% dust -dust val (ppmw), [] if needs to be solved for
% cc - cc val
% output:
%   out: fsca, fshade, grain radius (um), and dust conc (ppm)
persistent F
if isempty(F)
    X=load(Ffile);
    F=X.F;
end

ccflag=true;
if cc==0
    ccflag=false;
end

options = optimoptions('fmincon','Display','none','Algorithm','sqp');

% make all inputs column vectors
if ~iscolumn(R)
    R=R';
end
if ~iscolumn(R0)
    R0=R0';
end
if ~iscolumn(pshade)
    pshade=pshade';
end

out.x=NaN(4,1);

A=[1 1 0 0];
b=1;

%full dirty snow (dustmask) values
%fsca, fshade,grain size (um), dust (ppm)
fsca0=0.5;
fsca_range=[0 1];
fshade0=0.1;
fshade_range=[0 1];
r0=250;
r_range=[30 1200];
d0=10;
d_range=[0 1000];

try
    if ccflag %if there's canopy cover, set shade guess to cc
%         fshade0=cc;
%         fshade_range=[cc 1];
    end    
    if ~isempty(dust) %if dust values are provided, set them
        d0=dust;
        d_range=[dust dust];
    elseif ~dustmask %clean snow solution
        d0=0;
        d_range=[0 0];
    end
    
    x0=[fsca0 fshade0 r0 d0];
    lb=[fsca_range(1) fshade_range(1) r_range(1) d_range(1)];
    ub=[fsca_range(2) fshade_range(2) r_range(2) d_range(2)];
    X = fmincon(@SnowCloudDiff,x0,A,b,[],[],lb,ub,[],options);
    
    out.x=X;
catch ME
    
    warning([ME.message,' solver crashed, skipping']);
end

    function diffR = SnowCloudDiff(x)
        modelRefl=zeros(length(R),1);
        %x is fsca,radius,dust
        for i=1:length(R)
            %use radius,dust,solarZ, and band # for look up
            modelRefl(i)=F([x(3),x(4),solarZ,i]);
        end
        
        modelRefl=x(1).*modelRefl + x(2).*pshade + (1-x(1)-x(2)).*R0;
        diffR = norm(R - modelRefl);
    end
end