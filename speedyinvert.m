function out = speedyinvert(R,R0,solarZ,Ffile,pshade,dust_thresh,dust)
%stripped down inversion for speed
% input: 
%   R - Nx1 band reflectance as vector, center of bandpass
%   R0 - Nx1 band background reflectance
%   solarZ - solar zenith angle for flat surface, deg, scalar
%   Ffile, location of griddedInterpolant with 4 inputs: radius (um), 
% dust (ppm), solarZ (deg), 
%   and band for a specific sensor, e.g. LandSat 8 OLI or MODIS
%   pshade:  shade spectra (bx1)
%   dust_thresh - threshhold value for dust retrievals, e.g. 0.85
% dust -dust val (ppmw), [] if needs to be solved for
% output:
%   out: fsca, fshade, grain radius (um), and dust conc (ppm)
persistent F
if isempty(F)
   X=load(Ffile);
   F=X.F;
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

%solve for fsca/shade again
% knowns=false;
% 
% if ~isempty(dust)
%     knowns=true;
%     radius=solveS.radius;
%     dust=solveS.dust;
% end

out.x=NaN(4,1);

A=[1 1 0 0];
b=1;
   
try
    if ~isempty(dust)
        x0=[0.5 0.1 250 dust]; %fsca, fshade,grain size (um), dust (ppm)
        lb=[0 0 30 dust];
        ub=[1 1 1200 dust];
        X = fmincon(@SnowCloudDiff,x0,A,b,[],[],lb,ub,[],options); 
    else
    %try a clean snow solution
        x0=[0.5 0.1 250 0]; %fsca, fshade,grain size (um), dust (ppm)
        lb=[0 0 30 0];
        ub=[1 1 1200 0];
%         X=lsqnonlin(@SnowCloudDiff,x0,lb,ub,options);
        X = fmincon(@SnowCloudDiff,x0,A,b,[],[],lb,ub,[],options);
        X(4)=NaN; %dust is NaN unless...
        %if fsca is above threshold, re-solve for shade & dust 
        if X(1) >= dust_thresh
            x0=[X(1) 0.1 X(3) 0.1]; %fsca,fshade,grain size (um), dust (ppm)
            lb=[X(1) 0 X(3) 0];
            ub=[X(1) 1 X(3) 1000];
%             X=lsqnonlin(@SnowCloudDiff,x0,lb,ub,options);
            X = fmincon(@SnowCloudDiff,x0,A,b,[],[],lb,ub,[],options);
        end
    end
    out.x=X;
catch

      warning(...
      ['solver crashed, skipping: '...
      ' R=%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f;' ...
      'R0=%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f;' ...
      'solarZ:%2.1f'],...
       R(1),R(2),R(3),R(4),R(5),R(6),R(7),...
       R0(1),R0(2),R0(3),R0(4),R0(5),R0(6),R0(7),...
       solarZ);
%      warning('solver crashed, skipping');
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