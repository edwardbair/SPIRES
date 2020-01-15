function F=build_lt(sensor,bands,F,pshade)
%input: sensor, string, e.g. 'LandsatOLI' or 'MODIS'
%bands: 1xN vector indicating bands needed, e.g. 1:7
% output: gridded interpolant with inputs: 
%grain size (um), dust (conc. by 
% weight, ppm), solar zenith angle (deg), and band (N)

sT=SensorTable(sensor);

radius=30:10:1200;
dust=[0 0.1 1:10:1000];
solarZ=0:1:90;

lTbl=zeros(length(radius),length(dust),...
    length(solarZ),length(bands));
N=length(radius)*length(dust);

n=0;
tic;
for i=1:length(radius)
        for j=1:length(dust)
            parfor k=1:length(solarZ)
                    out=SnowCloudSpectralRefl('snow','cosZ',cosd(solarZ(k)),...
                       'radius',radius(i),...
                        'dust',dust(j)*1e-6,'wavelength',...
                         sT.CentralWavelength(bands),'waveu','um');
                    lTbl(i,j,k,:)=out.refl;
            end
            n=n+1;
               etime=toc;
               fprintf('pct done=%2.5f; time=%2.5f hr\n',...
                   n/N*100,etime/60/60);      
        end
end
[w,x,y,z]=ndgrid(radius,dust,solarZ,1:length(bands));
F=griddedInterpolant(w,x,y,z,lTbl,'pchip','nearest');

%sname='/raid/sandbox/snowhydro/nbair/SMARTS/SMARTS_295_Linux/smarts295.ext.txt';
% S=getSMARTSspectrum(sname);
% SolarT=table(S.waveL*1e-3,[S.HorzDirect,S.HorzDiffuse]*1e3,...
%     'VariableNames',{'wavelength','irradiance'});
% SolarT.Properties.VariableUnits={'um','W m^(-2) um^(-1)'};
% T=SensorTable('MODIS');
% T=T(1:7,:);
% [~,idx]=sort(T.CentralWavelength);
%                     R0T=table(T.CentralWavelength(idx),[thisR0,pshade],...
%                     'VariableNames',{'wavelength','reflectance'});
                   %would take 126,350 hr, use LUT 
%                 outS=invertSnowCloudIntgRefl(SolarT,thisR,unknowns,...
%                     'snow',...
%                     'R0',R0T,'cosZ',cosd(thisSolarZ),'bandPass',...
%                     [T.LowerWavelength(idx) T.UpperWavelength(idx)],...
%                     'waveu','um');