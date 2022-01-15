function [Tbl,TOAtbl,S] = prepSolarScale(refAtmos)
% [Tbl,TOAtbl,S] = prepSolarScale(refAtmos)
%prep data to build interpolants for SolarScale
%run SMARTS with combinations of solar zenith ('cosZ') and elevation ('altit')
%

%TOA will be the same for all
P = defaultSMARTSinput(refAtmos,'cosZ',cosd(48.19),'altit',3000);
S = SMARTSMain(P);
T = S.spectralTbl;
F = fit(T.waveL,T.TOA,'pchip');
I = integrate(F,T.waveL(end),T.waveL(1));
y = T.TOA/I;
TOAtbl = table(T.waveL,y,'VariableNames',{'wavelength','TOA'});

cosZ = linspace(cosd(85),1,20);
elev = 1000*(1:.5:7);
[c,e] = ndgrid(cosZ,elev);
c = c(:);
e = e(:);
[wv,cz,alt] = ndgrid(T.waveL,cosZ,elev);
origSize = size(wv);
wv = wv(:);
cz = cz(:);
alt = alt(:);
HorzGlobal = zeros(size(alt));
diffFrac = zeros(size(alt));

for k=1:length(c)
    P = defaultSMARTSinput(refAtmos,'cosZ',c(k),'altit',e(k));
    S = SMARTSMain(P);
    T = S.spectralTbl;
    
    hG = max(T.HorzGlobal,T.HorzDirect+T.HorzDiffuse);
    dF = T.HorzDiffuse./hG;
    t = hG>0 & ~isnan(dF);
    FG = fit(T.waveL(t),hG(t),'pchip');
    I = integrate(FG,T.waveL(end),T.waveL(1));
    hG = FG(T.waveL)/I;
    FD = fit(T.waveL(t),dF(t),'pchip');
    dF = FD(T.waveL);
    dF(dF>1) = 1;
    dF(dF<0) = 0;
    
    % fill these values into the output vectors
    t = cz==c(k) & alt==e(k);
    HorzGlobal(t) = hG;
    diffFrac(t) = dF;
    if mod(k,50)==0
        disp([k c(k) e(k)])
    end
end
Tbl = table(wv,cz,alt,HorzGlobal,diffFrac,...
    HorzGlobal.*(1-diffFrac),HorzGlobal.*diffFrac,...
    'VariableNames',...
    {'wavelength','cosZ','elev','Global','DiffuseFraction','Direct','Diffuse'});
Tbl.Properties.UserData = origSize;
fn = Tbl.Properties.VariableNames;
S = struct;
for k=1:5
    S.(fn{k}) = reshape(Tbl.(fn{k}),Tbl.Properties.UserData);
end
end