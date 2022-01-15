function S = assembleSPdata(chooseIdx)
% S = assembleSPdata(chooseIdx)
% use chooseIdx output from analyzeSPstats.m, partition data depending on
% qflag from chooseIdx in high confidence, medium, low and likely not snow
%
% return values in an explanatory structure

load('SPstats.mat','Rlsq','Rn','Rs')
load('newS.mat','newS')
load('atmosWeightAVIRIS-NG.mat','Ftrans','Fdf')

% quality flags
t2 = chooseIdx.qflag==2; % high confidence
tl = chooseIdx.qflag==1 & chooseIdx.whichMethod(:,1); % medium
tn = chooseIdx.qflag==1 & chooseIdx.whichMethod(:,2); % medium
t0 = chooseIdx.qflag==0; % low
tnot = chooseIdx.qflag==-1; % not snow?
% save output t values
S.thigh = t2;
S.tmedium = tl|tn;
S.tmaybe = t0;
S.tnot = tnot;

% x-axis
xl = [300 2580];
% background
fR0 = fit(newS.wavelength(:),newS.backReflectance(:),'pchip');
fR2 = fit([min(newS.wavelength) max(newS.wavelength)]',[0 0]','poly1');

% high, save mean of the values (use sqrt for radius, log for LAPfraction)
S.highConf.idx = Rlsq.idx(t2);
S.highConf.radius = mean(sqrt([Rlsq.radius(t2) Rn.radius(t2)]),2).^2;
S.highConf.fSCA = mean([Rlsq.fSCA(t2,1) Rn.fSCA(t2,1)],2);
S.highConf.LAPfraction = exp(mean([log(Rlsq.LAPfraction(t2)) log(Rn.LAPfraction(t2))],2));

% medium, save the one that met the fit criterion
S.mediumConf.idx = [Rlsq.idx(tl);Rn.idx(tn)];
S.mediumConf.radius = [Rlsq.radius(tl);Rn.radius(tn)];
S.mediumConf.fSCA = [Rlsq.fSCA(tl,1);Rn.fSCA(tn,1)];
S.mediumConf.LAPfraction = [Rlsq.LAPfraction(tl);Rn.LAPfraction(tn)];

% low, save the only one that's available
txl = t0 & ~isnan(chooseIdx.slope(:,1));
txn = t0 & ~isnan(chooseIdx.slope(:,2));
txs = t0 & ~isnan(chooseIdx.slope(:,3)); % spectral angle
S.lowConf.idx = [Rlsq.idx(txl);Rn.idx(txn);Rs.idx(txs)];
S.lowConf.radius = [Rlsq.radius(txl);Rn.radius(txn);Rs.radius(txs)];
S.lowConf.fSCA = [Rlsq.fSCA(txl,1);Rn.fSCA(txn,1);Rs.fSCA(txs,1)];
S.lowConf.LAPfraction = [Rlsq.LAPfraction(txl);Rn.LAPfraction(txn);Rs.LAPfraction(txs)];

% not snow
S.notSnow.idx = Rlsq.idx(tnot);
S.notSnow.idx(S.notSnow.idx==0) = [];

% plot the spectra
figure
subplot(2,2,1)
plot(newS.wavelength,newS.spectralSP(t2,:)')
xlabel('wavelength, nm')
title('high confidence')
xlim(xl)
subplot(2,2,2)
plot(newS.wavelength,newS.spectralSP(tl|tn,:)')
xlabel('wavelength, nm')
title('medium confidence')
xlim(xl)
subplot(2,2,3)
plot(newS.wavelength,newS.spectralSP(txl|txn|txs,:)')
xlabel('wavelength, nm')
title('low confidence')
xlim(xl)
subplot(2,2,4)
plot(newS.wavelength,newS.spectralSP(tnot,:)')
xlabel('wavelength, nm')
title('not snow?')
xlim(xl)

% examples of high, medium, low confidence
figure
subplot(1,3,1)
thigh = t2 & Rlsq.muS<.9;
thisIdx = find(thigh);
n = randperm(length(thisIdx),1);
showConf(thisIdx(n));
subplot(1,3,2)
tmed = (tl|tn) & Rlsq.muS<.9;
thisIdx = find(tmed);
n = randperm(length(thisIdx),1);
showConf(thisIdx(n));
subplot(1,3,3)
tlow = (txl|txn) & Rlsq.muS<.9;
thisIdx = find(tlow);
n = randperm(length(thisIdx),1);
showConf(thisIdx(n));

    function showConf(k)
        % return reflectances
        Plsq = fwdPrescription('snow','cosZ',newS.topography(k,2),'elevation',...
            newS.topography(k,1),'muS',Rlsq.muS(k),'fSCA',Rlsq.fSCA(k,1:2),'R0',fR0,...
            'waveu','nm','wavelength',newS.wavelength,'radius',Rlsq.radius(k),...
            'diffuseFraction',Fdf,'LAPfraction',Rlsq.LAPfraction(k),'LAP','dust');
        lsqR =  SPIReS_fwd(Plsq);
        Pn = fwdPrescription('snow','cosZ',newS.topography(k,2),'elevation',...
            newS.topography(k,1),'muS',newS.topography(k,3),'fSCA',Rn.fSCA(k,:),...
            'R0',{fR0,fR2},'LAPfraction',Rn.LAPfraction(k),'LAP','dust',...
            'waveu','nm','wavelength',newS.wavelength,'radius',Rn.radius(k),...
            'diffuseFraction',Fdf);
        nR =  SPIReS_fwd(Pn);
        % atmos weights
        wt = getAtmosProp(Ftrans,newS.wavelength,newS.topography(k,1),newS.topography(k,2));
        meas = newS.spectralSP(k,:).';
        % fit the results against measurements
        warnID = 'curvefit:fit:iterationLimitReached';
        warnStruct = warning('off',warnID);
        [F1,G1,~] = fit(lsqR,meas,'poly1','robust','bisquare','weights',wt); %#ok<ASGLU>
        [F2,G2,~] = fit(nR,meas,'poly1','robust','bisquare','weights',wt); %#ok<ASGLU>
        scatter(lsqR,meas,'b','*')
        hold on;
        scatter(nR,meas','r','+')
        plot([min(lsqR) max(lsqR)]',F1([min(lsqR) max(lsqR)]'),'b','linewidth',1)
        plot([min(nR) max(nR)]',F2([min(nR) max(nR)]'),'r-.','linewidth',1)
        xlabel('model reflectance')
        ylabel('measured reflectance')
        %turn warning back on
        warning(warnStruct);
        title(['superpixel ' num2str(k)])
        axis equal tight
    end
end

