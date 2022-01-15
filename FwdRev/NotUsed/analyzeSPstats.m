function chooseIdx = analyzeSPstats(file,specS,Ftran,Fdf)
%analyze output from fitAllSuperpixels

load(file); %#ok<LOAD>
assert(exist('Rlsq','var') && exist('Rn','var') && exist('Rs','var'),...
    'check input file ''%s''',file)

slim = SnowCloudLimits;

% add a zero shade endmember to Rlsq to make tables same size
Rlsq.fSCA = cat(2,Rlsq.fSCA,zeros(size(Rlsq.fSCA,1),1)); %#ok<NODEF>

% list as tables too
% Tlsq = struct2table(Rlsq);
% Tn = struct2table(Rn);
% Ts = struct2table(Rs);

% good values for each
glsq = Rlsq.normResiduals<=2 & Rlsq.idx>0 & Rlsq.radius>50 & Rlsq.radius<1200 & Rlsq.fSCA(:,1)>=.1;
gn = Rn.normResiduals<=2 & Rn.idx>0 & Rn.radius>50 & Rn.radius<1200 & Rn.fSCA(:,1)>=.1 & round(Rn.radius)~=slim.defaultSnowRadius;
gs = Rs.idx>0 & Rs.radius>50 & Rs.radius<1200 & Rs.fSCA(:,1)>=.1;
% allg = glsq & gn & gs; % good for all
anygood = glsq | gn | gs;
allg = glsq & gn;
Gcode = zeros(length(allg));
Gcode(glsq) = 1;
Gcode(gn) = 2;

% graph the comparisons
subplot(2,2,1)
gscatter(Rlsq.radius(allg),Rn.radius(allg),Gcode(allg))
title('''lsq'' vs ''norm''')
xlabel('radius')
axis equal tight
subplot(2,2,2)
scatter(Rlsq.radius(allg),Rs.radius(allg),'*');
title('''lsq'' vs ''spec''')
xlabel('radius')
axis equal tight
subplot(2,2,3)
gscatter(Rlsq.fSCA(allg,1),Rn.fSCA(allg,1),Gcode(allg))
title('''lsq'' vs ''norm''')
xlabel('fSCA')
axis equal tight
subplot(2,2,4)
scatter(Rlsq.fSCA(allg,1),Rs.fSCA(allg,1),'.')
title('''lsq'' vs ''spec''')
xlabel('fSCA')
axis equal tight

% compare spectralP, slope of fit
% idxG = Rlsq.idx(allg);
% t = allg & abs(sqrt(Rlsq.radius)-sqrt(Rn.radius))>5;
% tlsq = t & abs(sqrt(Rlsq.radius)-sqrt(Rs.radius))<abs(sqrt(Rn.radius)-sqrt(Rs.radius));
% tn = t & abs(sqrt(Rlsq.radius)-sqrt(Rs.radius))>=abs(sqrt(Rn.radius)-sqrt(Rs.radius));
fR0 = fit(specS.wavelength(:),specS.backReflectance(:),'pchip');
fR2 = fit([min(specS.wavelength) max(specS.wavelength)]',[0 0]','poly1');

% choose when one of the methods gives results
% column 1 is idx #, col2 shows ok slopes where IDs are
% 1 lsq, 2 norm, 3 spec but omitted if slope does not span 1
% col3 shows differences to 1.0, col 4 shows goodness
chooseIdx = struct;
chooseIdx.idx = zeros(length(allg),1);
chooseIdx.qflag = zeros(length(allg),1,'int8')-1; % set qual flag to -1
chooseIdx.whichMethod = false(length(allg),3);
chooseIdx.slope = nan(length(allg),3);
chooseIdx.good = nan(length(allg),3);

for k=find(anygood)'
    if glsq(k)
        Plsq = fwdPrescription('snow','cosZ',specS.topography(k,2),'elevation',...
            specS.topography(k,1),'muS',Rlsq.muS(k),'fSCA',Rlsq.fSCA(k,1:2),'R0',fR0,...
            'waveu','nm','wavelength',specS.wavelength,'radius',Rlsq.radius(k),...
            'diffuseFraction',Fdf,'LAPfraction',Rlsq.LAPfraction(k),'LAP','dust');
        lsqR =  SPIReS_fwd(Plsq);
    else
        lsqR = nan(size(specS.wavelength));
    end
    if gn(k)
        Pn = fwdPrescription('snow','cosZ',specS.topography(k,2),'elevation',...
            specS.topography(k,1),'muS',specS.topography(k,3),'fSCA',Rn.fSCA(k,:),...
            'R0',{fR0,fR2},'LAPfraction',Rn.LAPfraction(k),'LAP','dust',...
            'waveu','nm','wavelength',specS.wavelength,'radius',Rn.radius(k),...
            'diffuseFraction',Fdf);
        nR =  SPIReS_fwd(Pn);
    else
        nR = nan(size(specS.wavelength));
    end
    if gs(k)
        Ps = fwdPrescription('snow','cosZ',specS.topography(k,2),'elevation',...
            specS.topography(k,1),'muS',specS.topography(k,3),'fSCA',Rs.fSCA(k,:),...
            'R0',{fR0,fR2},'LAPfraction',Rs.LAPfraction(k),'LAP','dust',...
            'waveu','nm','wavelength',specS.wavelength,'radius',Rs.radius(k),...
            'diffuseFraction',Fdf);
        nS =  SPIReS_fwd(Ps);
    else
        nS = nan(size(specS.wavelength));
    end
    wt = getAtmosProp(Ftran,specS.wavelength,specS.topography(k,1),specS.topography(k,2));
    if mod(k,400)==0
        figure
        plot(specS.wavelength,[wt.*lsqR wt.*nR wt.*nS],'linewidth',1)
        hold on;
        plot(specS.wavelength,specS.spectralWeightedSP(k,:),'k','linewidth',1)
        hold off
        title(num2str(k))
    end
    
    % compare just the weighted linear fits
    % turn off warning about variable name change
    warnID = 'curvefit:fit:iterationLimitReached';
    warnStruct = warning('off',warnID);
    ok = false(1,3);
    slope = nan(size(ok));
    good = nan(size(ok));
    % which slopes have confidence intervals that span 1.0
    for n=1:3
        switch n
            case 1
                t = glsq(k);
                x = lsqR(:);
            case 2
                t = gn(k);
                x = nR(:);
            case 3
                t = gs(k);
                x = nS(:);
        end
        if t
            [F,G,~] = fit(x,specS.spectralSP(k,:)','poly1','robust','bisquare','weights',wt);
            ci = confint(F);
            if (ci(1,1)-1)*(ci(2,1)-1)<0
                ok(n) = true;
            end
            slope(n) = F.p1;
            good(n) = G.adjrsquare;
        end
    end
    if any(ok)
        chooseIdx.idx(k) = k;
        chooseIdx.qflag(k) = nnz(ok);
        chooseIdx.whichMethod(k,:) = ok;
        chooseIdx.slope(k,:) = slope;
        chooseIdx.good(k,:) = good;
    else
        chooseIdx.idx(k) = k;
        chooseIdx.qflag(k) = 0;
        bestworstchoice = find(abs(slope-1)==min(abs(slope-1)));
        chooseIdx.whichMethod(k,bestworstchoice) = true;
        chooseIdx.slope(k,bestworstchoice) = slope(bestworstchoice);
        chooseIdx.good(k,bestworstchoice) = good(bestworstchoice);
    end
end
%turn warning back on
warning(warnStruct);
end