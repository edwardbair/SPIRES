function [Tbl,plotTbl] = analyzeResults(inputStruct)
% Tbl = analyzeResults(inputStruct)
%Analyze results from synthesizeInversion

invMeth = inputStruct.solutionMethod;
if contains(invMeth,'lsq','IgnoreCase',true)
    valVar = 'resnorm';
else
    valVar = 'specAngle';
end
inputTbl = inputStruct.inputTbl;
oTbl = inputStruct.oTbl;
aTbl = inputStruct.albedoTbl;

% correct oTbl heading
vn = oTbl.Properties.VariableNames;
vn{6} = 'fSCA_soot';
oTbl.Properties.VariableNames = vn;

% variables to analyze
Variables = {'ssa','fSCA_dust','fSCA_soot','cleanA','dustyA','sootyA'};

for k=1:length(Variables)
    switch k
        case 1
            [M,G] = fit(inputTbl.ssa,oTbl.ssa,'poly1','robust','bisquare');
            Ci = confint(M);
            if Ci(1,1)<=1 && Ci(2,1)>=1
                bias = 0;
            else
                bias = M.p1-1;
            end
            plotTbl{k} = table(inputTbl.ssa,oTbl.ssa,'VariableNames',...
                {'inputSSA','solvedSSA'}); %#ok<*AGROW>
            rmse = G.rmse;
            rsquare = G.adjrsquare;
            if contains(invMeth,'lsq','IgnoreCase',true)
                F = fitdist(oTbl.resnorm,'gamma');
            else
                F = fitdist(oTbl.spectralAngle,'gamma');
            end
            gammaA = F.a;
            gammaB = F.b;
            rn = {['SSA ' invMeth]};
            thisTbl = table(rmse,rsquare,bias,gammaA,gammaB,'RowNames',rn);
        case {2,3}
            [M,G] = fit(inputTbl.fSCA,oTbl.(Variables{k}),'poly1','robust','bisquare');
            Ci = confint(M);
            if Ci(1,1)<=1 && Ci(2,1)>=1
                bias = 0;
            else
                bias = M.p1-1;
            end
            plotTbl{k} = table(inputTbl.fSCA,oTbl.(Variables{k}));
            switch k
                case 2
                    vn = {'input_fSCA','solved_fSCA_dust'};
                case 3
                    vn = {'input_fSCA','solved_fSCA_soot'};
            end
            plotTbl{k}.Properties.VariableNames = vn;
            rmse = G.rmse;
            rsquare = G.adjrsquare;
            switch k
                case 2
                    tvar = ['dust_' valVar];
                case 3
                    tvar = ['soot_' valVar];
            end
            F = fitdist(oTbl.(tvar),'gamma');
            gammaA = F.a;
            gammaB = F.b;
            rn = {[Variables{k} ' ' invMeth]};
            thisTbl = table(rmse,rsquare,bias,gammaA,gammaB,'RowNames',rn);
        case {4,5,6}
            xvar = Variables{k};
            yvar = [Variables{k} '_solve'];
            x = aTbl.(xvar);
            y = aTbl.(yvar);
            t = inputTbl.fSCA>=0.3;
            [M,G] = fit(x,y,'poly1','robust','bisquare');
            [M2,G2] = fit(x(t),y(t),'poly1','robust','bisquare');
            bias = zeros(2,1);
            Ci = confint(M);
            if Ci(1,1)>1 || Ci(2,1)<1
                bias(1) = M.p1-1;
            end
            Ci = confint(M2);
            if Ci(1,1)>1 || Ci(2,1)<1
                bias(2) = M2.p1-1;
            end
            plotTbl{k} = table(inputTbl.fSCA,x,y,'VariableNames',{'fSCA',xvar,yvar});
            rn = {[xvar ' ' invMeth ' all'],[xvar ' ' invMeth ' >thresh']};
            switch k
                case 4
                    tvar = valVar;
                    if strcmp(tvar,'specAngle')
                        tvar = 'spectralAngle';
                    end
                case 5
                    tvar = ['dust_' valVar];
                case 6
                    tvar = ['soot_' valVar];
            end
            v = oTbl.(tvar);
            F = fitdist(v,'gamma');
            gA(1) = F.a;
            gB(1) = F.b;
            F = fitdist(v(t),'gamma');
            gA(2,1) = F.a;
            gB(2,1) = F.b;
            thisTbl = table([G.rmse G2.rmse]',[G.adjrsquare G2.adjrsquare]',...
                bias,gA,gB,'VariableNames',{'rmse','rsquare','bias','gammaA','gammaB'},...
                'RowNames',rn);
    end
    if k==1
        Tbl = thisTbl;
    else
        Tbl = [Tbl; thisTbl];
    end
    
    
end

