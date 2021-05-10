function X = albedoFits(cosZ,radius,refl,threshZ)
% X = albedoFits(cosZ,radius,refl)
%fit values for power function, f(radius) separately for each cosZ
%then fit power function coefficients a f(cosZ)

%fit power function to each cosZ as f(radius)
a = zeros(size(cosZ));
b = zeros(size(a));
c = zeros(size(a));
wt = ones(size(a));
wt(cosZ<cosd(threshZ)) = 0;
for k=1:length(cosZ)
    F = fit(radius,refl(:,k),'power2');
    a(k) = F.a;
    b(k) = F.b;
    c(k) = F.c;
end
Fsize.a = a;
Fsize.b = b;
Fsize.c = c;
Fsize.exampleF = F;

%fit a,b,c coefficients to cosZ and select best fits (R2>=0.999)
FcosZ = cell(3,1);
GcosZ = cell(size(FcosZ));
cval = cell(1,3);
rsq = cell(size(cval));
rmse = cell(size(cval));
startVals = cell(size(cval));
P = zeros(3,3);
Q = zeros(size(P));
rowP = zeros(1,3);
rowQ = zeros(1,3);
for k=1:3
    switch k
        case 1
            var = a;
            f = 'rat22';
        case 2
            var = b;
            f = 'rat22';
        case 3
            var = c;
            f = 'rat21';
    end
    rng('shuffle')
    V = optFit(var,f,cosZ,wt);
    FcosZ{k} = V.Fout;
    GcosZ{k} = V.Gout;
    cval{k} = V.cval;
    rsq{k} = V.rsq;
    rmse{k} = V.rmse;
    startVals{k} = V.startVals;
    cname = coeffnames(V.Fout);
    ncoef = length(cname);
    for n=1:3
        rowP(n) = V.Fout.(cname{n});
    end
    if ncoef==5
        rowQ = [1 V.Fout.(cname{end-1}) V.Fout.(cname{end})];
    elseif ncoef==4
        rowQ = [0 1 V.Fout.(cname{end})];
    end
    P(k,:) = rowP;
    Q(k,:) = rowQ;
end
X.Fsize = Fsize;
X.FcosZ = FcosZ;
X.GcosZ = GcosZ;
X.cval = V.cval;
X.rsq = rsq;
X.rmse = rmse;
X.startVals = startVals;
X.P = P;
X.Q = Q;

end

function X = optFit(var,f,cosZ,wt)

numFit = 100;
thresh = 0.999;
fun = fittype(f);
cname = coeffnames(fun);
ncoef = length(cname);
startVals = rand(numFit,ncoef);
rsq = zeros(numFit,1);
rmse = zeros(size(rsq));
cval = zeros(ncoef,numFit);
Fhold = cell(numFit,1);
Ghold = cell(size(Fhold));

%eliminate values for which Weights are zero
t = wt==0;
var = var(~t);
cosZ = cosZ(~t);
wt = wt(~t);
if any(wt~=1)
    fitOpt = fitoptions(f,'Robust','bisquare','Weights',wt);
else
    fitOpt = fitoptions(f,'Robust','bisquare');
end

% try numfit estimates of coefficients
for k=1:numFit
    fitOpt.StartPoint = startVals(k,:);
    [F,G,~] = fit(cosZ,var,f,fitOpt);
    for n=1:ncoef
        cval(n,k) = F.(cname{n});
    end
    rsq(k) = G.adjrsquare;
    rmse(k) = G.rmse;
    Fhold{k} = F;
    Ghold{k} = G;
end
t = rsq>=thresh;
psq = prctile(rsq,95);
pse = prctile(rmse,5);
if nnz(t)==0
    warning('max value of rsq (%f) < thresh',max(rsq))
    tx = rsq>=psq & rmse<=pse;
    assert(nnz(tx)>0,...
        'no values of rsq above 95th percentile and rmse below 5th percentile')
else
    tx = t & rsq>=psq & rmse<=pse;
    assert(nnz(tx)>0,...
        'no values of rsq>=95th percentile and rmse<=5th percentile')
end
% reduce arrays to values that meet the criteria
rsq = rsq(tx);
rmse = rmse(tx);
cval = cval(:,tx);
startVals = startVals(tx,:);
k1 = find(rsq==max(rsq));
k2 = find(rmse==min(rmse));
Fhold = Fhold(tx,:);
Ghold = Ghold(tx,:);
if k1==k2
    Fout = Fhold{k1};
    Gout = Ghold{k1};
else
    [~,Irsq] = sort(rsq,'descend');
    [~,Irmse] = sort(rmse);
    [~,iSelect] = sort(sum([Irsq Irmse],2));
    Fout = Fhold{iSelect(1)};
    Gout = Ghold{iSelect(1)};
end

X.Fout = Fout;
X.Gout = Gout;
X.cval = cval;
X.rsq = rsq;
X.rmse = rmse;
X.startVals = startVals;
end