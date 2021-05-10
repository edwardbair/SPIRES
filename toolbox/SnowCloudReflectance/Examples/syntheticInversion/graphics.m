function graphics(pL,pS,logs)
% figh = graphics(Tl,Ts)
%   graphics from analyzeResults
%   inputs are one pair of lsqnonlin and spectralAngle tables

vL = pL.Properties.VariableNames;
vS = pS.Properties.VariableNames;
assert(isequal(vL,vS),'variable names differ, check input')
plot(pL.(vL{1}),pL.(vL{2}),'LineStyle','none','Marker','.','MarkerSize',10),
axis equal tight
hold on
plot(pS.(vS{1}),pS.(vS{2}),'LineStyle','none','Marker','+','MarkerSize',8);
axis equal tight
if logs
    set(gca,'YScale','log','XScale','log')
end
xlabel(vL{1})
ylabel(vL{2})
legend({'nonlinear lsq','spectral angle'},'Location','Best')
end

