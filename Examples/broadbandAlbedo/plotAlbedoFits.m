function [estRefl] = plotAlbedoFits(cosZ,radius,P,Q)
% [estRefl] = plotAlbedoFits(cosZ,radius,P,Q)
% plot estimated reflectance from fits
estRefl = zeros(length(radius),length(cosZ));
for k=1:length(cosZ)
    if cosZ(k)<0.07
        mu0 = 0.07;
    else
        mu0 = cosZ(k);
    end
    uvec = [mu0^2 mu0 1]';
    v = P*uvec./(Q*uvec);
    a = v(1);
    b = v(2);
    c = v(3);
    estRefl(:,k) = a*radius.^b+c;
end
[X,Y] = ndgrid(sqrt(radius),acosd(cosZ));
surf(X,Y,estRefl,'EdgeColor','none')
zlabel('broadband albedo')
ylabel('solar zenith, degrees')
xlabel('snow radius, {\mu}m')
rval = [50 300 800 1500];
set(gca,'XDir','reverse','XTick',sqrt(rval),'XTickLabel',...
    {num2str(rval(1)),num2str(rval(2)),num2str(rval(3)),num2str(rval(4))})
ylim([0 85])
xlim([sqrt(30) sqrt(1500)])
colorbar
end
