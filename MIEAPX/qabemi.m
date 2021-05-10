function Q = qabemi( t )
% Q = qabemi( t )
% Below-edge integrand for absorption efficiency

global CN GAM XX

sq = @(c) real(c).^2+imag(c).^2;

illcon = 1e-4;

x = t.^2;
sinth = 1-GAM^2*x/2;
csithp = sinth/CN;
costhp = sqrt(1-csithp.^2);
cw = CN*costhp;
z = GAM*cfitf2(t);
crm22 = (conj(z)-cw)./(z+cw);
cre22 = (CN*conj(z)-costhp)./((CN*z)+costhp);
crm = ((GAM*t)-cw)./((GAM*t)+cw);
cre = ((CN*(GAM*t))-costhp)./((CN*(GAM*t))+costhp);

% special case, high absorption
if(imag(CN)*XX>5)
    Q = 2*t.*(sq(crm)-sq(crm22)+sq(cre)-sq(cre22));
else % normal or low absorption
    cnre = real(CN);
    cnim = imag(CN);
    %well conditioned approximation
    if cnim<illcon
        fact = cnim*sqrt(1-(sinth/cnre).^2);
    else
        zeta = log(sq(costhp+complex(0,1)*csithp))/2;
        fact = imag(cw)-zeta.*sinth;
    end
    arg = 4*XX*fact;
    e = exp(-arg);
    if abs(arg)<1e-2
        oneme = arg.*(1+arg.*(-0.5+arg/6));
    else
        oneme = 1-e;
    end
    crm11 = (z-cw)./(z+cw);
    cre11 = ((CN*z)-costhp)./((CN*z)+costhp);
    rmsq = sq(crm);
    resq = sq(cre);
    Q = 2*t.*oneme.*((1-sq(crm22))./(1-e.*sq(crm11))+(1-sq(cre22))./...
        (1-e.*sq(cre11))-(1-rmsq)./(1-rmsq.*e)-(1-resq)./(1-resq.*e));
end
end