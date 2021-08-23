function Q = qaaemi( x )
% Q = qaaemi( x )
% Above-edge integrand for absorption efficiency

global CN GAM XX

sq = @(c) real(c).^2+imag(c).^2;

illcon = 1e-4;

sinth = 1.+0.5*GAM^2*x;
csithp = sinth/CN;
costhp = sqrt(1-csithp.^2);
z = GAM*cfitf1(sqrt(x));
zc = conj(z);
cw = CN*costhp;
crm22 = (zc-cw)./(z+cw);
cre22 = (zc*CN-costhp)./(z*CN+costhp);

% special case, high absorption
if imag(CN)*XX>5
    Q = 2-sq(crm22)-sq(cre22);
    
else % normal or low absorption
    cnre = real(CN);
    cnim = imag(CN);
    
    if cnim<illcon % well conditioned approximation
        fact = cnim*sqrt(1-(sinth/cnre).^2);
    else
        zeta = log(sq(costhp+complex(0,1)*csithp))/2;
        fact = imag(cw)-zeta.*sinth;
    end
    
    arg = 4*XX*fact;
    e = exp(-arg);
    if abs(arg)<1.e-2
        oneme = arg.*(1+arg.*(-0.5+arg/6));
    else
        oneme = 1-e;
    end
    crm11 = (z-cw)./(z+cw);
    cre11 = (z*CN-costhp)./(z*CN+costhp);
    Q = oneme.*((1-sq(crm22))./(1-e.*sq(crm11))+(1-sq(cre22))./(1-e.*sq(cre11)));
    
end