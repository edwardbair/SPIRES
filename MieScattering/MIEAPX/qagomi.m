function Q = qagomi(theta)
% Q = qagomi(theta)
% Integrand in geometric optics part of absorption efficiency
% theta - angle in radians, between 0 and pi/2

global CN XX

sq = @(c) real(c).^2+imag(c).^2;

illcon = 1e-4;
ci = complex(0,1);

sinth = sin(theta);
costh = cos(theta);
csithp = sinth/CN;
costhp = sqrt(1-csithp.^2);
cw = CN*costhp;
cr1 = (costh-cw)./(costh+cw);
cr2 = (costh*CN-costhp)./(costh*CN+costhp);

if imag(CN)*XX>5
    % specialcase: high absorption
    Q = sinth.*costh.*(2-sq(cr1)-sq(cr2));
else
    % normal or low absorption
    cnre = real(CN);
    cnim = imag(CN);
    
    if cnim<illcon
        % use well-conditioned approximation
        tmp = sqrt(1-(sinth/cnre).^2);
        v = cnim./tmp;
        fact = cnim*tmp;
    else
        v = imag(cw);
        zeta = log(sq(costhp+ci*csithp))/2;
        fact = v-zeta.*sinth;
    end
    
    arg = 4*XX*fact;
    e = exp(-arg);
    if abs(arg)<1.e-2
        oneme = arg.*(1+arg.*(-0.5+arg/6));
    else
        oneme = 1-e;
    end
    
    r1sq = sq(cr1);
    r2sq = sq(cr2);
    vt = cnre*imag(costhp)-cnim*real(costhp);
    costsq = costh.^2;
    xtra1 = 16*r1sq.*e.*(v.^2.*costsq)./((costsq-sq(cw)).^2+4*(v.^2.*costsq));
    xtra2 = 16*r2sq.*e.*(vt.^2.*costsq)./...
        ((sq(CN)*costsq-sq(costhp)).^2+4*(vt.^2.*costsq));
    
    Q = sinth.*costh.*oneme.*((1-r1sq-xtra1)./...
        (1-r1sq.*e)+(1-r2sq-xtra2)./(1-r2sq.*e));
end

end