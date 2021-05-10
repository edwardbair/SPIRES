function Q = qpbemi( t )
% Q = qpbemi( t )
% Below-edge integrand for w*g

global CN GAM XX

sq = @(c) real(c).^2+imag(c).^2;

illcon = 1e-4;
ci = complex(0,1);

x = t.^2;
cnsq = CN^2;
cimsq = complex(0,1)*(cnsq-1);
cc0 = 2*cnsq-1;
sinth = 1-GAM^2*x/2;
csithp = sinth/CN;
costhp = sqrt((1-csithp).*(1+csithp));
cw = CN*costhp;
cwp = cw+cimsq;
cp = cc0+cimsq*cw;
z = GAM*cfitf2(t);
cq = GAM*complex(x.*t,1/4);

% special case, high absorption
if imag(CN)*XX>5
    zc = conj(z);
    cqc = conj(cq);
    
    Q = 2*t.*(real(conj((cw-zc)./(cw+z)).*(cwp-cnsq.*zc)./(cwp+cnsq.*z))+...
        real(conj((cw-cnsq.*zc)./(cw+cnsq.*z)).*(cwp-cp.*zc)./(cwp+cp.*z))-...
        real(conj((x.*cw-cqc)./(x.*cw+cq)).*(x.*cwp-cnsq.*cqc)./(x.*cwp+cnsq.*cq))-...
        real(conj((x.*cw-cnsq.*cqc)./(x.*cw+cnsq.*cq)).*(x.*cwp-cp.*cqc)./(x.*cwp+cp.*cq)));
    
else % normal or low absorption
    cnre = real(CN);
    cnim = imag(CN);
    
    if cnim<illcon % well conditioned approximation
        fact = cnim*sqrt(1-(sinth./cnre).^2);
        
    else
        zeta = log(sq(costhp+ci*csithp))/2;
        fact = imag(cw)-zeta.*sinth;
    end
    
    arg = 4*XX*fact;
    e = exp(-arg);
    cwm = cw-cimsq;
    cm = cc0-cimsq*cw;
    ce = 2*costhp.^2-1+complex(0,2)*csithp.*costhp;
    
    for ii = 1:2
        zc = conj(z);
        crm11 = -(z-cw)./(z+cw);
        crm22 = (zc-cw)./(z+cw);
        c11m = (cnsq*z-cwm)./(cnsq*z+cwp);
        c22m = (cnsq*zc-cwp)./(cnsq*z+cwp);
        cre11 = -(cnsq*z-cw)./(cnsq*z+cw);
        cre22 = (cnsq*zc-cw)./(cnsq*z+cw);
        c11e = (cm.*z-cwm)./(cp.*z+cwp);
        c22e = (cp.*zc-cwp)./(cp.*z+cwp);
        
        if ii==1
            cff = (1+ci*zc)./(1-ci*z);
        elseif ii==2
            cff = (x+ci*zc)./(x-ci*z);
        end
        
        part{ii} = real(c22m.*conj(crm22))+real(c22e.*conj(cre22))-...
            e.*(real((cff+c22m).*(ce+c11m).*conj((1.+crm11).*(1+crm22))./...
            (1+e.*c11m.*conj(crm11)))+...
            real((cff+c22e).*(ce+c11e).*conj((1+cre11).*(1+cre22))./...
            (1+e.*c11e.*conj(cre11)))); %#ok<AGROW>
        
        
        if ii==1
            z = cq;
            cw = x.*cw;
            cwm = x.*cwm;
            cwp = x.*cwp;
        end
    end
    
    Q = 2*t.*(part{1}-part{2});
end

end