function Q = qpgomi(theta)
% Q = qagomi(theta)
% Integrand in the geometric optics expression for w*g, where w and g are
% the geometrical optics contributions to Q_sca and the asymmetry factor
% theta - angle in radians, between 0 and pi/

global CN XX

sq = @(c) real(c).^2+imag(c).^2;

illcon = 1e-4;
ci = complex(0,1);

sinth = sin(theta);
costh = cos(theta);
csithp = sinth./CN;
costhp = sqrt(1-csithp.^2);
cw = CN*costhp;
cr1 = (costh-cw)./(costh+cw);
cr2 = (costh*CN-costhp)./(costh*CN+costhp);

if imag(CN)*XX>5
    % specialcase: high absorption
    Q = sinth.*costh.*(2*costh.^2-1).*(-sq(cr1)-sq(cr2));
    
else
    % normal or low absorption
    cnre = real(CN);
    cnim = imag(CN);
    
    if cnim<illcon
        % use well-conditioned approximation
        fact = cnim*sqrt(1-(sinth/cnre).^2);
    else
        zeta = log(sq(costhp+ci*csithp))/2;
        fact = imag(cw)-zeta.*sinth;
    end
    
    arg = 4*XX*fact;
    e = exp(-arg);
    em2ith = complex(2*costh.^2-1,-2*sinth.*costh);
    r1sq = sq(cr1);
    r2sq = sq(cr2);
    cfac = e.*(2*costhp.^2-1+complex(0,2)*csithp.*costhp);
    Q = sinth.*costh.*real(em2ith.*(-r1sq-r2sq+cfac.*(sq(1-cr1.^2)./...
        (1+cfac.*r1sq)+sq(1-cr2.^2)./(1+cfac.*r2sq))));
end
end