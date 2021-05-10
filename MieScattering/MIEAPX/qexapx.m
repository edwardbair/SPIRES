function Qext = qexapx( cRefIn, SizPar, nTerms)
% Qext = qexapx( cRefIn, SizPar, nTerms)
% complex angular momentum approximation to Mie extinction efficiency
%
% Input
%   cRefIn - Complex refractive index (imaginary part may have either sign)
%   SizPar - Mie size parameter (circumference/wavelength)
%   nTerms - number of terms in CAM expansion (3 enough for 1% accuracy, 6
%       is the max)
%
% Ouput
%   Qext - Efficiency factor for extinction (cross-section divided by
%       projected area of particle)
%
% ACCURACY : Better than 1 per cent for size parameter > 20, compared to
% smoothed Mie results. On the other hand, comparing to unsmoothed Mie
% results which are full of ripple spikes often shows much larger errors.
% Thus, tabular comparisons with mie theory can be quite misleading; only
% graphs can reveal how good these asymptotic approximations really are.

% Fock coefficients, Lisle et al. (1985) Phys Rev Lett 55, 555
eipid3 = complex(0.5,sqrt(3)/2);
cm = [1.25512456*eipid3; 0.5322907*eipid3^2; 0.067717];

% check inputs
assert(all(SizPar(:)>=20),'minimum SizPar = %g, must be >=20',min(SizPar(:)))
assert(all(real(cRefIn(:))>1),'minimum real(cRefIn) = %g, must be >1',...
    min(real(cRefIn(:))))
assert(nTerms>=2 && nTerms<=6,'nTerms = %f, must be >=2 and <=6',nTerms)

% arguments into local variables
gam = (2./SizPar ).^( 1/3);
neg = imag(cRefIn)<0;
cn = cRefIn;
if any(neg)
    cn(neg) = conj(cRefIn(neg));
end
    
xx  = SizPar;

% all terms with powers of gam are edge ray contributions
Qext = 2+2*real(cm(1)).*gam.^2;
if nTerms==2
    return
end
cnsq = cn.^2;
capm = sqrt(cnsq-1);
Qext = Qext+imag((cnsq+1)./capm).*gam.^3;

% axial ray contribution
cnm1 = cn-1;
ctemp = 1+complex(0,0.5./xx).*(1./cnm1-cnm1./cn);
multiray = imag(cn).*xx<=4;
if any(multiray)
    % multiple reflected axial rays contribute
    z = (cnm1./(cn+1)).^2.*exp(complex(0,4)*cn.*xx);
    cc0 = cnm1/2;
    cphi = complex(0,0);
    zpowj = complex(1,0);
    for k=1:10
        zpowj = z.*zpowj;
        cterm = zpowj./(k-cc0);
        cphi = cphi+cterm;
        if all(abs(cterm./cphi)<1.e-6)
            break
        end
    end
    ctemp(multiray) = ctemp(multiray)-cc0(multiray).*cphi(multiray);
end

cc1 = (cnsq./((cn+1).*(cnsq-1))).*exp(complex(0,2).*cnm1.*xx);
Qext = Qext-(8./xx).*imag(cc1.*ctemp);

% end axial ray contribution

if nTerms>3
    Qext = Qext+(16/15)*real(cm(2)).*gam.^4;
end
if nTerms>4
    Qext = Qext-(gam.^5/6).*imag(cm(1).*((cnsq+1)./capm.^3).*(3+cnsq.*(-6+2*cnsq)));
end
if nTerms>5
    Qext = Qext+gam.^6.*((8/175)*real(cm(3))-397/480+...
        0.25*real((-1+cnsq.*(-1+cnsq^2))./capm.^2));
end

end