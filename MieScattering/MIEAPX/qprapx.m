function Qpr = qprapx( cRefIn, SizPar, Accur)
% Qpr = qprapx( cRefIn, SizPar, Accur)
% Mie radiation pressure efficiency from asymptotic formula
% (geometric optics plus edge corrections)
%
% Input
%   cRefIn - Complex refractive index (imaginary part may have either sign)
%   SizPar - Mie size parameter (circumference/wavelength)
%   Accur - Desired relative accuracy of result; if ACCUR = 10^(-n), then
%       the result will be accurate to at least n+1 significant digits
%       (usually more).
%   A good choice for production runs is 1.e-1. Values smaller than 1.e-2
%   waste lots of computer time and may cause many 'convergence failure'
%   messages.
%
% Ouput
%   Qpr - Efficiency factor for radiation pressure (cross-section divided by
%       projected area of particle)
%
% ACCURACY: Better than 1 per cent for size parameter > 50, compared to
% smoothed Mie results. On the other hand, comparing to unsmoothed Mie
% results which are full of ripple spikes often shows much larger errors.
% Thus, tabular comparisons with mie theory can be quite misleading; only
% graphs can reveal how good these asymptotic approximations really are.

% check inputs
assert(isequal(size(SizPar),size(cRefIn)),...
    'input vectors SizPar and cRefIn must be same size')
assert(all(SizPar>=10) && all(real(cRefIn)>=1),...
    'SizPar = %f, cRefIn = %f, too small',...
    min(SizPar), min(real(cRefIn)))

% imaginary part positive
if any(imag(cRefIn)<0)
    t = imag(cRefIn)<0;
    cRefIn(t) = conj(cRefIn(t));
end

Qpr = serialLoop(cRefIn,SizPar,Accur);
    function Q = serialLoop(C,x,A)
        global CN XX GAM
        N = numel(C);
        C = C(:);
        x = x(:);
        Q = zeros(N,1);
        for n=1:N
            % copy to globals
            CN = C(n);
            XX = x(n);
            GAM = (2/x(n))^(1/3);
            
            % geometric optics integral
            geomop = quadgk(@qpgomi,0,pi/2,'RelTol',A);
            
            % edge integrals
            beup = 1/GAM;
            aemax = (real(C(n))-1)/(0.5*GAM^2);
            aeup = min(aemax,6);
            edge = edgint(@qpaemi,@qpbemi,aeup,beup,A);
            
            Q(n) = 1-(geomop+0.5*GAM^2*edge);
        end
    end

end