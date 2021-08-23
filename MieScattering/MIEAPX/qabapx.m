function Qabs = qabapx( cRefIn, SizPar, Accur)
% Qabs = qabapx( cRefIn, SizPar, Accur)
% Mie absorption efficiency from asymptotic formula (geometric optics plus
% edge corrections)
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
%   Qabs - Efficiency factor for absorption (cross-section divided by
%   projected area of particle)
%
% ACCURACY : Better than 1 per cent for size parameter > 50-100, compared to
% smoothed Mie results. On the other hand, comparing to unsmoothed Mie
% results which are full of ripple spikes often shows much larger errors.
% Thus, tabular comparisons with mie theory can be quite misleading; only
% graphs can reveal how good these asymptotic approximations really are.

global CN GAM XX

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

Qabs = serialLoop(cRefIn,SizPar,Accur);
    function Q = serialLoop(C,x,A)
        zeroab = 1e-10;
        N = numel(C);
        Q = zeros(N,1);
        % check for near-zero absorption
        tz = imag(C).*x < zeroab;
        if any(tz)
            Q(tz) = 0;
            kN = find(~tz);
            if iscolumn(kN)
                kN = kN.';
            end
        else
            kN = 1:N;
        end
        
        for n=kN % (bypassing zero absorption)
            % copy to globals, passed to various integrands
            CN = C(n);
            XX = x(n);
            GAM = (2/XX)^(1/3);
            % geometric optics integral
            geomop = quadgk(@qagomi,0,pi/2,'RelTol',A);
            % edge integrals
            beup = 1/GAM;
            aemax = (real( CN )-1)/(0.5*GAM^2);
            aeup = min(aemax,6);
            edge = edgint(@qaaemi,@qabemi,aeup,beup,A);
            
            Q(n) = geomop+0.5*GAM^2*edge;
        end
    end
end