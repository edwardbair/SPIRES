function mieQvalsTest()
%run Wiscombe's test cases
%

%test values
TestXX = [0.099 0.101 100 10000 0.099 0.101 10 1000 1 100 10000 0.055 0.056 ...
    1 100 10000 1 100 10000]';
TestQE = [3.209674E-04 3.477160E-04 2.008102E+00 2.000289E+00 7.417859E-06 ...
    8.033542E-06 2.232265E+00 1.997908E+00 9.395198E-02 2.101321E+00 ...
    2.004089E+00 1.014910E-01 1.033467E-01 2.336321E+00 2.097502E+00 ...
    2.004368E+00 2.532993E+00 2.071124E+00 2.005914E+00]';
TestQS = [3.209674E-04 3.477160E-04 2.008102E+00 2.000289E+00 7.417859E-06 ...
    8.033542E-06 2.232265E+00 1.997908E+00 9.392330E-02 2.096594E+00 ...
    1.723857E+00 1.131687E-05 1.216311E-05 6.634538E-01 1.283697E+00 ...
    1.236574E+00 2.049405E+00 1.836785E+00 1.795393E+00]';
TestGQ = [-1.275386E-04 -1.381344E-04 1.005911E+00 1.000284E+00 1.074279E-08 ...
    1.211000E-08 2.001164E+00 1.688121E+00 1.733048E-02 1.821854E+00 ...
    1.564987E+00 5.558541E-09 6.193255E-09 1.274736E-01 1.091466E+00 ...
    1.046525E+00 -2.267961E-01 1.021648E+00 9.842238E-01]';
TestCR = complex([0 0 0 0 0.75 0.75 0.75 0.75 1.33 1.33 1.33 1.5 1.5 1.5 1.5 1.5 10 10 10],...
    [0 0 0 0 0 0 0 0 -1e-05 -1e-05 -1e-05 -1 -1 -1 -1 -1 -10 -10 -10]);
TestCR = TestCR(:);

nCases = length(TestXX);

disp('Test of miev0 compared to Wiscombe''s tabulated values in MVTstNew.f')
tstart = tic;
for n=1:nCases
    M = mieQvals(TestCR(n),TestXX(n));
    fprintf('case %d, ref index = (%g,%g), x=%g\n',...
        n,real(TestCR(n)),imag(TestCR(n)),TestXX(n))
    disp(M)
    fprintf('ratios: qext %g, qsca %g, qabs %g, gqsc %g\n\n',...
        ratio(M.Qext,TestQE(n)),ratio(M.Qsca,TestQS(n)),...
        ratio(M.Qabs,TestQE(n)-TestQS(n)),...
        ratio(M.g*M.Qsca,TestGQ(n)));
end
tend = toc(tstart);
fprintf('time for the %d tests = %g sec\n',nCases,tend);
end

function r = ratio(a,b)
small = eps('single');
big = realmax('single');
powmax = log10(big);
powmin = log10(small);
if a==0
    if b==0
        r = 1;
    else
        r = 0;
    end
elseif b==0
    if sign(a)==1
        r = realmax('single');
    else
        a = realmin('single');
    end
elseif abs(a)<small && abs(b)<small
    r = 1;
elseif abs(log10(abs(a))-log10(abs(b)))>=powmax
    r = big;
elseif abs(log10(abs(a))-log10(abs(b)))<=powmin
    r = small;
else
    r = abs(a/b);
end
if sign(a)~=sign(b)
    r = -r;
end
end