function Q = MieStuff(radius,lambda,N,nm,useParallel)
% Q = MieStuff(radius,lambda,N,nm,useParallel)
%front end to MatScat
% radius, in meters
% lambda, in meters
% N complex index of refraction
% nm index of refraction in medium (e.g. air, water)

% if not scalars, make into vectors
if isscalar(radius) && isscalar(N) && isscalar(lambda)
    % not using S or ANG (1st & 3rd arg)
    [~,C,~,an,bn] = calcmie(radius,N,nm,lambda,180);
    X = getEfficiencies(C, radius, 3);
    % rename just to go into the M structure
    Q.Qext = X.ext;
    Q.Qsca = X.sca;
    Q.Qabs = X.abs;
    Q.omega = X.sca/X.ext;
    Q.g = getAsymmetry(2*pi*radius/lambda,Q.Qsca,an,bn);
    Q.Qpr = X.ext-Q.g*X.sca;
else
    sizeProblem = size(radius);
    radius = radius(:);
    lambda = lambda(:);
    N = N(:);
    Qext = zeros(size(radius));
    Qsca = zeros(size(radius));
    Qabs = zeros(size(radius));
    g = zeros(size(radius));
    if useParallel
        parfor r=1:length(radius)
            [~,C,~,an,bn] = calcmie(radius(r),N(r),nm,lambda(r),180);
            X = getEfficiencies(C, radius(r), 3);
            Qext(r) = X.ext;
            Qsca(r) = X.sca;
            Qabs(r) = X.abs;
            g(r) = getAsymmetry(2*pi*radius(r)/lambda(r),X.sca,an,bn);
        end
    else
        for r=1:length(radius)
            [~,C,~,an,bn] = calcmie(radius(r),N(r),nm,lambda(r),180);
            X = getEfficiencies(C, radius(r), 3);
            Qext(r) = X.ext;
            Qsca(r) = X.sca;
            Qabs(r) = X.abs;
            g(r) = getAsymmetry(2*pi*radius(r)/lambda(r),X.sca,an,bn);
        end
    end
    % dirty fix here because of error in calcmie that occasionally
    % produces NaN - will need to use Wiscombe MIEV0 instead
    if any(isnan(g))
        if nnz(isnan(g))==numel(g)
            error('%d NaN values out of %d, cannot fix',...
                nnz(isnan(g)),numel(g))
        end
        warning('%d NaN values, fixed with mean',nnz(isnan(g)))
        t = isnan(g);
        Qext(t) = nanmean(Qext(:));
        Qsca(t) = nanmean(Qsca(:));
        Qabs(t) = nanmean(Qabs(:));
        g(t) = nanmean(g(:));
    end
    omega = Qsca./Qext;
    Qpr = Qext-g.*Qsca;
    Q.Qext = reshape(Qext,sizeProblem);
    Q.Qabs = reshape(Qabs,sizeProblem);
    Q.Qsca = reshape(Qsca,sizeProblem);
    Q.Qpr = reshape(Qpr,sizeProblem);
    Q.omega = reshape(omega,sizeProblem);
    Q.g = reshape(g,sizeProblem);
end
end