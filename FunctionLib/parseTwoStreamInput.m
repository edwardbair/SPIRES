function S = parseTwoStreamInput(R)
% S = parseTwoStreamInput(R)
% parse inputs for the twostream codes, values returned in structure S

% angles, omega, and g must be same size
% direct or diffuse reflectance?
if isempty(R.mu0)
    S.Direct = false;
    [S.omega,S.g] = checkSizes(R.omega,R.g);
else
    S.Direct = true;
    [S.mu0,S.omega,S.g] = checkSizes(R.mu0,R.omega,R.g);
end

% semi-infinite?
if (isscalar(R.tau) && isinf(R.tau)) ||...
        (~isscalar(R.tau) && all(isinf(R.tau(:))))
    S.semiInfinite = true;
else
    S.semiInfinite = false;
    [S.tau0,S.R0,~] = checkSizes(R.tau,R.R0,S.omega);
end
S.method = R.method;
if ~isempty(R.airmass)
    assert(~isempty(R.mu0),'if ''airmass'' is specified, so must cosine of incidence angle')
    [S.pathLength,~] = checkSizes(R.airmass,S.omega);
end

% if inputs are matrices, convert to vectors but save original size
if size(S.omega,1)~=1 && size(S.omega,2) ~=1
    S.origSize = size(S.omega);
    fn = fieldnames(S);
    for k=1:length(fn)
        if ~islogical(S.(fn{k}))
            S.(fn{k}) = S.(fn{k})(:);
        end
    end
end
end
    