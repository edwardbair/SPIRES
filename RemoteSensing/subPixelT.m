function [ Ttarget,varargout ] = subPixelT( Tbright,lambda,varargin )
%Ttarget = subPixelT( Tbright,lambda,'Tback',value,'fraction',value )
%[Ttarget,fraction] = subPixelT([Tb1 Tb2],[lam1 lam2],'Tback',value)
%[Ttarget,fraction,Tback] = subPixelT([Tb vector],[lam vector])
%[Ttarget,fraction,Tback] = subPixelT(TbMatrix,[lam1 lam2])
%
%For any option, emissivity may be specified 'emissivity', value,
%specified as a matrix Nx2 for N wavelengths, where the first column are
%the emissivities of the target material and the second column are the
%emissivities of the background.
%
%Solves for sub pixel temperature from thermal IR data.
%
%For all options, Tbright and Ttarget are in degrees Kelvin, lambda is in
%micrometers, and fractions are 0-1 (not percentages)
%
%In a pixel large enough to have different surfaces  of different
%temperatures, for example a pixel with a small fire or a pixel that is
%partly covered with snow, the pixel's infrared brightness temperature is a
%combination of the temperature of a "target" and that of the "background"
%modulated by the fraction of the pixel covered by the target.
%
%From the Planck equation. at any wavelength lambda
%   planck(lambda,Tbright) = fraction*emissivty(target)*planck(lambda,Ttarget) +
%       (1-fraction)*emissivity(background)*planck(lambda,Tback)
%
%Therefore one wants to solve for the target temperature. Because the
%Planck equation is nonlinear, in most situations analytic solutions are
%not possible.
%
%A variety of situations exist, depending on what is known a priori. In the
%case where the background and fraction are known, for example fractional
%snow cover and snow temperature, the solution for the local non-snow
%temperature is possible with just one thermal band
%   Ttarget = subPixelT(Tbright,lambda,'Tback',SnowTemp,'fraction',1-fSCA)
%If vectors of Tbright and lambda are given, then a least-squares solution
%is provided.
%
%With two or more wavelengths and two unknowns, for example a fire of
%unknown size and temperature but with the background temperature known,
%   [Tfire,fraction] = subPixelT([Tb1 Tb2],[lam1 lam2],'Tback',value)
%or if more than 2 wavelengths
%   [Tfire,fraction,Tback] = subPixelT([Tb vector],[lam vector])
%
%If there are just two wavelengths, but there are 2 or more pixels where
%Ttarget and Tback are the same, for example clouds over ocean, then
%   [T]target,fractions,Tback] = subPixelT([Tb matrix],[lam1 lam2])
%where the matrix is Nx2, N pixels 2 wavelengths

% parse input and figure out which scenario we're dealing with
p = inputParser;
addRequired(p,'Tbright',@isnumeric);
addRequired(p,'lambda',@isnumeric);
addParameter(p,'Tback',[],@isnumeric);
addParameter(p,'fraction',[],@isnumeric);
addParameter(p,'emissivity',1,@isnumeric);
parse(p,Tbright,lambda,varargin{:});
Tback = p.Results.Tback;
fraction = p.Results.fraction;
if isequal(size(p.Results.emissivity),[length(lambda) 2])
    emissivity = p.Results.emissivity;
else
    emissivity = repmat(p.Results.emissivity,length(lambda),2);
end

% check inputs
minTempK = 60;
assert(all(Tbright(:)>minTempK),'brightness temperatures must be in Kelvin');
assert(all(lambda(:)>2) && all(lambda(:)<50),'lambda must be in micrometers')
if ~isempty(Tback)
    assert(all(Tback(:)>minTempK),'Tback must be in Kelvin');
end
if ~isempty(fraction)
    assert(all(fraction(:)>0) && all(fraction(:)<1),...
        'fractions must be >0 and <1')
end

% convert wavelength to meters
lambda = lambda*1.e-6;

% which scenario
if isscalar(Tbright) && isscalar(lambda) &&...
        length(fraction)==1 && length(Tback)==1
    scenario = 1; % nonlinear equation of 1 variable
elseif isempty(fraction) && ~isempty(Tback) &&...
        numel(Tbright(:))>=2 && numel(lambda(:))>=2
    scenario = 2; % nonlinear simultaneous equation of 2 variables
end

switch scenario
    case 1 % nonlinear equation of 1 variable
        % Tback and fraction of pixel at Ttarget are known
        Ttarget = solveEq1Var(Tbright,lambda,Tback,fraction,emissivity);
        if nargout>1
            varargout{1} = fraction;
        end
    case 2 % nonlinear simultaneous equations of 2 variables, possibly overdetermined
        % this is the fire scenario in Dozier 1981 doi 10.1016/0034-4257(81)90021-3
        % Tback is known, but Ttarget and its fraction are not
        Tbright = Tbright(:);
        lambda = lambda(:);
        [Ttarget,fraction] = solveEq2Var(Tbright,lambda,Tback,emissivity);
        if nargout>1
            varargout{1} = fraction;
        end
end

end

function T = solveEq1Var(Tb,lambda,Tback,fraction,emissivity)
% start value from Rayleigh-Jeans equation
T0 = (Tb+(fraction-1)*Tback)/fraction;
% solve
[T,fval,exitflag,output] = fzero(@f1Var,T0);
% print warnings if errors
if exitflag~=1
    warning('fzero return value %d',exitflag)
    disp(fval)
    disp(output)
end

    function z=f1Var(T)
        z = planck(lambda,Tb)-(fraction*emissivity(1)*planck(lambda,T)+...
            (1-fraction)*emissivity(2)*planck(lambda,Tback));
    end
end

function [T,f] = solveEq2Var(Tb,lambda,Tback,emissivity)
assert(length(Tb)==length(lambda),'lengths of Tb and lambda must be same')
assert(length(Tb)>1,'need at least 2 values of Tb and lambda')
% starting guess
X0 = [mean(Tb)+(mean(Tb)-Tback) .5];

if length(Tb)==2 % exact solution possible
    options = optimoptions('fsolve','Display','off');
    [X,fval,exitflag,output] = fsolve(@f2Var,X0,options);
    if exitflag<1
        warning('exitflag %d from fsolve',exitflag)
        disp(fval)
        disp(output)
    end
else % overdetermined, solve for minimum
    % upper and lower bounds
    lb = [0;0];
    ub = [Inf;1];
    % options - display off unless exitflag not 1 or 2
    options = optimoptions('fmincon','Display','off');
    [X,fval,exitflag,output] = fmincon(@f2VarMin,X0,[],[],[],[],lb,ub,[],options);
    if exitflag<1 || exitflag>2
        warning('exitflag %d from fmincon',exitflag)
        disp(fval)
        disp(output)
    end
end
T = X(1);
f = X(2);

    function z=f2Var(x)
        Temp = x(1);
        fraction = x(2);
        z = zeros(size(Tb));
        for k=1:length(Tb)
            z(k) = planck(lambda(k),Tb(k))-...
                (fraction*emissivity(k,1)*planck(lambda(k),Temp)+...
                (1-fraction)*emissivity(k,2)*planck(lambda(k),Tback));
        end
    end

    function f=f2VarMin(x)
        z = f2Var(x);
        f = sum(z.^2);
    end
end