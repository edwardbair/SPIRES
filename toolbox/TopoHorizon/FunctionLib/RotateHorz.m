function [ Ir, rback, hmask ] = RotateHorz( I, varargin )
%[ Ir, rback, mask ] = RotateHorz( I, azimuth, {'nomask'})
%[ Ir ] = RotateHorz(I, rback)
%   rotate elevation grid for horizon computations
%   (uses the Image Processing Toolbox)
%
%INPUT - 1st option, to rotate for computing horizons
%   I - elevation grid to be rotated
%   azimuth - from south, + counter-clockwise
%   'nomask' - optional, default is to create data mask for horizon
%       processing to avoid running the horizon routine outside the grid
%       (this option typically applies when rotating two or more grids with
%       identical geometry)
%OUTPUT - 1st option
%   Ir - rotated grid, with zeros in the areas outside the grid
%   rback - structure needed to reverse rotation
%
%INPUT - 2nd option, to reverse rotation to get horizons back in original
%   coordinates
%   I - grid to be rotated, typically the information about the horizons
%       rather than the elevation grid itself
%   rback - structure with information needed to reverse the rotation
%OUTPUT - 2nd option
%   Ir - grid rotated back to original orientation

% parse the 2nd & 3rd arguments
hmask = [];
assert(~isempty(varargin),...
    'needs 2nd argument, either azimuth angle or structure to reverse rotation')
if isnumeric(varargin{1})
    azm = double(varargin{1}); % rotate the grid for horizon comps
    assert(isscalar(azm),'azimuth must be a scalar')
    forward = true;
elseif isstruct(varargin{1})
    RB = varargin{1};
    forward = false;
else
    error('2nd argument of unknown type, check')
end

if forward
    rback.azimuth = azm;
    if mod(azm,90)==0 % special case, multiple of 90
        rback.multipleof90 = true;
        if azm == 0 % do nothing
            m = 0;
            Ir = I;
        else
            m = -azm/90;
            Ir = rot90(I,m);
        end
        rback.mult = m;
    else % general case, so need to account for padding in order to reverse
        Ir = imrotate(I,-azm,'nearest');
        rback.multipleof90 = false;
        % OK, this is really cheating here, I should figure out how to do
        % it analytically
        mask = true; % default
        if length(varargin)>1
            if strcmpi(varargin{2},'nomask')
                mask = false;
            else
                warning('unknown 3rd argument, setting mask=true') %#ok<WNTAG>
            end
        end
        if mask
            Imask = imrotate(true(size(I)),-azm,'nearest');
            Itry = imrotate(Imask,azm,'nearest','crop');
            [r,c] = find(Itry);
            ur = unique(r);
            uc = unique(c);
            rback.isize = size(I);
            while length(ur)>rback.isize(1)
                ur(1) = [];
                if length(ur)>rback.isize(1)
                    ur(end) = [];
                end
            end
            while length(uc)>rback.isize(2)
                uc(1) = [];
                if length(uc)>rback.isize(2)
                    uc(end) = [];
                end
            end
            rback.row = [min(ur) max(ur)];
            rback.col = [min(uc) max(uc)];
            rback.minmax = [min(I(:)) max(I(:))];
            hmask = Imask;
        end
    end
else % reverse a previous rotation
    rback = []; % make sure we don't return something spurious
    if RB.multipleof90
        if RB.mult==0 % do nothing
            Ir = I;
        else
            Ir = rot90(I,-RB.mult);
        end
    else
        Ir = imrotate(I,RB.azimuth,'nearest','crop');
        Ir = Ir(RB.row(1):RB.row(2),RB.col(1):RB.col(2));
        if ~isequal(size(Ir),RB.isize)
            if isempty(Ir)
                Ir = NaN(RB.isize);
            else
                Ir = imresize(Ir,RB.isize);
            end
        end
        fixmask = Ir<RB.minmax(1) | Ir>RB.minmax(2) | isnan(Ir);
        if nnz(fixmask)
            Ir(fixmask) = mean(RB.minmax);
            Ir = inpaintCoherent(Ir,fixmask);
        end
    end
end

end