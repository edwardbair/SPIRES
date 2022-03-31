function [Rc,c]=normalizeReflectance(R,slope,aspect,solarZ,solarAzimuth,varargin)
%perform topographic correction on reflectances such that all surfaces are
%normalized to a flat and level surface
%input:
%R - reflectance cube, m x n x z - w/ m x n as height and width and z as band #
%slope and aspect, e.g. from consolidateTopography
% optional: model type 'cosine', e.g. eq 1&2 below or 'rotation'
% 'rotation', e.g. 
% output:
%corrected reflectance
% Reference: Tan et al 2013 http://dx.doi.org/10.1016/j.rse.2013.05.013
% c - correction factor

%b/c of Dozier's strange convention for aspect
if azimuthPreference
    aspect = 180 - aspect;
end

model='cosine';
if nargin==6
    model=varargin{1};
end


%Tan et al eqs 1&2
ic=(cosd(solarZ).*cosd(slope)+sind(solarZ).*...
            sind(slope).*cosd(solarAzimuth-aspect));

switch model
    case 'cosine'     
    c=cosd(solarZ)./ic;
    Rc=R.*c;
    %fix negative values
    Rc(Rc<0)=0;
    % scale spectra with max > 1 by max
    Rc=scaleMultiBandCube(Rc);
    case 'rotation'
        Rc=NaN(size(R));
        sz=size(Rc);
        Rc=reshape(Rc,[sz(1)*sz(2) sz(3)]);
        for i=1:size(R,3) %band dependent, Tan et al eq 5
            B=R(:,:,i);
            B=B(:);
            A=[ic(:) ones(size(ic(:)))];
            t=~isnan(ic(:)) & ~isnan(B);
            x=A(t,:)\B(t);
            y=B-x(1).*(ic(:)-cosd(solarZ(:)));
            y(y<0)=0;
            y(y>1)=1;
            Rc(:,i)=y;
        end
        Rc=reshape(Rc,sz);
end

end

    
    
