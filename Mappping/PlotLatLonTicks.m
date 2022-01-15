function [xt,yt,lontick,lattick,lonLabel,latLabel]=PlotLatLonTicks(Proj,R,N,interval,varargin)
% tick marks for lat/lon
% [xt,yt,lontick,lattick]=PlotLatLonTicks(Proj,R,N,interval)
% Proj projection structure
% R referencing matrix
% N grid size (vector of length 2)
% interval in degrees (scalar or 2-element vector, lat-lon interval)
% variable arguments (1 or 2)
%   'x' supresses plots on x-axis
%   'y' supresses plots on y-axis

% axes to suppress
numArg = 4;
optarg = nargin-numArg;
supressX = false;
supressY = false;
if optarg>0
    for k=1:optarg
        if strcmpi(varargin{k},'x')
            supressX = true;
        elseif strcmpi(varargin{k},'y')
            supressY = true;
        end
    end
end

if length(interval)==1
    latInterval = interval;
    lonInterval = interval;
else
    latInterval = interval(1);
    lonInterval = interval(2);
end

I=(1:N(1))';
J=(1:N(2))';
[x,~]=pix2map(R,N(1)*ones(N(2),1),J);
[~,y]=pix2map(R,I,ones(N(1),1));
[~,lon]=projinv(Proj,x,y(end)*ones(length(x),1));
lontick=lonInterval*ceil(min(lon/lonInterval)):lonInterval:max(lon);
xt=interp1(lon,x,lontick);
[lat,~]=projinv(Proj,x(1)*ones(length(y),1),y);
lattick=latInterval*ceil(min(lat/latInterval)):latInterval:max(lat);
yt=interp1(lat,y,lattick);
% lattick and lontick must be strings instead
lonLabel = cell(1,length(lontick));
latLabel = cell(1,length(lattick));
for k=1:length(lonLabel)
    if lontick(k)>0
        lonLabel{k} = [num2str(lontick(k)) 'E'];
    else
        lonLabel{k} = [num2str(-lontick(k)) 'W'];
    end
end
for k=1:length(latLabel)
    if lattick(k)>0
        latLabel{k} = [num2str(lattick(k)) 'N'];
    else
        latLabel{k} = [num2str(-lattick(k)) 'S'];
    end
end
if supressX
    lonLabel = '';
end
if supressY
    latLabel = '';
end
set(gca,'XTick',xt,'YTick',yt,'XTickLabel',lonLabel,...
    'YTickLabel',latLabel,'FontName','Calibri')
end