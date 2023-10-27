function yy=taperVector(y,Nd,endval)
    % taper vector using a spline with 2 points
    %input:
    %y - vector, Nx1
    %Nd - number of index values to taper, Nx1
    %yy - tapered vector, Nx1
    yy=y;
    if all(isnan(yy))
        return;% dont do anything if all nans
    else
        yy(end-Nd:end)=endval; %preallocate w/ endvals
        if y(end-Nd) ~= endval && ~isnan(y(end-Nd))
            X=[length(y)-Nd length(y)];
            Y=[0 [y(end-Nd) endval] 0];
            pp=spline(X,Y);
            yy(X(1):X(2))=fnval(pp,(X(1):X(2)));
        end
    end
end