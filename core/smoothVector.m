function yy=smoothVector(x,y,w,p)
    %smooths a vector using a weighted smoothing spline
    %input: x, e.g, days, locations, Nx1, no nans;
    % y, points to smooth, Nx1, can include nans
    % w, weights, Nx1
    % p, smoothing parameter for spline, scalar
    %output
    % yy, smoothed data points, Nx1
    y=double(y);
    x=double(x);
    w=double(w);

    t=~isnan(y);
    if sum(y,'omitnan') == 0 || nnz(t)<2 || nnz(w>0) < 2  
                yy = zeros(size(y));
    else
        F=fit(x,y,'smoothingspline',...
            'weights',w,'exclude',~t,...
            'SmoothingParam',p);
        yy=F(x);
    end
end