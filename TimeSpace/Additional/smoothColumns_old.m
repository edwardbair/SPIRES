function newColumns = smoothColumns(columns,weight)

% limits to stay within
limits = [nanmin(columns(:)) nanmax(columns(:))];

% some columns all NaN or all zero or all weight = 0, leave those as is
wsum = sum(weight,1);
t = ~isnan(columns);
nsum = sum(t,1);
tokay = wsum>0 & nsum>0;
tzero = nansum(columns,1)==0;

Z = double(columns);
nCols = size(Z,2);
columnLength = size(Z,1);
zzero = zeros(columnLength,1);
parfor c=1:nCols
    if tzero(c) % all zeros
        zhat = zzero;
    elseif tokay(c)
        y = double(columns(:,c));
        t = ~isnan(y);
        if nnz(y(t))>columnLength/4;
            x = (1:columnLength)';
            F = fit(x(t),y(t),'smoothingspline');
            if any(~t) % smoothed original values, interpolated
                G = griddedInterpolant(x(t),F(x(t)),'pchip','nearest');
                zhat = G(x);
            else % no NaNs in original, so just interpolate with the spline
                zhat = F(x);
            end
        else
            zhat = Z(:,c);
            zhat(~t) = 0;
        end
    else
        zhat = Z(:,c);
    end
    Z(:,c) = zhat;
end
newColumns = cast(truncateLimits(Z,limits),'like',columns);
end