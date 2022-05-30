function output_cube=taperDataCube(cube,endvals,mask)
%taper data cube to endvals using a spline, only start and end are used
%all other values are ignored
%input:
%cube - m x n x d cube
%endvals - mxn set of endvals
%mask - skip false values
  % (MATLAB is column-major order)
        N=size(cube);
        cube = reshape(cube,N(1)*N(2),N(3))';
        mask = reshape(mask,N(1)*N(2),1)';
        endvals = reshape(endvals,N(1)*N(2),1)';

        % boundaries of the mask
        fcol = find(mask,1,'first');
        if isempty(fcol) % if the mask is all false
            output_cube = reshape(cube',N(1),N(2),N(3));
            return
        end
        lcol = find(mask,1,'last');
        iCube = double(cube(:,fcol:lcol));
        mask = mask(fcol:lcol);
        endvals = endvals(fcol:lcol);
        limits = [nanmin(iCube(:)) nanmax(iCube(:))];
        
        % by column
        sCube = nan(size(iCube));
        x = (1:size(iCube,1))';
        
        k=gcp('nocreate');
        if ~isempty(k)
            for c=1:size(iCube,2)
                if mask(c)
                    y = iCube(:,c);
                    ev = endvals(c);
                    t=~isnan(y);
                    y=y(t);
                    if length(y)<2 || (y(1)-ev)==0
                        sCube(:,c) = zeros(size(x));
                    else
                        pp=spline([x(1) x(end)] ,[0 [y(1) ev] 0]);
                        sCube(:,c)=fnval(pp,x(1:end));
                    end
                end
            end
        else
            error('start a parpool');
        end
        
        % back into original
        cube(:,fcol:lcol) = cast(truncateLimits(sCube,limits),'like',cube);
        
        output_cube = reshape(cube',N(1),N(2),N(3));