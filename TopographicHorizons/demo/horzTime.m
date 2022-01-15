function [timevector,varargout] = horzTime(poolsize,Z,R)
% [timevector [,A,H,D]] = horzTime(poolsize,Z,R)
%time the horizon program for a geographic grid

timevector = zeros(1+2*length(poolsize),3);

m = 1;
for k=1:1+length(poolsize)
    p = k-1;
    if k==1
        tic
        [A,H,D] = horizonAllDirections(Z,R);
        elapsedTime = toc;
        timevector(m,:) = [1 elapsedTime k];
        nout = max(nargout,1) - 1;
        for n=1:nout
            switch n
                case 1
                    varargout{n} = A;
                case 2
                    varargout{n} = H;
                case 3
                    varargout{n} = D;
            end
        end
    else
        parpool(poolsize(p))
        tic
        [A,H,D] = horizonAllDirections(Z,R,'parallel','rotate'); %#ok<ASGLU>
        elapsedTime = toc;
        m = m+1;
        timevector(m,:) = [poolsize(p) elapsedTime k];
        tic
        [A,H,D] = horizonAllDirections(Z,R,'parallel','profile'); %#ok<ASGLU>
        elapsedTime = toc;
        m = m+1;
        timevector(m,:) = [poolsize(p) elapsedTime -k];
        delete(gcp);
    end
end

end