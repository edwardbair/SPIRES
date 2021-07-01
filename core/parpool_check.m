function pool=parpool_check(varargin)
%opens up a parpool
%optional input:
% poolsize - number of workers
% with no arg given poolsize=number of cores
p=inputParser;
defaultPoolSize=feature('NumCores');
addOptional(p,'poolsize',defaultPoolSize,@isnumeric);
parse(p,varargin{:});
poolsize=p.Results.poolsize;
q=gcp('nocreate');
distcomp.feature('LocalUseMpiexec',false);

%1)if there's no pool, create one.
%2) if there's an existing pool 
% but it's the wrong size delete it and open a new pool
%3) if there's an existing pool of the correct size do nothing
pool=[];
if ~isempty(q)
    if q.NumWorkers ~= poolsize
        delete(q)
        pool=openpool(poolsize);
    elseif q.NumWorkers == poolsize
        pool=q;
    end
else
    pool=openpool(poolsize);
end
%open new pool
    function pool = openpool(poolsize)
        try
           pool = parpool(poolsize);
        catch
            warning('not able to open parallel pool')
        end
    end
end