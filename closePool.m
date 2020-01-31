function closePool()
% close the MATLAB pool
if ~matlabpool('SIZE')
    matlabpool('close');
end
end