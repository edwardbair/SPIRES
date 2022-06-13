function [good,bad]=checkIntegrityBatch(filelist)
%takes about 3 hr w 30 cores for 7162 hdf files
%using 120 cores - < 1/3rd of that
good=cell(size(filelist));
bad=cell(size(filelist));
parfor i=1:length(filelist)
    try
        GetMOD09GA(filelist{i},'allbands');
        good{i}=filelist{i};
        fprintf('%s passed\n',filelist{i})
           
    catch
        bad{i}=filelist{i};
        fprintf('%s failed\n',filelist{i})
    end
end