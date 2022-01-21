outloc='/tmp/MODISTopography';
tilelist=HMAtiles;
for i=1:length(tilelist)
    fprintf('working on %s\n',tilelist{i});
    createTopoHorizonsMODIS(GMTED.Z,GMTED.R,tilelist{i},outloc);
    fprintf('done w/ %s\n',tilelist{i});
end