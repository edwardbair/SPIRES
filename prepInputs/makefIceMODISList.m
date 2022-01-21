function makefIceMODISList(S,tilelist,outloc,glims_res)
% run makefIce for a list of MODIS tiles
% S - (global) RGI shapefile as matlab structure
% tilelist - cell vector of tiles to process
% outloc
for i=1:length(tilelist)
    tile=tilelist{i};
    [~,~,RRf]=sinusoidProjMODtile(tile);
    RR=RRf.RasterReference_500m;
    tic;
    fprintf('working on %s\n',tile)
    fice=makefIce(S,RR,glims_res);
    save(fullfile(outloc,sprintf('%s.mat',tile)),'fice','RR');
    toc;
    fprintf('done w %s\n',tile)
end