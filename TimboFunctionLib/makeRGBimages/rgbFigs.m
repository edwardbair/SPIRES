function [  ] = rgbFigs( landsatFolder,RsHomeDir,sceneID, outputFolder )
%RGBFIG Summary of this function goes here
%   Make a figure that has false color landsat 8 scene, and colored images
%   of each mask result. Show pixels of snow and cloud that are within 5%
%   of each other in a different color

[ masks, Rs ] = OLIsnowCloudmasks(landsatFolder,RsHomeDir,sceneID,'cirrus');

[ cfRGB, cmapCF ] = maskRGBs(masks.CFmask);
[ bqaRGB, cmapBQA ] = maskRGBs(masks.BQA);
[ rtRGB, cmapRT ] = maskRGBs(masks.RTmask,'cirrus');
[ fRGB ] = landsat8fRGB(Rs);

fig = figure;
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);

%UL - fRBG image
axes(ha(1));
imagesc(fRGB)
freezeColors
set(gca,'XTickLabel',[],'YTickLabel',[])

%UR - RTmask
axes(ha(2));
imagesc(rtRGB)
colormap(cmapRT)
freezeColors
set(gca,'XTickLabel',[],'YTickLabel',[])

%LL - CFmask
axes(ha(3));
imagesc(cfRGB)
colormap(cmapCF)
freezeColors

%LR - BQA - change this to the legend?
axes(ha(4));
imagesc(bqaRGB)
colormap(cmapBQA)
freezeColors
set(gca,'XTickLabel',[],'YTickLabel',[])

% force each plot to be square
axesHandles = findobj(get(fig,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')

dpi=300; 
name=[sceneID '_fullScene_AGU.png'];
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 6];
fig.PaperPositionMode = 'manual';
print(fullfile(outputFolder,name),'-dpng',['-r',num2str(dpi)],'-opengl')
end

