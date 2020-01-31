function compareSnow( F,S )
%compare different snow fractions

scrsiz = get(0,'ScreenSize')
figure('Position',[1,1,scrsiz(3),scrsiz(4)])
for k=1:size(F.snow_fraction,3)
    subplot(1,3,1)
    imagesc(F.snow_fraction(:,:,k))
    freezeColors('nancolor',[1 1 1])
    colorbar
    axis equal tight
    subplot(1,3,2)
    imagesc(S.snow_fraction(:,:,k))
    freezeColors('nancolor',[1 1 1])
    colorbar
    axis equal tight
    subplot(1,3,3)
    imagesc(F.snow_fraction(:,:,k)-S.snow_fraction(:,:,k))
    freezeColors('nancolor',[1 1 1])
    colorbar
    axis equal tight
    pause(.5)
end

end