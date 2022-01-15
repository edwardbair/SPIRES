function et = Demo_rotate(Z)
%shows missing points after rotation
%
tic; % start the timer

I = reshape(1:20,5,4);
Ir = imrotate(I,45,'crop');
Ir(Ir==0) = NaN;
% eliminate all-NaN rows and columns
for k=size(Ir,2):-1:1
    if all(isnan(Ir(:,k)))
        Ir(:,k) = [];
    end
end
for k=size(Ir,1):-1:1
    if all(isnan(Ir(k,:)))
        Ir(k,:) = [];
    end
end
% image the original and rotated, labeling the values
figure('Name','Fig. 3 Cells lost in rotation')
subplot(2,2,1)
imagesc(I)
axis equal tight;
set(gca,'XTick','')
set(gca,'YTick','')
[x,y] = meshgrid(1:size(I,2),1:size(I,1));
for k=1:numel(I)
    if I(k)>10
        text(x(k),y(k),num2str(I(k)));
    else
        text(x(k),y(k),num2str(I(k)),'Color',[1 1 1]);
    end
end
subplot(2,2,3)
imagesc(Ir)
% freezeColors('nancolor',[1 1 1])
axis equal tight
set(gca,'XTick','')
set(gca,'YTick','')
[x,y] = meshgrid(1:size(Ir,2),1:size(Ir,1));
for k=1:numel(Ir)
    if ~isnan(Ir(k))
        if Ir(k)>10
            text(x(k),y(k),num2str(Ir(k)))
        else
            text(x(k),y(k),num2str(Ir(k)),'Color',[1 1 1])
        end
    end
end

% compute points lost with each rotation
ang = 0:90;
zlost = zeros(size(ang));
n = 2;
nZ = numel(Z);
for a=ang
    if a>0 && a<90
        Zr = imrotate(Z,a,'crop');
        nR = nnz(~isnan(Zr) & Zr>0);
        zlost(n) = (nZ-nR)/nZ;
        n = n+1;
    end
end
subplot(2,2,[2 4])
plot(ang,zlost,'k','LineWidth',1)
xlabel('angle of rotation')
ylabel('fraction cells lost')
set(gca,'YDir','reverse')
xlim([0 90])
set(gca,'XTick',0:15:90)

et = toc;

fprintf('This code %s reproduces Fig. 3\n',mfilename)

end