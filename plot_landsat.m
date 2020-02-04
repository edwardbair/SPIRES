function plot_landsat(Rdir,out,subset)
    % plot up landsat results;
    %input - Rdir- directory where OLI SR bands 1-7 live
    %subset, empty [] for none or [row1 row2;col1 col2]
    %out - output struct from run_scagd_lands sat w/ fields
    %fsca, grainradius, dust
    R=getOLIsr(Rdir,[]);
    if ~isempty(subset)
        Rsub=zeros([size(out.fsca) size(R.bands,3)]);
        for i=1:size(R.bands,3)
            Rsub(:,:,i)=R.bands(subset(1,1):subset(1,2),...
            subset(2,1):subset(2,2),i);
        end
        R.bands=Rsub;
    end
    
    f1=figure('Position',[0 0 1500 800],'Color',[0.6 0.6 0.6]);
    ha=tight_subplot(2, 2, [0.02 0.01], 0.02, [0 0.03]);
    
    for j=1:4
        axes(ha(j));
        ax=gca;
        
        if j==1
            xx=squeeze(R.bands(:,:,[3 2 1]));
        elseif j==2
            xx=out.fsca;
        elseif j==3
            xx=out.grainradius;
        elseif j==4
            xx=out.dust;
        end        
        if j==1
            image(xx);
        else
            imagesc(xx);
        end
        axis image;
        
        ax.XAxis.Color = 'none';
        ax.YAxis.Color = 'none';
        set(ax,'XTick',[],'YTick',[],'YDir','reverse',...
            'Color',[0.6 0.6 0.6]);
%         freezeColors('nancolor',[0.6 0.6 0.6]);
        
        if j>1
            c=colorbar('Location','EastOutside','Color','w');
            c.Label.Color=[1 1 1];
            c.FontSize=15;
            if j==2
                c.Label.String='fsca, canopy adj. ';
                caxis([0 1]);
            elseif j==3
                c.Label.String='grain radius, \mum';
                caxis([30 quantile(xx(:),0.95)]);
            elseif j==4
                c.Label.String='dust conc, ppmw';
                caxis([0 max(quantile(xx(:),0.95),1)]);
            end
        end 
    end