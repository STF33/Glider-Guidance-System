function sub_plot_TSxsct(fgn,cmp,XX,ZZ,T,S,cl1,cl2,ylim2);

figure(fgn); clf;
set(gcf,'Position',[1573         560         949         751]);
axes('Position',[0.08 0.2 0.85 0.7]);

pcolor(XX,ZZ,T); shading flat;
colormap(cmp);
hold on
caxis([cl1, cl2]);

contour(XX,ZZ,S,[35:0.2:37],'k');
contour(XX,ZZ,S,[36 36],'k','linewidth',1.6);
contour(XX,ZZ,S,[35 35],'Color',[0.2 0.2 0.2],'linewidth',1.6);

set(gca,'Tickdir','out',...
        'Color',[0 0 0],...
        'xlim',[0 max(max(XX))],...
        'ylim',[ylim2 0],...
        'xtick',[0:200:1800],...
        'ytick',[-4000:200:0]);

clb=colorbar('SouthOutside');
set(clb,'Position',[0.1 0.1 0.8 0.025],...
        'Fontsize',13,...
        'Ticks',[0:1:40],...
        'Ticklength',0.025);



return
