% Plot RMSE maps
%
function sub_plot_rmse_v2(nfg,HH,LON,LAT,Iocn,rmse,stl,c1,c2);

if isempty(c1) | isempty(c2)
  c1=0;
  c2=0.5;
end
nint=250;
CMP=create_colormap7(nint,c1,c2);
cmp=CMP.colormap;
cnt=CMP.intervals;

nk=round(length(cnt)*0.1);
for ik=1:nk
  cmp(ik,:)=[1 1 1];
end
cmp = smooth_colormap(cmp,23);
cmp = smooth_colormap(cmp,23);
cmp = smooth_colormap(cmp,23);



A=HH*nan;
A(Iocn)=rmse;

Lmsk = HH*0;
Lmsk(HH<0)=1;
lcmp=[0 0 0; 1 1 1];


figure(nfg); clf;
set(gcf,'Position',[252 466 1036 805]);
axes('Position',[0.09 0.22 0.86 0.7]);
hold on;
pcolor(LON,LAT,Lmsk); shading flat;
colormap(lcmp);
caxis([0 1]);
freezeColors;

axis('equal');
set(gca,'tickdir','out',...
    'xlim',[-98 -80],...
    'xtick',[-100:2:-70],...
    'ylim',[18 31],...
    'ytick',[18:34],...
    'Fontsize',12);

pcolor(LON,LAT,A); shading flat;
caxis([c1 c2]);
colormap(cmp);

contour(LON,LAT,HH,[-200 -200],'k-','Color',[0.4 0.4 0.4],...
        'Linewidth',1.6);

clb=colorbar('SouthOutside');
set(clb,'Position',[0.18 0.1 0.66 0.025],...
        'Fontsize',12,...
        'Ticklength',0.025);

title(stl);




return
