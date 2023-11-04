% Plot vertical section of U
function sub_plotUsection_Yucatan(Tr,fgn,kb,ctl,xlim0,ylim0,clim1,clim2,dcc);

ZZ = Tr(kb).ZZ;
dX = Tr(kb).Dst;
Vav= Tr(kb).Vav;
Hs = Tr(kb).Hs;

XX = cumsum(dX)*1e-3; % distance, km


figure(fgn); clf;
set(gcf,'Position',[1661         542         806         793]);
axes('Position',[0.1 0.1 0.75 0.8]);
hold on;

% Colormap - asumed abs. negative lim value is 3 times smaller than positive lim value
ncc=100;
%clr1 = colormap_purple(ncc);
cl1=[0.5,0,0.35];
cl2=[1,1,1];
clr1 = mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0,0.3,0.9];
clr2 = mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[0.5,1,0.5];
clr3 = mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.8,0];
clr4 = mix_2colors(cl1,cl2,ncc);

cl1=cl2;
cl2=[1,0.3,0];
clr5 = mix_2colors(cl1,cl2,ncc);


cmp = [clr1;clr2;clr3;clr4;clr5];
cmp = smooth_colormap(cmp,5);

pcolor(XX,ZZ,Vav); shading interp
colormap(cmp);
caxis([clim1 clim2]);
set(gca,'tickdir','out',...
        'xlim',[0 xlim0],...
        'ylim',[ylim0 0],...
        'xtick',[0:25:xlim0],...
        'ytick',[-2600:200:0],...
        'color',[0 0 0]);

contour(XX,ZZ,Vav,[-1:dcc:1.6],'k-','linewidth',1);
contour(XX,ZZ,Vav,[0 0],'k-','linewidth',1.8);


% Draw bottom
Hs=Hs(:);
XX=XX(:);
BZ = [ylim0;Hs;ylim0];
XZ = [XX(1);XX;XX(end)];
patch(XZ,BZ,[0 0 0]);

%keyboard

%axes('Position',[0.9 0.1 0.07 0.8]);
clb=colorbar;
set(clb,'Position',[0.88 0.1 0.025 0.8],...
        'Fontsize',14,...
        'Ticks',[clim1:0.3:clim2],...
        'Ticklength',0.03);

title(ctl);

return
