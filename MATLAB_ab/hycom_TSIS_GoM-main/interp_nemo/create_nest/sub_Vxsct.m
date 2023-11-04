  function sub_Vxsct(nf,XL,ZZ,Vav,Hs,tst);
%
% Plot vertical x-section of normal velocity
% across Yuc. channel
%

  cl1=flipud(colormap_cold(50));
  cl2=colormap_hot(50);
  cl3=colormap_green(50);
  cl4=colormap_gray(50);
  cmp=[cl1;cl4;cl3;cl2];
  c1=-0.5;
  c2=1.5;
  cnt=[c1:(c2-c1)/length(cmp):c2];


  figure(nf); clf;
  pcolor(XL,ZZ,Vav); shading flat;
  caxis([c1 c2]);
  colormap(cmp);
  hold on;
  plot(XL,Hs,'k-','linewidth',2);

%  contour(XL,ZZ,Vav,[0 0],'k','Linewidth',1.8);
%  contour(XL,ZZ,Vav,[-0.8:0.2:-0.2],'k--','linewidth',1);
%  contour(XL,ZZ,Vav,[0.2:0.2:2],'k-','linewidth',1);

  hght=[];
  lngth=[];
  mint=10;
  fsz=12;
  bxc='w';
%  posc=[0.1 0.045 0.8 0.05];
  posc=[0.92 0.1 0.8 0.08];
  mbx=mint;
  aend=0;
  [az,axc] = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);


  title(tst);
  set(gca,'tickdir','out','fontsize',14);
  set(gcf,'Color',[1 1 1]);






