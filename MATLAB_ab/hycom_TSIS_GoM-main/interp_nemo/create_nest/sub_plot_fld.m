function sub_plot_fld(var,F,zz,lon,lat);
% Check mean fields U,V,T,S

[kk,mm,nn]=size(F);
lr=1;

switch (var),
 case('t'),
  c1=4;
  c2=28;
 case('s'),
  c1=34;
  c2=37;
 case({'u','v'})
  c1=-0.5;
  c2=0.5;
 case('rho');
  c1=23;
  c2=27;
end


% Plot OB:
je=2;
ie=nn-1;
OBs=squeeze(F(:,je,:));
OBe=squeeze(F(:,:,ie));

stt1=sprintf('%s South OB',var);
stt2=sprintf('%s East OB',var);
stt3=sprintf('Layer %i, %s ',lr,var);

figure(10); clf;
axes('Position',[0.07 0.55 0.88 0.4]);
pcolor(lon,zz,OBs); shading flat;
caxis([c1 c2]);
colorbar
set(gca,'xlim',[-88 -78],...
	'ylim',[-5000 0],...
	'color',[0 0 0]);
title(stt1);

axes('Position',[0.07 0.08 0.88 0.4]);
pcolor(lat,zz,OBe); shading flat;
caxis([c1 c2]);
colorbar
set(gca,'xlim',[min(lat) max(lat)],...
	'ylim',[-5000 0],...
	'color',[0 0 0]);
title(stt2);

figure(11); clf;
a=squeeze(F(lr,:,:));
pcolor(a); shading flat;
caxis([c1 c2]);
title(stt3);
colorbar

return