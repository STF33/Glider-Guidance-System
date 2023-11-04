% Plot 2 spatial fields 
% To compare interpolation
% indicate NEMO grid box

function sub_plot2flds(nf,LNN,LTN,A01,LON,LAT,AAi,LONN,LATN,stl1,stl2,fldnm,zz0);

ln1=min(min(LONN));
ln2=max(max(LONN));
lt1=min(min(LATN));
lt2=max(max(LATN));

lnmn = min(min(LON));
lnmx = max(max(LON));
ltmn = min(min(LAT));
ltmx = max(max(LAT));

figure(nf); clf;
set(gcf,'Position',[1373         622        1110         709]);
axes('Position',[0.07 0.2 0.4 0.7]);
pcolor(LNN,LTN,A01); shading flat;
hold on;
% NEMO domain:
plot([ln1 ln2],[lt1 lt1],'k-');
plot([ln2 ln2],[lt1 lt2],'k-');
plot([ln1 ln2],[lt2 lt2],'k-');
plot([ln1 ln1],[lt1 lt2],'k-');

c1=0;
c2=32;
if strncmp(fldnm,'saln',4);
		c1=30;
		c2=37;
  if zz0<=-800
    c1=34.8;
    c2=35.5;
  elseif zz0>-800 & zz0<-100
    c1=34.8;
    c2=36.5;
  end
elseif strncmp(fldnm,'temp',4);
		if zz0>=-200
				c1=15;
				c2=30;
		elseif zz0<-200 & zz0>-1000
				c1=7;
				c2=18;
  else
    c1=3;
    c2=6;
		end
end
caxis([c1 c2]);

colorbar('SouthOutside');
axis('equal');
%  set(gca,'xlim',[950 1500],...
%          'ylim',[1060 1350]);
set(gca,'xlim',[lnmn lnmx],...
								'ylim',[ltmn ltmx]);
set(gca,'xtick',[],...
								'ytick',[]);
title(stl1);

axes('Position',[0.53 0.2 0.4 0.7]);
pcolor(LON,LAT,AAi); shading flat;
hold on;
% NEMO domain:
plot([ln1 ln2],[lt1 lt1],'k-');
plot([ln2 ln2],[lt1 lt2],'k-');
plot([ln1 ln2],[lt2 lt2],'k-');
plot([ln1 ln1],[lt1 lt2],'k-');

caxis([c1 c2]);
colorbar('SouthOutside');
axis('equal');
%  set(gca,'xlim',[1 1401],...
%          'ylim',[1 891]);
set(gca,'xlim',[lnmn lnmx],...
								'ylim',[ltmn ltmx]);
set(gca,'xtick',[],...
								'ytick',[])
%          'xlim',[1 ah2],...
%          'ylim',[1 ah1]);
title(stl2);


return
