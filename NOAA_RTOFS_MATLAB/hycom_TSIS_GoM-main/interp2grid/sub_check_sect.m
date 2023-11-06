function sub_check_sect(F,Fi,ZM,ZMi,Hs);
%
% Plot 2 sections for checking
%
%
%keyboard
[ll,msb]=size(F);
[xx,dmm]=meshgrid([1:msb],[1:ll]);
figure(11); clf;
pcolor(xx,ZM,F); shading flat;
colorbar
caxis([min(min(F)) max(max(F))]);

%keyboard

[ll,msb]=size(Fi);
[xxi,zzi]=meshgrid([1:msb],ZMi);
figure(12);clf;
pcolor(xxi,zzi,Fi); shading flat;
hold on;
plot(xxi,Hs,'r-');
colorbar
%caxis([5 25]);
caxis([min(min(F)) max(max(F))]);
title('Interpolated into z-layers');


return