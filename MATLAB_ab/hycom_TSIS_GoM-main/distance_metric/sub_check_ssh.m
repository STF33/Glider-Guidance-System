% Check ssh/LC contours
% for MHD calculation
%
function sub_check_ssh(xn,yn,xh1,yh1,LCEN,LCE,MHD_LCLCE,Nlc1,Nlc2,ihc,irc,...
                       HH,LON,LAT,fnmb,flg_lce1);

Nlc1 = length(LCEN);
Nlc2 = length(LCE);

mhd_fcst = min(min(MHD_LCLCE));
[j0,i0] = find(MHD_LCLCE==mhd_fcst,1);

figure(fnmb); clf;
set(gcf,'Position',[1705         682         801         639]);
hold on;
plot(xn,yn,'r.');
plot(xh1,yh1,'b.');
for ilce=1:Nlc1
	xen = LCEN(ilce).XY(irc).X;
	yen = LCEN(ilce).XY(irc).Y;
	Neddy=i0-1;
	if flg_lce1(ilce) == 1
		plot(xen,yen,'r.');
	else
		plot(xen,yen,'k.');
	end
end
for ilce=1:Nlc2
	nXY = length(LCE(ilce).XY);  % how many records in LCE# ilce
	if ihc>nXY; continue; end;
	xef = LCE(ilce).XY(ihc).X;
	yef = LCE(ilce).XY(ihc).Y;
	Neddy=j0-1;
	if ilce==Neddy
		plot(xef,yef,'b.');
	else
		plot(xef,yef,'c.');
	end
%keyboard
end
contour(LON,LAT,HH,[0 0],'k');
axis('equal');
set(gca,'xlim',[-98 -80],...
				'ylim',[17 31]);


plot(-97,30.5,'r.');
text(-96.8,30.5,'NEMO MHD contours');
plot(-97,30,'k.');
text(-96.8,30,'NEMO not MHD');
plot(-87,19,'b.');
text(-86.5,19,'HYCOM MHD contours');
plot(-87,18,'c.');
text(-86.5,18,'HYCOM not MHD');


return
