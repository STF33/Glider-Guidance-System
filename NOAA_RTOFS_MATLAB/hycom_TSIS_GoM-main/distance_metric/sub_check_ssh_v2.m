% Check ssh/LC contours
% for MHD calculation
% for modified code of mhd_LCLCEcntr_OSSEfcst
% with all possible combinations of the LCEs
%
function sub_check_ssh_v2(xn,yn,xh1,yh1,LCEN,LCE,MHD_LCLCE,icomb,ihc,irc,...
                       flg_lce1,HH,LON,LAT,fnmb);

Nlc1 = sub_numbLCE(LCEN,irc);
Nlc2 = sub_numbLCE(LCE,ihc); % # of LCE on this date 

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

% LCEs - best combination of LCEs
for ilce=1:Nlc2
	nXY = length(LCE(ilce).XY);  % how many records in LCE# ilce
	if ihc>nXY; continue; end;
	xef = LCE(ilce).XY(ihc).X;
	yef = LCE(ilce).XY(ihc).Y;
  if isempty(xef); continue; end;
	Neddy=j0-1;
  plot(xef,yef,'c.');
%keyboard
end

[EComb, icomb] = sub_lce_comb(Nlc2);

jlce = j0;
[Xhc,Yhc,flg] = sub_combineLCLCE_EComb(xh1,yh1,LCE,ihc,EComb,jlce);
plot(Xhc,Yhc,'b.');

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
