% Anlysis LC  of NEMO
% LC and LCE extracted in 
% extr_lc_ssh_nemo.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_calc=0;  %=1 - calc shortest dist to LCE and max dist LC-YC

pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat  = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
fmatout = sprintf('%sNEMO_LCcontour.mat',pthmat);

fprintf('Loading %s\n',fmatout);
load(fmatout);

dx=0.1;
[XI,YI] = meshgrid([-98:dx:-80],[14:dx:31]);
Acl = dx*dx*100;  % approximate area

% YC point:
xLC = -85;
yLC = 22;

% Calculate shortest distance from the LC 
TM = LCXY.TM;
nrec = length(TM);
nlce = length(LCE);
if f_calc>0;
		%
		% First, find closest eddy
		for irc=1:nrec
				if mod(irc,100)==0
						fprintf('Finding closest LCE, %4.1f%% done ...\n',irc/nrec*100);
				end
				X=LCXY.XY(irc).X;
				Y=LCXY.XY(irc).Y;

				DST=nan;
				for ie=1:nlce
						nee = length(LCE(ie).XY);
						if irc>nee, 
								DST(ie)=nan;
								continue;
						end

						xe=LCE(ie).XY(irc).X;
						ye=LCE(ie).XY(irc).Y;
						if isempty(xe),
								DST(ie)=nan;
								continue;
						end

		%    plot(xe,ye);
						x0=mean(xe);
						y0=mean(ye);
		%    plot(x0,y0,'*');
						IN = inpolygon(YI,XI,ye,xe);
						I = find(IN==1);
						Aeddy = length(I)*Acl;
						if Aeddy<100;
								DST(ie)=nan;
								continue;
						end

						D=distance_spheric_coord(Y,X,y0,x0);
						DST(ie)=min(D);
				end
				
				I=find(~isnan(DST));
				if isempty(I)
						LCXY.ClosestEddy(irc)=nan;
				elseif length(I)==1;
						LCXY.ClosestEddy(irc)=I;
				else
						[dmm,rI]=sort(DST);
						LCXY.ClosestEddy(irc)=rI(1);
				end 


		end

		% 
		% Calculate shortest distance to LCE
		for irc=1:nrec
				X=LCXY.XY(irc).X;
				Y=LCXY.XY(irc).Y;
				XX0=mean(X);
				Ilc=find(X<XX0); %consider only points west of middle long
				XE=X(Ilc);
				YE=Y(Ilc);


				ie=LCXY.ClosestEddy(irc);
				if isnan(ie),
						Dlce=0;
				else
						xe=LCE(ie).XY(irc).X;
						ye=LCE(ie).XY(irc).Y;
						x0=mean(xe);
						y0=mean(ye);
						I=find(xe>x0);  % consider only half closest to LC, eastern part of LCE
						xE=xe(I);
						yE=ye(I);
						nl = length(xE);
						clear D
						for jj=1:nl
								D(:,jj)=abs(XE-xE(jj))+abs(YE-yE(jj));
						end;
						[jM,iM]=find(D==min(min(D)));

						Dlce = distance_spheric_coord(yE(iM),xE(iM),YE(jM),XE(jM))*1e-3;
				end
				LCXY.Dist2LCE(irc)=Dlce;

				dLC = distance_spheric_coord(Y,X,yLC,xLC)*1e-3;
				maxD = max(dLC);
				LCXY.MaxD_LC2YC(irc)=maxD;


				fprintf('Distance: irec=%i  %8.2f km\n',irc,Dlce);
		end

		fprintf('Saving %s\n',fmatout);
		save(fmatout,'LCXY');

end


% Max dist to LC farhest end from Yuc Ch. in NEMO
MXD = LCXY.MaxD_LC2YC;
TM=LCXY.TM;
Td=TM-TM(1)+1;
Td=Td(:);
DV=datevec(TM);

cc=0;
for iyr=2011:2012
  for im=1:12
    i1=min(find(DV(:,1)==iyr & DV(:,2)==im));
    if ~isempty(i1),
      cc=cc+1;
      ttck(cc,1)=i1;
      tlbl{cc}=sprintf('%2.2i',im);
    end
  end
end


% Plot Dist LC, min D to LCE
btx ='anls_lc_nemo.m';

figure(1); clf;
axes('Position',[0.08 0.5 0.85 0.42]);
plot(MXD,'Color',[0. 0.4 0.8],'Linewidth',2.5);
ylm3=1.05*max(MXD);
ybm3=0.9*min(MXD);
nrc = length(MXD);
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[ybm3 ylm3],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

title('||LC - YC||_{supr}  NEMO LC 0.17 contours, 2011/2012');
%lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
%set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');
bottom_text(btx,'pwd',1,'Position',[0.02 0.2 0.4 0.04]);


Dlc = LCXY.Dist2LCE;
  
figure(2); clf;
axes('Position',[0.08 0.5 0.85 0.42]);
plot(Dlc,'Color',[0.9 0.4 0],'Linewidth',2.5);
ylm3=1.05*max(Dlc);
ybm3=0.9*min(Dlc);
nrc = length(Dlc);
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[ybm3 ylm3],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

title('||LC - LCE||_{infm}  NEMO LC 0.17 contours, 2011/2012');
%lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
%set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');
bottom_text(btx,'pwd',1,'Position',[0.02 0.2 0.4 0.04]);


% Scale by 1st barocl R
R1=250;   % Chelton = ~240km in subtropics/tropics, 1998
sDlc=Dlc/R1;
sMXD=MXD/R1;
figure(3); clf;
axes('Position',[0.08 0.5 0.85 0.42]);
hold on
plot(sMXD,'Color',[0. 0.4 0.8],'Linewidth',2.5);
plot(sDlc,'Color',[0.9 0.4 0],'Linewidth',2.5);
ylm3=1.05*max(sMXD);
ybm3=0;
set(gca,'tickdir','out',...
        'xlim',[1 nrc],...
        'ylim',[ybm3 ylm3],...
        'xtick',ttck,...
        'xticklabel',tlbl,...
        'xgrid','on',...
        'ygrid','on');

title('max(|LC-YC|) and min(|LC-LCE|) scaled by 1st brcl R,  NEMO LC 0.17 contours, 2011/2012');
%lgd = legend('HYCOM osse0','HYCOM osseE','HYCOM free','HYCOM AllSat');
%set(lgd,'position',[0.75 0.8 0.2 0.14]);
xlabel('Months');
bottom_text(btx,'pwd',1,'Position',[0.02 0.2 0.4 0.04]);






