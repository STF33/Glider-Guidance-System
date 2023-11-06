% Plot Yucatan Vol Flux
% for free forecast runs 3-mo
%
% calculated in
% calc_volFlux_GLBb.m - Global reanalysis
% and calc_volFlux_v3.m for HYCOM-TSIS
% Transport saved by years

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all


esim1='PIES'; % hindcast data assimilative runs with deep PIES obs
esim2='GLBu191'; % 
%esim3='noPIES'; % hindcast data assimilative runs with deep PIES obs
%esim4='freerun'; % hindcast data assimilative runs with deep PIES obs

%f_fcst=0; %=1 - plot 3-mo forecast runs

pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';

xsct_name = 'Yucatan';

ys=2009;
ye=2009;
im1=1;
im2=12;


% Get hindcast and Global reanalysis:
TbH=[];
TH=[];
TbG=[];
TG=[];
TM=[];
for iyr=ys:ye
  fmat1 = sprintf('%s%s_VolTrt_%s_%4.4i-%4.4i.mat',...
		 pthmat,esim1,xsct_name,iyr,iyr);
  fmat2 = sprintf('%s%s_VolTrt_%s_%4.4i-%4.4i.mat',...
		 pthmat,esim2,xsct_name,iyr,iyr);


  fprintf('Loading %s\n',fmat1);
  load(fmat1);
  aa=TRP(1).TrLayers;
  T1=nansum(aa,2); % layer averaged
  Tb=TRP(1).TBtp;
  aa=TRP(2).TrLayers;
  T2=nansum(aa,2); % layer averaged
  Tb=Tb+TRP(2).TBtp;
  T=T1+T2;
  Time=TRP(1).Time;
  Time=Time(:);
  Tb=Tb(:);
  T=T(:);
% Missing Jan 1, some years:
  dj1=datenum(iyr,1,1);
  i1=find(Time==dj1);
  if isempty(i1),
    Time=[dj1;Time];
    Tb=[NaN;Tb];
    T=[NaN;T];
  end
  
  TM=[TM;Time];
  TbH=[TbH;Tb];
  TH=[TH;T];
  
  fprintf('Loading %s\n',fmat2);
  load(fmat2);
  aa=TRP(1).TrLayers;
  T1=nansum(aa,2); % layer averaged
  tbg=TRP(1).TBtp;
  aa=TRP(2).TrLayers;
  T2=nansum(aa,2); % layer averaged
  tbg=tbg+TRP(2).TBtp;
  tg=T1+T2;
  
  tbg=tbg(:);
  tg=tg(:);
  TbG=[TbG;tbg];
  TG=[TG;tg];
  
end

TMh=TM;
DV=datevec(TM);
YRh=[];
for iyr=ys:ye
  I=find(DV(:,1)==iyr);
  dmm=[iyr:1/length(I):iyr+0.99999]';
  YRh=[YRh;dmm];
end





% Get free runs 3-mo segments:
% Day 1st missing - take from hindcast
rsg=0;
%if f_fcst==1
for iyr=ys:ye
  TM=[];
%  TbH=[];
  Trp=[];
%  TbG=[];
%  TG=[];
  for imo=im1:im2
    fmat1=sprintf('%sfcst_VolTrt_Yucatan_%4.4i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat1);
    load(fmat1);
    
    if ~exist(fmat1,'file'); continue; end;

    aa=TRP(1).TrLayers;
    T1=nansum(aa,2); % layer averaged
%    Tb=TRP(1).TBtp;  % barotropic transp, for checking
    aa=TRP(2).TrLayers;
    T2=nansum(aa,2); % layer averaged
%    Tb=Tb+TRP(2).TBtp;
    T=T1+T2;
    Time=TRP(1).Time;
    Time=Time(:);
%    Tb=Tb(:);
    T=T(:);
% Missing 1st day of simulation:
    dj1=datenum(iyr,imo,1);
    i1=find(Time==dj1);
    if isempty(i1),
      ihnd=find(TMh==dj1);
      tr0=TH(ihnd);
      Time=[dj1;Time];
%      Tb=[NaN;Tb];
      T=[tr0;T];
    end
  
    TM=Time;
%    TbH=Tb;
    Trp=T;

    DV=datevec(TM);
    YR=[];    
    ndays=length(TM);
    dt=TM-TM(1);
    dyr=TM(1)-datenum(DV(1),1,1)+1;
    ndyr=datenum(DV(1),12,31)-datenum(DV(1),1,1)+1;
    YR=[iyr+dyr/ndyr:1/ndyr:iyr+dyr/ndyr+(ndays-1)/ndyr]';

    rsg=rsg+1;
    TFCST(rsg).YR=YR;
    TFCST(rsg).Trnsp=Trp;
    TFCST(rsg).dv=datevec(Time(1));
    
  end
end
%end


% tick labels for months:
tmm=[datenum(ys,1,1):datenum(ye,12,31)+90];
yrm=[];
DV=datevec(tmm);
for iyr=DV(1,1):DV(end,1)
  I=find(DV(:,1)==iyr);
  dmm=[iyr:1/length(I):iyr+0.99999]';
  yrm=[yrm;dmm];
end

cc=1;
yr=DV(1,1);
mo=DV(1,2);
xtk(1,1)=1;
xlb{1}='01';
for imm=1:length(TMh);
  if (mo==12);
    yr=yr+1;
    mo=1;
  else
    mo=mo+1;
  end
  I=min(find(DV(:,1)==yr & DV(:,2)==mo));
  if ~isempty(I);
    cc=cc+1;
    xtk(cc)=yrm(I);
    xlb{cc}=sprintf('%2.2i',mo);
  end
end
  
  


figure(1); clf;
axes('Position',[0.08 0.47 0.85 0.44]);
hold on;
plot(YRh,TH*1e-6,'k','linewidth',2);
plot(YRh,TG*1e-6,'k','Color',[0.6 0.6 0.6],'linewidth',1.6);

txp{1}='hndcst';
txp{2}='GLBb';


%lg=legend('TSIS','GLBb','Location','SouthOutside');

CLR=colormap(jet(12));

for ik=1:rsg
  YR=TFCST(ik).YR;
  TF=TFCST(ik).Trnsp*1e-6;
  clr=CLR(ik,:);
  plot(YR,TF,'Color',clr);  
  yy=TFCST(ik).dv(1);
  mm=TFCST(ik).dv(2);
  txp{ik+2}=sprintf('%2.2i/%4.4i',mm,yy);
end

%set(lg,'Fontsize',12,'Position',[0.8 0.35 0.07 0.05]);
set(gca,'tickdir','out',...
	'xtick',xtk,...
	'xticklabel',xlb,...
	'ytick',[0:5:35],...
	'xlim',[ys max(YR)],...
	'ylim',[0 36],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',12);
stx{1}=sprintf('TSIS=%3.1f Sv',nanmean(TH*1e-6));
stx{2}=sprintf('GLBb=%3.1f Sv',mean(TG*1e-6));
text(YR(2),5,stx,'Fontsize',12);

lts=sprintf('Forecasts, Yucatan Transport, Sv, %i-%i',ys,ye);
title(lts);
xlabel('Months');

axes('Position',[0.08 0.1 0.8 0.3]);
hold on;
cff=1;
y0=1;
ntxt=length(txp);
for ik=1:rsg
  if mod(ik,10)==0;
    cff=cff+1;
    y0=1;
  end
  x0=(cff-1)*0.2;
  y0=y0-1/10;
  
  if ik==1
    clr=[0 0 0];
  elseif ik==2
    clr=[0.6 0.6 0.6];
  else
    clr=CLR(ik-2,:);
  end
  
  plot([x0 x0+0.02],[y0 y0],'-','Color',clr,'linewidth',1.8);
  text(x0+0.03,y0,txp{ik},'Fontsize',12);
end
set(gca,'Visible','off',...
	'xlim',[-0.05 1],...
	'ylim',[-0.05 1]);


btx='plot_fcstVolFlux.m';
bottom_text(btx,'pwd',1);


