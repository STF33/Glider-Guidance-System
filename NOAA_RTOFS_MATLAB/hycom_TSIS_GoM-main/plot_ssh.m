% Plot ssh subtacted area average SSH from TSIS analysis
% submit:  matlab -nodesktop -nosplash < plot_ssh5.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

% ----------------
% Flags
% ---------------
sfig=logical(0);
f_uv=0;   % =1 - plot u,v vectors
f_cont = 0; %=1 - consecutive numbering through all years 

ys=2009;
ye=2009;
dday = 1; % day stepping for plotting
ds=1;
%ms=305;
%dJ1 = datenum(ys,1,1);
%nday0=dJ1+ds-1;

de=2;
id1=365;
id2=365;

rg=9806;  % convert pressure to depth, m
huge=1e20;


pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/%s/data_anls/',expt);
btx='plot_ssh.m';

% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);


% GoM region:
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
clear XM YM

c1=-0.5;
c2=0.5;

nc=50;
clB=flipud(colormap_blue(nc));
clR=colormap_red(nc);
clB(nc-3,:)=[1 1 1];
clB(nc-2,:)=[1 1 1];
clB(nc-1,:)=[1 1 1];
clB(nc,:)  =[1 1 1];
clR(4,:)   =[1 1 1];
clR(3,:)   =[1 1 1];
clR(2,:)   =[1 1 1];
clR(1,:)   =[1 1 1];
cmp=[clB;clR];
cmp=smooth_colormap(cmp,5,3);
nint=length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % fake colormap to color gray nan shallow regions
%


cntr=0;
for year=ys:ye
  pthd=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i/',year);

if (mod(year,4)==0 & id2>=365 ), id2=366; end

%ff=figure('Visible','off');

for iday=id1:dday:id2
  sday=sprintf('%3.3i',iday);
  hr=0;
  fina=sprintf('%sarchv.%4.4i_%3.3i_%2.2i.a',pthd,year,iday,hr);
  finb=sprintf('%sarchv.%4.4i_%3.3i_%2.2i.b',pthd,year,iday,hr);
  fin=fina;

  ie = exist(fin,'file');
  
  if ~ie
    fprintf('Missing: %s\n',fin);
    continue;
  end
  

  cntr=cntr+1;
  dJ1=datenum(year,1,1);
  nday = dJ1+iday-1;
  dnmb = nday+hr/24;
  DV = datevec(dnmb);
  
  fprintf('Reading %s\n',fina);
  
%  date_str=datestr(nday,29);
%  disp(date_str);
  fprintf('\n iday=%i, %4.4i/%2.2i/%2.2i:%2.2ih\n',iday, DV(1:4));
  date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:4));
  
  fld = 'srfhgt';
  [F,nn,mm,ll] = read_hycom(fina,finb,fld);
  F(F>huge)=nan;
  ssh=squeeze(F)./(1e-3*rg);  % ssh m
%
% Subtract anomaly:
  dmm=ssh;
  dmm(IN==0)=nan;
  dmm(HH>-200)=nan;
  sshM=nanmean(nanmean(dmm));
  ssh=ssh-sshM;

  clf
  pcolor(LON,LAT,ssh); shading flat;
  caxis([c1 c2]);
  colormap(cmp);
  hold on;
  contour(LON,LAT,HH,[-200 -200],'Color',[0.8 0.8 0.8]);
  contour(LON,LAT,HH,[-5000:1000:-1000],'Color',[0.7 0.7 0.7]);
%  title(['SSH, ',regn,' - ',expt,' Year: ',int2str(year),'; D= ',sday]);
  stt=sprintf('TSIS-IAS0.03, SSH %s',date_str);
  title(stt);
  
%keyboard
  axis('equal');
  set(gca,'xlim',[-98 -56.08],'ylim',[7.0025 31.926]);
  set(gca,'tickdir','out','Color',[0 0 0]);
%  set(gcf,'Color',[1 1 1],'InvertHardCopy','off');


  if f_uv
    fld = 'u-vel.';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    u=squeeze(F(1,:,:));  % 1st layer U
    fld = 'v-vel.';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    v=squeeze(F(1,:,:));  % 1st layer V
  end   % if f_uv

% ==================
% Plot ssh
% =================
  contour(LON,LAT,ssh,[0:0.1:1],'w','linewidth',1.);
  contour(LON,LAT,ssh,[-1:0.1:-0.01],'w--','linewidth',1.);

  chb = colorbar;
  set(chb,'Fontsize',16);

  bottom_text(btx,'pwd',1,'Position',[0.08 0.12 0.4 0.05]);
  drawnow
  
  if sfig==1
%    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
    fnm=sprintf('ssh%s-%4.4i',expt,cntr);
    fout=[pthout,fnm];
    fprintf('Saving %s\n',fout);
    set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
    print('-dpng','-r150',fout);
  end
%keyboard

  end;  % for iday
end; %% year

%exit



