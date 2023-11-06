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
s_fig=1;
%f_uv=0;   % =1 - plot u,v vectors

expt='PIES';
ys=2009;
ye=2009;

id1=254;
id2=254;
dday=2;

rg=9806;  % convert pressure to depth, m
huge=1e20;


pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig_ssh/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/%s/data_anls/',expt);
btx='plot_ssh_v2.m';

if s_fig==1
  fprintf('Fig saved is ON\n');
else
  fprintf('Fig saved is OFF\n');
end


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

cntr=0;
for year=ys:ye
  pthd=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i_%s/',year,expt);

if (mod(year,4)==0 & id2>=365 ), id2=366; end

%ff=figure('Visible','off');

for iday=id1:dday:id2
  tic;
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
  
  fld = 'temp';
  [F,nn,mm,ll] = read_hycom(fina,finb,fld,'r_layer',1);
  F(F>huge)=nan;
  T=squeeze(F);
%

  date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:4));
  stt=sprintf('HYCOM_TSIS IAS0.03, SSH %s',date_str);
  if cntr==1
    figure('Position',[1064 229 1421 1113]); clf;
  end
  fnb=1;
  c1=26;
  c2=31;
  pos=[0.12 0.1 0.82 0.82];
  xl1=-98;
  xl2=-81;
  yl1=18;
  yl2=30.;
  clb=3;
  CMP = colormap_PBYR(200,c1,c2);
  cmp = CMP.colormap;

  sub_plot_2fld(pos,T,LON,LAT,HH,c1,c2,stt,IN,xl1,xl2,yl1,yl2,clb,'cmp',cmp);
  set(gcf,'Position',[1000 481 1383 835]);
  
  bottom_text(btx,'pwd',1,'Position',[0.08 0.06 0.4 0.05]);

keyboard  
%  if s_fig==1
%    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
%    fnm=sprintf('tsis-%s_ssh_%4.4i%2.2i%2.2i',expt,DV(1:3));
%    fout=[pthfig,fnm];
%    fout=sprintf('%stsis-%s_ssh_%4.4i',pthfig,expt,cntr);
%    fprintf('Saving %s\n',fout);
%    set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
%    print('-dpng','-r150',fout);
%  end
%keyboard

  fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);
  end;  % for iday
end; %% year

%exit



