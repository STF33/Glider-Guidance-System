% Plot ssh subtacted area average SSH from TSIS analysis
% HYCOM2.3-TSIS IAS0.03 nested into GLBb0.08 GOFS3.1
%
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

expt='freerun';
ys=2019;
ye=2019;

id1=3;
id2=60;
dday=4;

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
%  pthd=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i_%s/',year,expt);
  pthd=sprintf('/nexsan/people/ddmitry/hycom/TSIS/IASx0.03/output_gofs3.1_freerun/%4.4i/',year);

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

  date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:4));
  stt=sprintf('HYCOM2.3_TSIS IAS0.03, GOFS3.1 nest, SSH %s',date_str);
  if cntr==1
    figure('Position',[1064 229 1421 1113]); clf;
  end
  fnb=1;
  c1=-0.5;
  c2=0.5;
  xl1=-98.08;
  xl2=-56.08;
  yl1=7.0025;
  yl2=31.927;

  sub_plot_ssh2(fnb,ssh,LON,LAT,HH,c1,c2,stt,IN,xl1,xl2,yl1,yl2);
%  set(gcf,'Position',[1000 481 1383 835]);
  set(gcf,'Position',[69 167 1055 633]);
 
 
  bottom_text(btx,'pwd',1,'Position',[0.08 0.06 0.4 0.05]);
  
  if s_fig==1
%    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
%    fnm=sprintf('tsis-%s_ssh_%4.4i%2.2i%2.2i',expt,DV(1:3));
%    fout=[pthfig,fnm];
    fout=sprintf('%stsis-%s_ssh_%4.4i',pthfig,expt,cntr);
    fprintf('Saving %s\n',fout);
%    set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
    print('-dpng','-r150',fout);
  end
%keyboard

  fprintf('Processed 1 rec, %6.4f min\n\n',toc/60);
  end;  % for iday
end; %% year

%exit



