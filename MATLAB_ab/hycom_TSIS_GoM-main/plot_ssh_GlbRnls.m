% SSH demeaned
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig=1;
rg=9806;  % convert pressure to depth, m
huge=1e20;

ys=2009;
ye=2009;
id1=136;
id2=365;
dday=15;


pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_ssh/';
ftopo = sprintf('%sGLBb0.08_T11_TSIS_IAS.mat',pthmat);
pthglb = '/nexsan/hycom/GLBv0.08_53X/data/';
pthglbtopo = '/nexsan/hycom/GLBv0.08_53X/topo/';
btx='plot_ssh_GlbRnls.m';

if s_fig==1
  fprintf('Fig saved is ON\n');
else
  fprintf('Fig saved is OFF\n');
end


%HG=load(ftopo);
%HHg=HG.Hglb;
%LONg=HG.LONglb;
%LATg=HG.LATglb;
%I1=HG.Indx(1);
%I2=HG.Indx(2);
%J1=HG.Jndx(1);
%J2=HG.Jndx(2);
%di=HG.di;
%dj=HG.dj;

% Find location of the TSIS grid:
% TSIS IAS grid:
% Read topography:
ftopots = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopots,'mdepth'));
LAT = nc_varget(ftopots,'mplat');
LON = nc_varget(ftopots,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;


for yr=ys:ye
  for iday=id1:dday:id2
    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    mo=DV(2);
    dd=DV(3);
%expt = sub_exptGLBb(dnmb);
    expt=sub_exptGLBb53X(dnmb);
    sexpt='53.X'; % GOFS3.1 Reanalysis, 1994-2015/12/31

    hr=0;
    fnm=sprintf('hycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',expt,yr,mo,dd,hr);
    fin=sprintf('%s%4.4i/%s',pthglb,yr,fnm);
    fprintf('Reading %s\n',fin);

% Find TSIS IAS grid in GLBb grid:
    if ~exist('HHg','var');
      lon=squeeze(nc_varget(fin,'lon'));
      lat=squeeze(nc_varget(fin,'lat'));
      [LONg,LATg]=meshgrid(lon,lat);
      D=distance_spheric_coord(LAT(1,1),LON(1,1),LATg,LONg);
      [j1,i1]=find(D==min(min(D)));
      fprintf('Closest point in GLBb is %8.4f m away from TSIS(1,1)\n',min(min(D)));
      D=distance_spheric_coord(LAT(1,end),LON(1,end),LATg,LONg);
      [j2,i2]=find(D==min(min(D)));
      fprintf('Closest point in GLBb is %8.4f m away from TSIS(1,end)\n',min(min(D)));
      D=distance_spheric_coord(LAT(end,end),LON(end,end),LATg,LONg);
      [j3,i3]=find(D==min(min(D)));
      fprintf('Closest point in GLBb is %8.4f m away from TSIS(end,end)\n',min(min(D)));
      D=distance_spheric_coord(LAT(end,1),LON(end,1),LATg,LONg);
      [j4,i4]=find(D==min(min(D)));
      fprintf('Closest point in GLBb is %8.4f m away from TSIS(end,1)\n',min(min(D)));

      di=i2-i1+1;
      dj=j3-j1+1;
    %
    % Note topo grid has different lon starting index  
      tglb=sprintf('%sdepth_GLBv0.08_11.nc',pthglbtopo);
      lng=nc_varget(tglb,'Longitude');
      ip=find(lng>180);
      lng(ip)=lng(ip)-360;
      dll=abs(lng-LON(1,1));
      it1=find(dll==min(dll),1);
      HHg=squeeze(nc_varget(tglb,'bathymetry',[0,j1-1,it1-1],[1,dj,di]));
      dmm=nc_varget(tglb,'Longitude',[i1-1],[di]);
      dmm=nc_varget(tglb,'Longitude');
      I=find(isnan(HHg));
      HHg=HHg*-1;
      HHg(I)=100;

      [mm,nn]=size(HHg);
      LONg=LONg(j1:j3,i1:i2);
      LATg=LATg(j1:j3,i1:i2);

    end

    % GoM region in the subset array
    GOM=[3   201
	 4   224
	 5   260
	37   289
       122   300
       187   296
       206   263
       220   229
       215   199
       194   198
       173   189
       139   179
       103   157
	82   130
	29   133
	 7   166
	 3   189];
    [XM,YM]=meshgrid([1:nn],[1:mm]);
    IN = inpolygon(XM,YM,GOM(:,1),GOM(:,2));

    % Mask off Pacific:
    Pf=[1   124
	40   125
	94   109
       164    59
       197    24
       214    18
       231    28
       239    28
       248    24
       269     1
       1 1];
    INp=inpolygon(XM,YM,Pf(:,1),Pf(:,2));

    ssh=squeeze(nc_varget(fin,'surf_el',[0,j1-1,i1-1],[1,dj,di]));
    % Subtract anomaly:
    dmm=ssh;
    dmm(IN==0)=nan;
    dmm(HHg>-200)=nan;
    sshM=nanmean(nanmean(dmm));
    ssh=ssh-sshM;
    ssh(INp==1)=nan;
    HHg(INp==1)=nan;

    date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:3),hr);
    stt=sprintf('GLBv0.08_%3.3i, SSH %s',expt,date_str);
    fnb=1;
    c1=-0.5;
    c2=0.5;
    sub_plot_ssh(fnb,ssh,LONg,LATg,HHg,c1,c2,stt,IN);
    set(gcf,'Position',[1000 481 1383 835]);

    bottom_text(btx,'pwd',1,'Position',[0.08 0.06 0.4 0.05]);

    drawnow

    if s_fig==1
      fnm=sprintf('GLBv008-%3.3i_ssh_%4.4i%2.2i%2.2i',expt,DV(1:3));
      fout=[pthfig,fnm];
      fprintf('Saving %s\n',fout);
      set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
      print('-dpng','-r150',fout);
    end
  end  % dday
end    % year
