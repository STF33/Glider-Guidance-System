% Interpolate hycom TSIS onto 
% NAS common grid
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

f_write=1;  % creatte netcdf fi
expt='noPIES';
ys=2011;
ye=2011;

rg=9806;  % convert pressure to depth, m
huge=1e25;

fprintf(' %s %i-%i\n',expt,ys,ye);

pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig_ssh/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
pthgrid='/Net/mars/ddmitry/hycom/hycom_TSIS/data_project/';

btx='interp_hycom2nas_grid.m';

% Read HYCOM topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

[XM,YM]=meshgrid([1:n],[1:m]);

% NAS common grid
%
fgrid=sprintf('%scommon_gom_domain.mat',pthgrid);
NAS=load(fgrid);
nlat=NAS.lat;
nlon=NAS.lon;
ynas1=min(min(nlat));
ynas2=max(max(nlat));
xnas1=min(min(nlon));
xnas2=max(max(nlon));
ZMnas=-NAS.z;


% NAS domain:
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);

XX=LON(jts1:jts2,its1:its2);
YY=LAT(jts1:jts2,its1:its2);
HHs=HH(jts1:jts2,its1:its2);
[msb,nsb]=size(HHs);

knas=length(ZMnas);
lnas=length(ZMnas);
[mnas,nnas]=size(nlon);

p_map=0;
if p_map==1
  figure(10); clf;
  contour(LON,LAT,HH,[0 0],'k');
  hold on;
  plot([xnas1 xnas1],[ynas1 ynas2],'r');
  plot([xnas2 xnas2],[ynas1 ynas2],'r');
  plot([xnas1 xnas2],[ynas1 ynas1],'r');
  plot([xnas1 xnas2],[ynas2 ynas2],'r');
  
  plot([xt1 xt1],[yt1 yt2],'g');
  plot([xt2 xt2],[yt1 yt2],'g');
  plot([xt1 xt2],[yt1 yt1],'g');
  plot([xt1 xt2],[yt2 yt2],'g');
end


for year=ys:ye
  pthd=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i_%s/',...
	       year,expt);

  pthout=sprintf('/nexsan/people/ddmitry/GOM_TSIS0.03/hindcast/%s/%4.4i/',...
		  expt,year);
  if ~exist(pthout,'dir')
    system(sprintf('mkdir -pv %s',pthout));
  end
  
  id1=1;
  id2=365;
  if (mod(year,4)==0 & id2>=365 ), id2=366; end

  for iday=id1:id2
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

    tic;

    dJ1=datenum(year,1,1);
    nday = dJ1+iday-1;
    dnmb = nday+hr/24;
    DV = datevec(dnmb);

    fprintf('Reading %s\n',fina);

  %  date_str=datestr(nday,29);
  %  disp(date_str);
    fprintf('\n iday=%i, %4.4i/%2.2i/%2.2i:%2.2ih\n',iday, DV(1:4));
    date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:4));

% Prepare Land mask
    if ~exist('Lmask','var')
      LL=HHs*0;
      LL(HHs<0)=1;
      Ioc=find(LL==1); % ocean hycom tsis
      Ild=find(LL==0); % land points, hycom
      Lmask=interp2(XX,YY,LL,nlon,nlat);
      Lmask=round(Lmask);
      HHs(Ild)=0.;
      Hnas=interp2(XX,YY,HHs,nlon,nlat);
      HHs(Ild)=100;
      Hnas(Lmask==0)=100;
      Ib=find(Hnas>=0 & Lmask==1);
      if ~isempty(Ib),
	Hnas(Ib)=-1; 
      end

  % Create Land/bottom mask indices:
     for kk=1:knas
       I=find(Hnas>=ZMnas(kk));
       ILAND(kk).I=I;
     end

    end


    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F/rg;
    F(F>1e10)=0;
    dH=squeeze(F); % note this is not U,V thickness, center grid
	       % there are thck-u & thck-v 
    ih=find(dH<1e-1);
    dHs=dH(:,jts1:jts2,its1:its2);
    
    [ZM1,ZZ1] = sub_zz_zm(fina, finb,HH,'thknss',dH);
    ZM1(isnan(ZM1))=-100;
    ZZ=ZZ1(:,jts1:jts2,its1:its2);
    ZM=ZM1(:,jts1:jts2,its1:its2);

  
% Interpolate     
% Sea surface height:
    fld = 'srfhgt';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    ssh=squeeze(F)./(1e-3*rg);  % ssh m
    ssh=ssh(jts1:jts2,its1:its2);

    if ~exist('IFL','var');
      IFL=[];
    end

    fprintf('Interpolating ssh \n');
    [ssh,IFL]=sub_fill_land(ssh,IFL);

    sshi=interp2(XX,YY,ssh,nlon,nlat);
    il=ILAND(1).I;
    sshi(il)=huge;

    
% ------------------------    
%      Temperature:  
% ------------------------    
    fld = 'temp';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    F=squeeze(F);  % 
    T=F(:,jts1:jts2,its1:its2);

%
% Interpolate onto NAS grid:    
    Tzx=sub_interp3D_nas(T,IFL,XX,YY,HHs,ZMnas,ZM,nlon,nlat,ILAND,fld);

    f_sct=0;
    if f_sct==1
      ii=400;
      ZMh=squeeze(ZM(:,:,ii));
      F1=squeeze(T(:,:,ii));
      x1=XX(1,ii);
      dd=abs(nlon(1,:)-x1);
      ii2=find(dd==min(dd),1);
      F2=squeeze(Tzx(:,:,ii2));
      F2(F2>1e20)=nan;
      Hs=Hnas(:,ii2);
      sub_check_sect(F1,F2,ZMh,ZMnas,Hs);
      bottom_text(btx,'pwd',1);
    end

% 
% ------------------------    
%      Salinity
% ------------------------    
    fld = 'salin';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    F=squeeze(F);  % 
    S=F(:,jts1:jts2,its1:its2);
    
% Interpolate onto NAS grid:    
    Szx=sub_interp3D_nas(S,IFL,XX,YY,HHs,ZMnas,ZM,nlon,nlat,ILAND,fld);
%
    f_sct=0;
    if f_sct==1
      ii=400;
      ZMh=squeeze(ZM(:,:,ii));
      F1=squeeze(S(:,:,ii));
      x1=XX(1,ii);
      dd=abs(nlon(1,:)-x1);
      ii2=find(dd==min(dd),1);
      F2=squeeze(Szx(:,:,ii2));
      F2(F2>1e20)=nan;
      Hs=Hnas(:,ii2);
      sub_check_sect(F1,F2,ZMh,ZMnas,Hs);
      bottom_text(btx,'pwd',1);
    end

% 
% ------------------------    
% Need to collocate with V into 
% center of grid cell
% and add barotropic
%     U-vel
% ------------------------    
    [F,n,m] = read_hycom(fina,finb,'u_btrop');
    F(F>1e6)=nan;
    F=squeeze(F);
    Ub=F(jts1:jts2,its1:its2);

    fld = 'u-vel';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    F=squeeze(F);  % 
    U=F(:,jts1:jts2,its1:its2);
    
% Adding barotropic U    
    for kk=1:ll
      U(kk,:,:)=squeeze(U(kk,:,:))+Ub;
    end
    
% Collocate with rho-points:
    fprintf(' Collocating U\n');
%    Uclc=U;
    U(isnan(U))=0;
    Uclc=U;
    for ii=2:nsb-1
      dh1=squeeze(dHs(:,:,ii-1));
      dh2=squeeze(dHs(:,:,ii));
      dh12=0.5*(dh1+dh2);
      
      dh2=squeeze(dHs(:,:,ii));
      dh3=squeeze(dHs(:,:,ii+1));
      dh23=0.5*(dh1+dh2);
      
      u12=squeeze(U(:,:,ii));
      u23=squeeze(U(:,:,ii+1));
      u0=(u12.*dh12+u23.*dh23)./(dh12+dh23);
      a=dh12+dh23;
      Im=find(a<1e-3);
      u0(Im)=0;
      Uclc(:,:,ii)=u0;
    end
    
      
% Interpolate onto NAS grid:    
    Uzx=sub_interp3D_nas(Uclc,IFL,XX,YY,HHs,ZMnas,ZM,nlon,nlat,ILAND,fld);
    f_sct=0;
    if f_sct==1
      ii=400;
      ZMh=squeeze(ZM(:,:,ii));
      F1=squeeze(U(:,:,ii));
      x1=XX(1,ii);
      dd=abs(nlon(1,:)-x1);
      ii2=find(dd==min(dd),1);
      F2=squeeze(Uzx(:,:,ii2));
      F2(F2>1e20)=nan;
      Hs=Hnas(:,ii2);
      sub_check_sect(F1,F2,ZMh,ZMnas,Hs);
      bottom_text(btx,'pwd',1);
    end
%
% ------------------------    
%     V-vel
% ------------------------    
    [F,n,m] = read_hycom(fina,finb,'v_btrop');
    F(F>1e6)=nan;
    F=squeeze(F);
    Vb=F(jts1:jts2,its1:its2);

    fld = 'v-vel';
    [F,nn,mm,ll] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    F=squeeze(F);  % 
    V=F(:,jts1:jts2,its1:its2);

% Adding barotropic V
    for kk=1:ll
      V(kk,:,:)=squeeze(V(kk,:,:))+Vb;
    end

% Collocate with rho-points:
    fprintf(' Collocating V\n');
    V(isnan(V))=0;
    Vclc=V;
    for jj=2:msb-1
      dh1=squeeze(dHs(:,jj-1,:));
      dh2=squeeze(dHs(:,jj,:));
      dh12=0.5*(dh1+dh2);
      
      dh2=squeeze(dHs(:,jj,:));
      dh3=squeeze(dHs(:,jj+1,:));
      dh23=0.5*(dh1+dh2);
      
      v12=squeeze(V(:,jj,:));
      v23=squeeze(V(:,jj+1,:));
      v0=(v12.*dh12+v23.*dh23)./(dh12+dh23);
      a=dh12+dh23;
      Im=find(a<1e-3);
      v0(Im)=0;
      Vclc(:,jj,:)=v0;
    end
      
% Interpolate onto NAS grid:    
    Vzx=sub_interp3D_nas(Vclc,IFL,XX,YY,HHs,ZMnas,ZM,nlon,nlat,ILAND,fld);
    
    f_sct=0;
    if f_sct==1
      ii=400;
      ZMh=squeeze(ZM(:,:,ii));
      F1=squeeze(V(:,:,ii));
      x1=XX(1,ii);
      dd=abs(nlon(1,:)-x1);
      ii2=find(dd==min(dd),1);
      F2=squeeze(Vzx(:,:,ii2));
      F2(F2>1e20)=nan;
      Hs=Hnas(:,ii2);
      sub_check_sect(F1,F2,ZMh,ZMnas,Hs);
      bottom_text(btx,'pwd',1);
    end
% Check if U,V are collocated:
%colloc=check_collocatedU(HH,U,V);
%clc=0;
%utotal = check_totalU(fina,finb,clc);    
 

% Create netcdf file and write fields:    
    if  f_write==1
      flcdf=sprintf('%shycom0.03_%s_%4.4i_%3.3i.nc',pthout,expt,year,iday);
      sub_write_netcdf(flcdf,nlon,nlat,ZMnas,sshi,Tzx,Szx,Uzx,Vzx,dnmb,huge);
    end

   fprintf('1 day processed %6.2f min\n\n',toc/60);
   
  end % day
end  % year
