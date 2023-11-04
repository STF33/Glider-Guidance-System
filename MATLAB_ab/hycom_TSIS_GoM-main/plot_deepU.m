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
s_fig = 1;

ys=2009;
ye=2010;
dday = 1; % day stepping for plotting
ds=1;
%ms=305;
%dJ1 = datenum(ys,1,1);
%nday0=dJ1+ds-1;

de=2;
id1=2;
id2=365;

rg=9806;  % convert pressure to depth, m
huge=1e20;

pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/frames_deepU/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/GoM/GOMl0.04/%s/data_anls/',expt);
btx = 'plot_deepU.m';

% Read the topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
m=mm;
n=nn;
HH(isnan(HH))=100;

[DX,DY]=sub_dx_dy(LON,LAT);

Lmsk = HH*0;
Lmsk(HH<0) = 1;


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

IDp = find(IN==1 & HH<-1000);
IB  = find(IN==0 | HH>-1000);

HH(IB)=nan;

cmpL = [0 0 0; .4 .4 .4];

c1=0.;
c2=0.3;
CMP = create_colormap_WBYR(200,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint= length(cnt);
%



ff=figure('Visible','off');
cntr=0;
for year=ys:ye
  pthd=sprintf('/Net/gleam/dmitry/hycom/TSIS/IASx0.03/output/%i/',year);

  id2=365;
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
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    ssh=squeeze(F)./(1e-3*rg);  % ssh m
  %
  % Subtract anomaly:
    dmm=ssh;
    dmm(IN==0)=nan;
    dmm(HH>-200)=nan;
    sshM=nanmean(nanmean(dmm));
    ssh=ssh-sshM;
    
    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge) = nan;
    dH = F./rg;

% Note that in the raw binary, U, V are on u,v grid
% h is on p-grid
% !!! NOTE: in the instanteneous output files archv u-vel, v-vel 
% are BAROCLINIC (not total) components
% need to add u barotropic !!!
    fld = 'v_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    vb=squeeze(F);

  % Baroclinic V vel:  
    fld='v-vel.';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    vbcl = F;

  % Get U vel 
    fld = 'u_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    ub = squeeze(F);

  % Baroclinic U vel:  
    fld='u-vel.';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    ubcl=F;

%      clear U V

    dH(dH<0.1)=nan;
    
    if ~exist('iz0','var')
      zz0=-2500;
      ZZ(1,:,:)=HH*0;
      for kk=1:l
	ZZ(kk+1,:,:)=ZZ(kk,:,:)-dH(kk,:,:);
	zmm = squeeze(ZZ(kk+1,:,:));
	zmm(HH>zz0)=nan;
	zmm(isnan(HH))=nan;
	zmn = nanmean(nanmean(zmm));
	if zmn<=zz0
	  iz0=kk;
          Zlvl = zmn;
	  break;
	end
      end
    end
    
      

% Near-bottom speed:
    clear Ub Vb
    Ub = ub+squeeze(ubcl(iz0,:,:));
    Vb = vb+squeeze(vbcl(iz0,:,:));

    Sb = sqrt(Ub.^2+Vb.^2);
    Sb(IB) = nan; % OB
    Sb(HH>Zlvl | isnan(HH))=nan; 
    Ub(HH>Zlvl | isnan(HH))=nan;
    Vb(HH>Zlvl | isnan(HH))=nan;

% ==========================
%   Plot deep flow:
% ==========================
    clf;
%      contour(lon,lat,HH,[-100 -100],'Color',[0.7 0.7 0.7],'linewidth',1);
    pcolor(LON,LAT,Lmsk); shading flat;
    colormap(cmpL);
    freezeColors;

    hold on;

    pcolor(LON,LAT,Sb); shading flat;
    caxis([c1 c2]);
    colormap(cmp);

    contour(LON,LAT,HH,[-100 -100],'Color',[0.8 0.8 0.8]);
    contour(LON,LAT,HH,[-5000:1000:-1000],'Color',[0.6 0.6 0.6]);

    axis('equal');
    set(gca,'tickdir','out',...
            'xlim',[-97.5 -83.8],...
	    'ylim',[19 30]);

% Plot vectors:
    fprintf('Drawing vectors ...\n');
    Ib=find(~isnan(Ub));
    dltx=25;
    cf=0.25;
    beta=30;
    v_col=[0 0. 0];
%    v_col=[1 0 0];
    lwd=1.;
    scl=0.25;
    for iib=1:dltx:length(Ib)
      I0=Ib(iib);
      [j,i] = ind2sub(size(HH),I0);
      if ~isnan(Sb(j,i)) & Sb(j,i)>0.02
        u=Ub(j,i);  % m/s
	v=Vb(j,i);
	ss=sqrt(u.^2+v.^2);
	u=u/ss;
	v=v/ss;
	x1=LON(j,i);
	x2=x1+u*scl;
	y1=LAT(j,i);
	y2=y1+v*scl;
	draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
      end;
    end;  % indices loop

% Plot ssh
  contour(LON,LAT,ssh,[0:0.1:0.9],'Color',[1 1 0.7]);
  contour(LON,LAT,ssh,[-0.5:0.05:-0.01],'k--','Color',[1 1 0.7]);

%  title(['SSH, ',regn,' - ',expt,' Year: ',int2str(year),'; D= ',sday]);
  stt=sprintf('SSH(+0.1/-0.05), TSIS-IAS0.03, U, %im, %s',round(abs(Zlvl)),date_str);
  title(stt,'Fontsize',11);
  

  chb = colorbar;
  set(chb,'Position',[0.91 0.12 0.02 0.78],...
	  'TickLength',0.032);

  bottom_text(btx,'pwd',1,'Fontsize',4);
%  drawnow
  
  if s_fig==1
%    fnm=sprintf('goml%s-SSH%4.4i_%s',expt,year,sday);
    fnm=sprintf('tsis_ias003_deepU-%4.4i',cntr);
    fout=[pthfig,fnm];
    fprintf('Saving %s\n',fout);
    set(gcf,'InvertHardcopy','off','Color',[1 1 1]);
    print('-dpng','-r250',fout);
  end
  
  fprintf('1 step processed: %6.3f min\n\n',toc/60);
%keyboard

  end;  % for iday
end; %% year

%exit

