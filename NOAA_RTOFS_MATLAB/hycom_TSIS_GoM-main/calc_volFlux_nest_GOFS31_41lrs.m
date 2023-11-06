% New HYCOM-TSIS:
% HYCOM2.3, with 41 vertical layers
% nested within GOFS3.1 GLBb0.08 expt 73.7
% 
%
% Calculate Mass flux budget in a closed volume
% - saves transport by layers
%
% All output - averaged in time
% NOTE: u,v, h are not collocated in the nest files, 
% they are in U-, V-, P-grid
%  HYCOM grid:
%
%          ------------   | V(i,j+1) ------------------
%          |                             |
%          |                             |
%        --- U(i,j)      *P(i,j)        ---  U(i+1,j)
%
%          |                             |
%          |                             |
%          ------------   | V(i,j) ------------------
%
%
% Layer thickness (dP) changes in time
% So, mean U is calculated:
% mean(U(level=k))=sum(U(t=1,k)*dP(t=1,k)+U(t=2,k)*dP(t=2,k)+...)/(sum(dP(t=1,k)+dP(t=2,k)+...)
% mean layer thicknesses is saved as well
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
%addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

% Note nest files are created for 4-day intervals
% use only for archm files - U&V collocated and U is total
HR=0;
esim='gofs31_nest';
if strncmp(esim,'fcst',4),
 ffcst=1;
else
  ffcst=0;
end


% Years to calculate:
ys = 2019;
ye = ys;  % do by 1 year at a time
im1 = 1; % only for forecasts - start month of fcst, 
         % otherwise ignored
im2 = 12;


pltx=0;
sfig=0;
f_mat = 2;  % = 0 - do not save,
            % = 1 - save matlab file with transport, etc.;
            % = 2 - load saved mat file and continue from where the last day is 
            %      in the mat file till the last day = dy2
	    % = 3 - load and plot what there is in the mat file
	    %       1 year only, for all years - use plot_YucTransport_v2.m

xsct_name = 'Yucatan';
%xsct_name = 'Yuc2';
%xsct_name = 'FlStr';
%xsct_name = 'GoMOB';
%xsct_name = 'DomainOB'; % OB of the whole domain


% -------------------------
% Directories:
% -------------------------
pthmat = '/Net/kronos/ddmitry/hycom/TSIS/datamat3/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
btx = 'calc_volFlux_nest_GOFS31_41lrs.m';


rg=9806;  % convert pressure to depth, m
huge=1e10;

mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(ys,4)==0; mday(2)=29; end;

YRPLT=[];
if ffcst==0
  cc=0;
  for iyr=ys:ye
    jd1=datenum(iyr,1,1);
    id1=1;
    id2=365;
    ic=mod(iyr,4);
    if ic==0 & id2==365, id2=366; end;
    if iyr>ys
      id1=1;
    end
    for iday=id1:id2
      dnmb=jd1+iday-1;
      DV = datevec(dnmb);
			imo = DV(2);

			if imo<im1 | imo>im2; continue; end

      cc=cc+1;
      YRPLT(cc,1)=iyr;
      YRPLT(cc,2)=iday;
      YRPLT(cc,3)=dnmb;
    end
  end
else
% Dates for forecast:
  jd1=datenum(ys,1,1);
  id11=datenum(ys,imo,1)-jd1+1; % start Yr. day
  dnmb1 = datenum(ys,imo,1);
  if ys==2009 & imo==6
    dnmb1 = datenum(ys,imo,16);
  end
  dnmb2=datenum(ys,imo,1)+100;
  dv=datevec(dnmb2);
  ye=dv(1);
  ime=dv(2);
  dnmb2=datenum(ye,ime,1)-1;
  jd2=datenum(ye,1,1);
  dv2=datevec(dnmb2);
  ye=dv2(1);
  id22=dnmb2-jd2+1;
  ndays = dnmb2-dnmb1+1;
 
  cc=0;
  dnmb=dnmb1-1;
  for idd=1:ndays
    dnmb=dnmb+1;
    cc=cc+1;
    DV=datevec(dnmb);
    jd1=datenum(DV(1),1,1);
    yday=dnmb-jd1+1;
    YRPLT(cc,1)=DV(1);
    YRPLT(cc,2)=yday;
    YRPLT(cc,3)=dnmb;
  end
  
end

fprintf('Data extraction: %4.4i/%2.2i - %4.4i/%2.2i\n',...
	YRPLT(1,1),YRPLT(1,2),YRPLT(end,1),YRPLT(end,2));
fprintf('Data extraction: %s - %s\n',...
	datestr(YRPLT(1,3)),datestr(YRPLT(end,3)));
  


if ~exist(pthmat,'dir')
  create_directories(pthmat);
end;
if ~exist(pthfig,'dir')
  sst=sprintf('mkdir -pv %s',pthfig);
  system(sst);
end;

% ------------------------------------------------------

UFLX = struct;
UFLX.esim = esim;
TRP = struct;

fprintf('Caclulating Volume Transp, Y=%i, esim=%s \n',...
	ys,esim);
% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
% Read topography:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
alat = nc_varget(ftopo,'mplat');
elon = nc_varget(ftopo,'mplon');
[mm,nn]=size(HH);
JD=mm;
ID=nn;
m=mm;
n=nn;
HH(isnan(HH))=100;
IJDM = ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

[DX,DY]=sub_dx_dy(elon,alat);

% Specify start-end points for the sections (vertices)
%
%  NW 
%   *-------*  NE
%   |       |
%   |       |
%   |       |
%  *-------* SE
%  SW

switch (xsct_name),
  case ('Yucatan');    % calculate flux in the Yucatan Channel
% Indices of the 2 sides: Yucatan Channel transport:
% Sum over Southern and Eastern sides:
    SWE = [-86.89, 21.32; -84.24, 21.32]; % Weast- East section
    SSN = [SWE(2,1), SWE(2,2); SWE(2,1), 22.0];
    splt = [1;2];   % sides where the transport is plotted, 
                    % 1 - Southern, 
                    % 2 - eastern, 3 - NOrthern, 4 - Western
    yl1=24;
    yl2=35;
    dym=2;

  case ('GoMOB');  % calculate flux along the OB 
% For the OB, sum over Sout, East, and North boundaries
% to get net flux into the GoM domain
    SWE = [-88.12, 18.12967; -76.44, 18.12967]; % Weast- East section
%    SSN = [SWE(2,1), SWE(2,2); SWE(2,1), 31.586];
    SSN = [SWE(2,1), SWE(2,2); SWE(2,1), 20.13246];
    splt = [1;2;3];  % sides where the transport is plotted, 1 - Southern, 2 - eastern, 3 - NOrthern, 4 - Western
    yl1=-1.5;
    yl2=1.5;
    dym=0.25;
  case ('Yuc2'); % calculate flux along the OB (inside the domain 1 pnt from relax. zone)
% For the OB, sum over South, East, 
% and North boundaries to get net flux into the GoM domain
    SWE = [-88.2, 18.4710; -76.84, 18.470]; % Weast- East section
    SSN = [SWE(2,1), SWE(2,2); SWE(2,1), 31.586];
    splt = [1;2;3];   % sides where the transport is plotted, 1 - Southern, 2 - eastern, 3 - NOrthern, 4 - Western
    yl1=-1.5;
    yl2=1.5;
    dym=0.25;
 case('DomainOB');
% OBs:
  iSW = 1296;
  jSW = 12;  % step few grid points  from the bndry in the domain
  iSE = n-11;
  jSE = jSW;
  iNE = iSE;
  jNE = m-11;
  iNW = iSW;
  jNW = jNE;
 case('DomainOBnp'); % step in the domain n points
% OBs:
  np=3;
  iSW = 252;
  jSW = 1+np;  % step few grid points  from the bndry in the domain
  iSE = n-np;
  jSE = jSW;
  iNE = iSE;
  jNE = m-np;
  iNW = iSW;
  jNW = jNE;
 case('DomainOB5'); % step in the domain n points
% OBs:
  np=5;
  iSW = 1296;
  jSW = np;  % step few grid points  from the bndry in the domain
  iSE = n-np;
  jSE = jSW;
  iNE = iSE;
  jNE = m-np;
  iNW = iSW;
  jNW = jNE;
    
end;


% Find indices:
if ~exist('iSW','var');
  clear IY
  for k=1:2
    if sign(elon(1))~=sign(SWE(k,1)); error('ERR: Need to convert longitudes'); end;
    DST = distance_spheric_coord(alat,elon,SWE(k,2),SWE(k,1));
    [j,i]=find(DST==min(min(DST)));
    IY(k,:)=[i,j];
  end
  DST = distance_spheric_coord(alat,elon,SSN(2,2),SSN(2,1));
  [j,i]=find(DST==min(min(DST)));
  IY(3,:)=[i,j];

  iSW = IY(1,1);
  jSW = IY(1,2);  % step few grid points  from the bndry in the domain
  iSE = IY(2,1);
  jSE = IY(2,2);
  iNE = IY(3,1);
  jNE = IY(3,2);
  iNW = iSW;
  jNW = jNE;
end


NB  = 4;     % # of boundaries that flux is calculated through, 
             % if NB<4, the last OBs (>NB) are not counted

fmat = sprintf('%s%s_VolTrt_%s_%4.4i-%4.4i.mat',...
     pthmat,esim,xsct_name,ys,ye);

% ------------------------------------------------------
ip1=1;
if f_mat<=2,   % extract data
  BB.esim= esim;
  BB.iSW = iSW;
  BB.jSW = jSW;
  BB.iSE = iSE;
  BB.jSE = jSE;
  BB.iNE = iNE;
  BB.jNE = jNE;
  BB.iNW = iNW;
  BB.jNW = jNW;
  BB.HH  = HH;
  BB.NB  = 4;
  BB.fplt= 0;
  [UFLX,TRP] = sub_UFLX_TRP(BB);
  
  
  cc=0;
  yrOLD=YRPLT(1,1);
  if f_mat==2
    if ~exist(fmat,'file'),
      error('MAT file does not exist, change option f_mat to 1');
    end
    load(fmat);
    TM=TRP(1).Time;
    dd=TM(end);
    dv=datevec(dd);
    cc=length(TM);  % do over the last record
  %    dy1=dy1+1;  % next day
    fprintf('\n Appending to the existing transport\n');
    fprintf('\n Last saved record %s\n\n', ...
	    datestr(dd));
    cntr=UFLX(1).counter;
    if (cntr ~= cc), 
      error('UFLX and TRP arrays counters do not agree');
    end
%    ip1=cc+1;
    TMy=datenum(YRPLT(1,1),1,1) + YRPLT(:,2)-1;
    ip1=find(TMy==dd)+1;
    if isempty(ip1),
      error('Couldnot find date to start ...');
    end
    
    yrOLD=dv(1);
%    keyboard
  end;

  nrc=size(YRPLT,1);
  for ip=ip1:nrc
    iyr=YRPLT(ip,1);
    iday=YRPLT(ip,2);
%    icyc=YRPLT(ip,3);
    dJ1=datenum(iyr,1,1);
    dnmb=dJ1+iday-1;
    DV=datevec(dnmb);
    dd=dnmb;
%keyboard
    fprintf('Reading: %s; esim=%s\n',datestr(dnmb),esim);

  % Update at new year: 
  % Save mat file
    if iyr~=yrOLD
      yrOLD=iyr;

      if cc>0 & f_mat>0
				fprintf('Saving mat file: %s\n',fmat);
				save(fmat,'TRP','UFLX');
      end

      cc=0;
      
      BB.fplt= 0;
      [UFLX,TRP] = sub_UFLX_TRP(BB);
      
    end

		pthi=sprintf('/Net/kronos/ddmitry/hycom/TSIS/nest_files_gofs3.1/%4.4i/',iyr);
		fina = sprintf('%sarchm.GLBb008ias003topo_%i_%3.3i.a',pthi,iyr,iday);
		finb = sprintf('%sarchm.GLBb008ias003topo_%i_%3.3i.b',pthi,iyr,iday);

    if ~exist(fina,'file')
      fprintf('Missing %s\n',fina);
      continue
    end


    cc=cc+1;
    tic;
    
  %  fld='montg1';
  %  [F,n,m,l] = read_hycom(fina,finb,fld);
  %  F(F>huge)=nan;
  %  Mt=squeeze(F);

  % sea surface height (ssh) in the archv. binary is
  % "srfhgt" to get it in m: ssh = srfhgt/(thref*9806), thref=1e-3
    fld='srfhgt';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    SSH=squeeze(F)./(1e-3*rg);  % ssh m
    for kb=1:NB
      OB = UFLX(kb).OB;
      js1= UFLX(kb).IJ(1,2);
      js2= UFLX(kb).IJ(end,2);
      is1= UFLX(kb).IJ(1,1);
      is2= UFLX(kb).IJ(end,1);

      clear ssh
      if strcmp(OB,'S') | strcmp(OB,'N')
				j0=UFLX(kb).IJ(1,2);
				II=UFLX(kb).IJ(:,1);

				for ip=1:length(II)
					i0=II(ip);
					ssh(ip)=0.5*(SSH(j0,i0)+SSH(j0-1,i0));
				end;

      elseif strcmp(OB,'E') | strcmp(OB,'W');
				i0=UFLX(kb).IJ(1,1);
				II=UFLX(kb).IJ(:,2);
				for ip=1:length(II)
					j0=II(ip);
					ssh(ip)=0.5*(SSH(j0,i0)+SSH(j0,i0-1));
				end

      end
      UFLX(kb).SSH=ssh;
    end;

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    THK = F./rg;

    for kb=1:NB
      OB = UFLX(kb).OB;
      js1= UFLX(kb).IJ(1,2);
      js2= UFLX(kb).IJ(end,2);
      is1= UFLX(kb).IJ(1,1);
      is2= UFLX(kb).IJ(end,1);

  % For V normal to the bndry:
  % Layer thickness is collocated with V point
  % from j and j-1 grid cells
      clear h3
      if strcmp(OB,'S') | strcmp(OB,'N')
				j0=UFLX(kb).IJ(1,2);
				II=UFLX(kb).IJ(:,1);
        for ip=1:length(II)
          i0=II(ip);
          h3(:,ip) = THK(:,j0,i0);
        end;

      end
  % For U normal to the bndry:
  % Layer thickness is collocated with U point
  % from i and i-1 grid cells
      if strcmp(OB,'E') | strcmp(OB,'W');
				i0=UFLX(kb).IJ(1,1);
				II=UFLX(kb).IJ(:,2);
				for ip=1:length(II)
					j0=II(ip);
					h3(:,ip)=THK(:,j0,i0);
				end
      end   % if strcmp
      UFLX(kb).thk=squeeze(h3);  % layer thickness at the left & right edges of the cell normal to V
      if ~isfield(UFLX,'sumthk') | cc==1;  UFLX(kb).sumthk=0; end;
      UFLX(kb).sumthk= UFLX(kb).sumthk+UFLX(kb).thk;  % sum thicknesses
    end  % for kb



  %  =================
  % Create array of interface depths from layer thicknesses:
  %  =================
    clear ZZ ZM
    for ks=1:NB
      Dsec=UFLX(ks).thk;
      [a1,a2]=size(Dsec);
      ZZ=zeros(a1,a2);
      Dsec(isnan(Dsec))=0;
      for ii=1:a2
				for k=1:a1
					dz=abs(Dsec(k,ii));
					ZZ(k+1,ii)=ZZ(k,ii)-dz;
					if dz<1e-3, ZZ(k+1,ii)=nan; end;
				end
      end
      UFLX(ks).ZZ=ZZ;
    end

  % Note that in the raw binary, U, V are on u,v grid
  % h is on p-grid
  % !!! NOTE: in the archive mean files archm u-vel, v-vel 
  % are TOTAL components
  % NO need to add u barotropic !!!
  %
  %keyboard
    fld = 'v_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    v_btrop=squeeze(F);

  % Baroclinic V vel:  
    fld='v-vel.';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    v_vel = F;

  % Get U vel 
    fld = 'u_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    u_btrop = squeeze(F);

  % Baroclinic U vel:  
    fld='u-vel.';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>huge)=nan;
    u_vel=F;

% Check if U& V are collocated
% and if these are anomalies or full U:
%    clc=check_collocatedU(HH,u_vel,v_vel);
%    utotal=check_totalU(fina,finb,clc);

%keyboard
    for kb=1:NB
      OB = UFLX(kb).OB;
      js1=UFLX(kb).IJ(1,2);
      js2=UFLX(kb).IJ(end,2);
      is1=UFLX(kb).IJ(1,1);
      is2=UFLX(kb).IJ(end,1);

      if strcmp(OB,'S') | strcmp(OB,'N')
				bt=squeeze(v_btrop(js1:js2,is1:is2));
				bc=squeeze(v_vel(:,js1:js2,is1:is2));
      elseif strcmp(OB,'E') | strcmp(OB,'W')
				bt=squeeze(u_btrop(js1:js2,is1:is2));
				bc=squeeze(u_vel(:,js1:js2,is1:is2));
      end

      if size(bt,1)>1, bt=bt'; end;
      for k=1:l
				UFLX(kb).V(k,:)=bc(k,:); % total baroclinic 
      end
      UFLX(kb).UBRT = bt;  % keep barotropic U 
    end

  % For plotting purposes add extra layer to show the bottom:
    for kb=1:NB
      V=UFLX(kb).V;  % keep 
      [a1,a2]=size(V);
      for ii=1:a2
        V(a1+1,ii)=V(a1,ii);
      end;

      if ~isfield(UFLX,'Vav') | isempty(UFLX(kb).Vav)
    	  UFLX(kb).Vav=V*0;
      end
    end


  % Calculate flux:
  % Create Depth array of interface depths:
    Ttot=0;
    for kb=1:NB
      js1 = UFLX(kb).IJ(1,2);
      js2 = UFLX(kb).IJ(end,2);
      is1 = UFLX(kb).IJ(1,1);
      is2 = UFLX(kb).IJ(end,1);

      ni1 = is2-is1+1;
      nj1 = js2-js1+1;

      Hs  = UFLX(kb).Hs;
      thk = UFLX(kb).thk;
      V   = UFLX(kb).V;
      SSH = UFLX(kb).SSH;
      Vav = UFLX(kb).Vav;
      Ubt = UFLX(kb).UBRT;
      sthk= UFLX(kb).sumthk;
  %    SSH=SSH*0;  % ignore SSH

  % Note that positive influx is "IN" the volume
  % so at the Northern & Eastern bndry, posit. flux is opposite to +V/U
      ob=UFLX(kb).OB;
      if strcmp(ob,'N')>0 | strcmp(ob,'E')>0
      	sgn=-1;
      else
      	sgn=1;
      end

      VF=V*nan;
      D1=Hs*0;
      ni=max([ni1,nj1]);

  %keyboard
  % ===========================
  % FLUX CALCULATION:
  % ===========================
     clear AdZ
      for ii=1:ni
				i0=UFLX(kb).IJ(ii,1);
				j0=UFLX(kb).IJ(ii,2);
				dx=distance_spheric_coord(alat(j0,i0),elon(j0,i0),alat(j0,i0+1),elon(j0,i0+1));
        dnn=0;
				for k=1:l
					h3=thk(k,ii);  % cell thickness =average between i and i-1 
												 % (or j and j-1), gives H at U (or V) point
					if k==1, h3=h3+SSH(ii); end;
					A=h3*dx;
					D1(ii)=D1(ii)+h3;
					vf=sgn*V(k,ii)*A;
					VF(k,ii)=sgn*V(k,ii)*A;
					dnn=dnn+A;
				end
				AdZ(ii)=dnn; % area time depth at grid point ii
			 end
			 TRbtp=sgn*(AdZ.*UFLX(kb).UBRT);
			 Tbtp=nansum(TRbtp);
 % calc. flux through East. Carib.
			 if kb==2 & strncmp(xsct_name,'DomainOB',5) 
			 TrCB=nansum(TRbtp(1:60))*1e-6;
			 TrNA=nansum(TRbtp(61:end))*1e-6;
      end

      
      %keyboard
      Vav=Vav(1:l,:)+V.*thk;  % U*thickness
      VdT = Vav./sthk; % vav(layer=k) = sum(u(k)*thik(k))/sum(thk(time,k))
      VdT(end+1,:)=VdT(end,:);
      Vav(end+1,:)=Vav(end,:);
      UFLX(kb).Vav=Vav;
      UFLX(kb).THKav=sthk./cc;  % averaged layer thickness
      UFLX(kb).counter=cc;  % # of fields averaged

      Ltr=nansum(VF,2); % transport by layers
      TRP(kb).TrLayers(cc,:)=Ltr;
%      TRP(kb).TDav(cc,:)=nansum(VF,1);  % transport averaged over the total depth
      TRP(kb).VvelAv = VdT;    % mean U, l. thick is considered 
      Ttot=Ttot+nansum(TRP(kb).TrLayers(cc,:));
      TRP(kb).DepthPrf = Hs;
      TRP(kb).Time(cc) = dd;
      TRP(kb).TBtp(cc) = Tbtp; % transport from Barotropic U only
  %keyboard
      fprintf('====> KB=%i, NET Transport T=%5.2f Sv\n',kb, Ttot*1e-6);
      fprintf('====> Max Uav=%5.2f, Mean dP=%8.1f\n\n',...
	      nanmax(nanmax(VdT)),nanmean(nanmean(sthk/cc)));
    end;

    
    fprintf('====> cc=%i,   %s done, Transports (Sv): \n',cc,datestr(dd));
    fprintf('====> S: %5.2f; E: %5.2f; N: %5.2f; W: %5.2f \n',...
	    nansum(TRP(1).TrLayers(cc,:))*1e-6,...
	    nansum(TRP(2).TrLayers(cc,:))*1e-6,...
	    nansum(TRP(3).TrLayers(cc,:))*1e-6,...
	    nansum(TRP(4).TrLayers(cc,:))*1e-6);
    fprintf('=== Brtrp Trt S: %5.2f; E: %5.2f; N: %5.2f; W: %5.2f \n',...
	    TRP(1).TBtp(cc)*1e-6,TRP(2).TBtp(cc)*1e-6,...
	    TRP(3).TBtp(cc)*1e-6,TRP(4).TBtp(cc)*1e-6);
    if strncmp(xsct_name,'DomainOB',5) 
      fprintf('===   Brtp Trt East:   Caribbean = %4.1f, NAtl  = %4.1f\n',...
	      TrCB,TrNA);
    end
    fprintf(' Processing 1 rec: %6.4f min\n',toc/60);
    fprintf('------------------------------------\n\n');

  %  TD=nansum(VF);
%  keyboard

    if f_mat>0 & mod(cc,30)==0
      fprintf('Saving mat file: %s\n',fmat);
      save(fmat,'TRP','UFLX');
    end

  end;  % for ip
  
  if f_mat>0
    fprintf('Saving mat file: %s\n',fmat);
    save(fmat,'TRP','UFLX');
  end
      
elseif f_mat==3
% Loads 1st year only 
  iyr=YRPLT(1,1);
  
  fprintf('Loading mat file: %s\n',fmat);
  load(fmat);

end   % if f_mat
%exit

tp=-1.e-6*TRP(3).TBtp; % approximate Yuc. trpt










