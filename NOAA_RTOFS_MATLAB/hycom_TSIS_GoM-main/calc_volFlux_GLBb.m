% HYCOM GLBb reanalysis 19.1
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
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

HR=0;
esim='GLBu191'; % 

% Years to calculate:
ys = 2011;
ye = 2011;

pltx=0;
sfig=0;
f_mat = 1;  % = 0 - do not save,
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



rg=9806;  % convert pressure to depth, m
huge=1e10;

mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(ys,4)==0; mday(2)=29; end;
cc=0;
for iyr=ys:ye
  id1=1;
  id2=365;
  ic=mod(iyr,4);
  if ic==0 & id2==365, id2=366; end;
  if iyr>ys
    id1=1;
  end

  for iday=id1:id2
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
  end
end


% -------------------------
% Directories:
% -------------------------
ptht   = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
pthmat = '/Net/mars/ddmitry/hycom/hycom_TSIS/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/hycom_TSIS/fig/';
pthtopo = '/nexsan/archive/GLBu0.08_191/topo/';
btx = 'calc_volFlx_GLBb.m';


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
ftopo = sprintf('%sdepth_GLBu0.08_07b.nc',pthtopo);
HH  = -1*squeeze(nc_varget(ftopo,'depth'));
HH(isnan(HH))=100;
alat = nc_varget(ftopo,'Latitude');
elon = nc_varget(ftopo,'Longitude');
[mm,nn]=size(HH);
JD=mm;
ID=nn;
m=mm;
n=nn;
HH(isnan(HH))=100;
IJDM = ID*JD;


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
    SWE = [-86.83, 21.32; -84.24, 21.32]; % Weast- East section
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
end;


% Find indices:
if ~exist('iSW','var');
  clear IY
  for k=1:2
    if sign(elon(1))~=sign(SWE(k,1)); 
      error('ERR: Need to convert longitudes'); 
    end;
    dLG=abs(elon-SWE(k,1));
    dLT=abs(alat-SWE(k,2));
    i=find(dLG==min(dLG),1);
    j=find(dLT==min(dLT),1);
    if elon(i)>SWE(1,1) & elon(i)<SWE(2,1),
      if k==1,
 i=i-1;
      else
 i=i+1;
      end
    end
    
      
%    DST = distance_spheric_coord(alat,elon,SWE(k,2),SWE(k,1));
%    [j,i]=find(DST==min(min(DST)));
    IY(k,:)=[i,j];
  end
  dLG=abs(elon-SSN(2,1));
  dLT=abs(alat-SSN(2,2));
  i=find(dLG==min(dLG),1);
  j=find(dLT==min(dLT),1);
%  DST = distance_spheric_coord(alat,elon,SSN(2,2),SSN(2,1));
%  [j,i]=find(DST==min(min(DST)));
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
  BB.fplt= 0; % =0 - do not plot transect
  [UFLX,TRP] = sub_UFLX_TRPcollocated(BB);
  
  
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
      [UFLX,TRP] = sub_UFLX_TRPcollocated(BB);
      
    end


    pthi=sprintf('/nexsan/archive/GLBu0.08_191/data/%4.4i/',iyr);
    fina=sprintf('%shycom_GLBu0.08_191_%i%2.2i%2.2i00_t%3.3i.nc',...
   pthi,DV(1:3),HR);
    
    if ~exist(fina,'file')
      fprintf('Missing %s\n',fina);
      continue
    end

    cc=cc+1;
    tic;
    
    SSH=squeeze(nc_varget(fina,'surf_el'));

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

  %  =================
  % Create array of interface depths from layer thicknesses:
  %  =================
    ZZ=nc_varget(fina,'depth');
    nl=length(ZZ);
    dZ(1,1)=0.5*(ZZ(2)-ZZ(1));
    for kk=2:nl-1
      dZ(kk,1)=0.5*abs(ZZ(kk+1)-ZZ(kk-1));
    end
    dZ(nl)=0.5*abs(ZZ(nl)-ZZ(nl-1));

    for ks=1:NB
      UFLX(ks).ZZ=ZZ;
    end


  % Note that in the raw binary, U, V are on u,v grid
  % h is on p-grid
  % !!! NOTE: in the instanteneous output files archv u-vel, v-vel 
  % are BAROCLINIC (not total) components
  % need to add u barotropic !!!
  %  In nc z-interpolated U,V H collocated
  % U = total U
  %keyboard

  % Get velocities along the sides of the box
    for kb=1:NB
      OB = UFLX(kb).OB;
      js1=UFLX(kb).IJ(1,2);
      js2=UFLX(kb).IJ(end,2);
      is1=UFLX(kb).IJ(1,1);
      is2=UFLX(kb).IJ(end,1);

      dj=js2-js1+1;
      di=is2-is1+1;
      U = squeeze(nc_varget(fina,'water_u',...
           [0 0 js1-1 is1-1],[1 nl dj di]));
      V = squeeze(nc_varget(fina,'water_v',...
           [0 0 js1-1 is1-1],[1 nl dj di]));
      
      % Barotropic U & V:
      dmm=U';
      dmm(isnan(dmm))=0;
      I=dmm*0+1;
      I(dmm==0)=0;
      Hbz=I*dZ; % Bottom depth in z cells
      Ubrt=dmm*dZ./Hbz;
      
      dmm=V';
      dmm(isnan(dmm))=0;
      I=dmm*0+1;
      I(dmm==0)=0;
      Vbrt=dmm*dZ./Hbz;

      UFLX(kb).sumDZ=Hbz; % integrated Z cells to bottom 
      
      if strcmp(OB,'S') | strcmp(OB,'N')
        UFLX(kb).V=V;
        UFLX(kb).UBRT=Vbrt';
      elseif strcmp(OB,'E') | strcmp(OB,'W')
        UFLX(kb).V=U;
        UFLX(kb).UBRT=Ubrt';
      end

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
    Ttot=0;
    for kb=1:NB
      js1 = UFLX(kb).IJ(1,2);
      js2 = UFLX(kb).IJ(end,2);
      is1 = UFLX(kb).IJ(1,1);
      is2 = UFLX(kb).IJ(end,1);

      ni1 = is2-is1+1;
      nj1 = js2-js1+1;

      Hs  = UFLX(kb).Hs;
%      thk = UFLX(kb).thk;
      V   = UFLX(kb).V;
      SSH = UFLX(kb).SSH;
      Vav = UFLX(kb).Vav;
      Ubt = UFLX(kb).UBRT;
%      sthk= UFLX(kb).sumthk;
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

      Dst=UFLX(kb).Dst;
      if isempty(Dst)
        Dst=zeros(ni,1);
      end
  %keyboard
  % ===========================
  % FLUX CALCULATION:
  % ===========================
     clear AdZ
      for ii=1:ni
        i0=UFLX(kb).IJ(ii,1);
        j0=UFLX(kb).IJ(ii,2);
        if ii<ni
          i1=UFLX(kb).IJ(ii+1,1);
          j1=UFLX(kb).IJ(ii+1,2);
        else
          i1=UFLX(kb).IJ(ii-1,1);
          j1=UFLX(kb).IJ(ii-1,2);
          
        end
          
        if Dst(ni)==0;
          dx=distance_spheric_coord(alat(j0),elon(i0),...
             alat(j1),elon(i1));
          Dst(ii)=dx;
        end

        if ii>1 & ii<ni
          dx=Dst(ii);
        else
          dx=0.5*dx;
        end
        
        dnn=0;
        vv=V(:,ii);
        dz=dZ;
        dz(1)=dz(1)+SSH(ii);
        VF(:,ii)=sgn*vv.*dz.*dx;
        in=find(isnan(vv));
        dz(in)=nan;
        AdZ(ii)=nansum(dz.*dx); % area of w colmn * dx
      end
      TRbtp=sgn*(AdZ.*UFLX(kb).UBRT);
      Tbtp=nansum(TRbtp);
% calc. flux through East. Carib.
      if kb==2 & strncmp(xsct_name,'DomainOB',5) 
        TrCB=nansum(TRbtp(1:60))*1e-6;
        TrNA=nansum(TRbtp(61:end))*1e-6;
      end

      
      for ii=1:ni
        Vav(1:nl,ii)=Vav(1:nl,ii)+V(:,ii);  % 
      end
     
      if isempty(UFLX(kb).Dst)
        UFLX(kb).Dst = Dst;
      end 
      UFLX(kb).Vav=Vav;
      UFLX(kb).dZ = dZ;
%      UFLX(kb).THKav=sthk;  % averaged layer thickness
      UFLX(kb).counter=cc;  % # of fields averaged

      TRP(kb).VvelAv = Vav/cc;    % mean U vertical cross-section
      TRP(kb).TrLayers(cc,:)=nansum(VF,2);
      Ttot=Ttot+nansum(TRP(kb).TrLayers(cc,:));
      TRP(kb).DepthPrf = Hs;
      TRP(kb).Time(cc) = dd;
      TRP(kb).TBtp(cc) = Tbtp; % transport from Barotropic U only
  %keyboard
      fprintf('====> KB=%i, NET Transport T=%5.2f Sv\n',kb, Ttot*1e-6);
      fprintf('====> Max Uav=%5.2f \n\n',...
       nanmax(nanmax(abs(TRP(kb).VvelAv))));
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
  fmat = sprintf('%s%s_VolTrtBoxi%3.3ij%3.3i_%4.4i-%i.mat',...
        pthmat,esim,iSW,jSW,iyr,CYC);
  
  fprintf('Loading mat file: %s\n',fmat);
  load(fmat);

end   % if f_mat
%exit

tp=-1.e-6*TRP(3).TBtp; % approximate Yuc. trpt










