% 0.04 HYCOM GOMu reanalysis 50.1
% HYCOM + NCODA Gulf of Mexico 1/25° Reanalysis
% Tidal forcing is included!
%
% netcdf, on z levels, collocated
%
% Calculate Mass flux budget in a closed volume
% - saves transport by layers
%
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_gom04/plot_binary_output;
startup

clear all
close all

HR=0;
esim='GoMu501'; % 

% Years to calculate:
ys = 2011;
ye = 2011;

pltx=0;
sfig=0;
f_mat = 1;  % = 0 - do not save,
            % = 1 - save matlab file with transport, etc.;
            % = 2 - load saved mat file and continue from where the last day is 
            %      in the mat file till the last day = dy2

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
pthtopo = '/nexsan/archive/GOMu0.04_501/topo/';
btx = 'calc_volFlx_GoMreanalysis.m';

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
ftopo = sprintf('%sdepth_GOMu0.04_03i.nc',pthtopo);
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


% Daily average fields to get rid of tides:
	ihr = 0;
	[UFLXd,TRPd] = sub_UFLX_TRPcollocated(BB);

  tic;
	
	for HR=0:3:21
%    pthi=sprintf('/nexsan/GLBu0.08/GLBu0.08_191/data/%4.4i/',iyr);
		pthi=sprintf('/nexsan/archive/GOMu0.04_501/data/netcdf/%4.4i/',iyr);
		fina=sprintf('%shycom_gomu_501_%i%2.2i%2.2i00_t%3.3i.nc',...
								pthi,DV(1:3),HR);
	
		if ~exist(fina,'file')
			fprintf('Missing %s\n',fina);
			continue
		end

    ihr = ihr+1;
    dnmbh = datenum(DV(1),DV(2),DV(3),HR,0,0);	
		[UFLXd,TRPd] = sub_calcUFLX_znc(fina,UFLXd,NB,TRPd,ihr,alat,elon,dnmbh);

		for kb=1:NB
      if ihr==1
  			SM(kb).UBRT   = UFLXd(kb).UBRT*0;
	  		SM(kb).Vav    = UFLXd(kb).Vav*0;
		  	SM(kb).VvelAv = TRPd(kb).VvelAv*0;
			end
			SM(kb).UBRT     = SM(kb).UBRT+UFLXd(kb).UBRT;
			SM(kb).Vav      = UFLXd(kb).Vav;
%      SM(kb).VvelAv   = SM(kb).VvelAv+TRPd(kb).VvelAv;
    end
	end

% Daily average:
  nhrs = UFLXd(1).counter;
  if nhrs<3
    continue;
  end
  for kb=1:NB
    SM(kb).TrLayers = TRPd(kb).TrLayers(1,:)*0;
    SM(kb).TBtp     = 0;
  end

  for kb=1:NB
    for ihr=1:nhrs
      SM(kb).TrLayers = SM(kb).TrLayers+TRPd(kb).TrLayers(ihr,:);
      SM(kb).TBtp     = SM(kb).TBtp+TRPd(kb).TBtp(ihr);
    end
  end


	cc=cc+1;
  Ttot=0;
  fprintf('\n   =====      Daily Mean    =====\n');
  for kb=1:NB
% Daily means:
    Vav   = SM(kb).Vav/nhrs;
    UBRT  = SM(kb).UBRT/nhrs;
    TrLayers = SM(kb).TrLayers/nhrs;
    TBtp  = SM(kb).TBtp/nhrs;
    dZ    = UFLXd(kb).dZ;
    Hs    = UFLXd(kb).Hs;

    if ~isfield(UFLX,'Vav') | cc==1;  UFLX(kb).Vav=Vav*0; end;  
		UFLX(kb).Vav=UFLX(kb).Vav+Vav;
		UFLX(kb).dZ = dZ;
		UFLX(kb).counter=cc;  % # of fields averaged
    UFLX(kb).Dst = UFLXd(kb).Dst;
    UFLX(kb).ZZ = UFLXd(kb).ZZ;
    Vav_all     = UFLX(kb).Vav/cc;  % overall mean U section

		TRP(kb).VvelAv = Vav_all;    % mean U section
		Ttot=nansum(TrLayers);
    TRP(kb).TrLayers(cc,:) = TrLayers;
		TRP(kb).DepthPrf = Hs;
		TRP(kb).Time(cc) = dd;
		TRP(kb).TBtp(cc) = TBtp; % transport from Barotropic U only
%keyboard
		fprintf('====> KB=%i, NET Transport T=%5.2f Sv\n',kb, Ttot*1e-6);
		fprintf('====> Max Uav=%5.2f \n',...
    nanmax(nanmax(abs(Vav_all))));
	end;
%keyboard
	
	
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
	fprintf(' Processing 1 day: %6.4f min\n',toc/60);
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
		









