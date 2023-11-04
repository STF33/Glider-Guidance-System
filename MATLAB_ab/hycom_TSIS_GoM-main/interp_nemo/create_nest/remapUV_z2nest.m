% After relax files with T/S fields are created
% see create_climat/relax_hycom_tsis.csh
% interpolate U & V into hycom vertical grid
%
% T and S (density) are interpolated into HYCOM
% layers using relax.com, this script
% dumps a nest files *.[ab]  with T,S,dP, etc and empty U,V
%
% U,V fields are interpolated into the vert. layers from the 
% Z-level files on NEMO vertical grid and HYCOM
% horizontal grid
%
% HYCOM can use mean file or instanteneous files
% in blkdat.input have to specify nestfq - nesting freq. (days)
% and bnstfq - barotrop. field nesting
% Mean archives: total vel = U_depth_av(barotropic)+U_depth_anomalies(bclinic)
%   "u-vel." and "v-vel." = total vel.
%   all variables are collocated: no C-grid, Land masking is 
%   the same for all variables and follows the P-grid
%
% Inst. archives: barotropic and baroclinic U,V not added to total U
% and "u-vel." is u baroclinic, "v-vel" is v baroclinic
%  Need to use different land mask for U,V, rho/dp 
%
% hycom 2.2.18 distiniguishes between mean and inst. archives
% by looking at fields inside *.b file
% "kemix", "kebtrop" are only in mean archive file
%
% Inst. archive files are created *archv*
% These will be used to generate restart files
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear
close all

dnmb = datenum(2011,7,9);  % interpolation date
DV = datevec(dnmb);

fprintf('Creating nest files archv %s\n\n',datestr(dnmb));

% For creating daily mean: need to average fields - archm files
% for instant. archive/nest files - no averaging,
% =1 - uvel = ubtrop+uanom (archm); =0 - uvel = uanom (arhv); 
f_avrg = 0; 

f_plt_section = 0;  % plot averaged fields, along OBs
esim   = '001';     
%esim   = '002';    % test 
thbase = 34.0;
thref  = 1e-3;   % reference specific vol, m3/kg
hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
rg     = 9806;
Nday   = 0;     % Averaging time, has to be an odd number!
Tv     = 00;    % IAS topography
gg     = 9.806;
ssh2m  = 1/gg; % convert HYCOM srfhgt to ssh in m
adj_btrop = 0;  % =0 - no btrop adjustmnt
                  % >0 - adjust U,V Btrop to match mean Yuc. transport

nlrs = 30;     % depth layers in HYCOM TSIS

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthintrp  = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';
pthrelax  = '/Net/kronos/ddmitry/hycom/TSIS/relax/';
pthmat    = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/datamat/';


%pthin0  = '/Net/kronos/morey/nest_files/realtime/';
%pthrlx0 = sprintf('/Net/kronos/ddmitry/hycom/goml0.04/relax/%s/SCRATCH/',esim); 
%pthnest0= sprintf('/nexsan/people/ddmitry/hycom/GOMl0.04/nest/%s/',esim);
%ptht    = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';

% HYCOM grid/ TOPO files:
% Read HYCOM topo:
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
HH  = -1*(nc_varget(ftopo,'mdepth'));
LAT = nc_varget(ftopo,'mplat');
LON = nc_varget(ftopo,'mplon');
HH(isnan(HH))=100;
[mm,nn]=size(HH);
[JDM,IDM] = size(HH);
IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


%ftopo = sprintf('%sdepth_GOMl0.04_%2.2i.a',pthtopo,Tv); % 
%flga  = sprintf('%sregional.grid.a',pthtopo);
%flgb  = sprintf('%sregional.grid.b',pthtopo);
%[HH,LAT,LON]  = read_topo(ftopo,flga,flgb);
if (min(min(HH))>0)
  error('Topo has to be <0 for ocean points ');
end


% Get land mask
% for u,v, p points:
pthnest = '/Net/kronos/ddmitry/hycom/TSIS/nest_files/';
fnsta = sprintf('%sarchv.2011_107_00.a',pthnest);
fnstb = sprintf('%sarchv.2011_107_00.b',pthnest);
LMASK = sub_lmask(fnsta,fnstb);

% Get Pbavg
%Pbavg = sub_pbavg;

iyr = DV(1);
iday = dnmb-datenum(iyr,1,1)+1;

fprintf('\n====> Processing %s\n\n',datestr(dnmb));

frlxa = sprintf('%srelax.m%2.2i%2.2i%4.4i.a',pthrelax,DV(3),DV(2),DV(1));
frlxb = sprintf('%srelax.m%2.2i%2.2i%4.4i.b',pthrelax,DV(3),DV(2),DV(1));

[TT,SS,dH] = sub_read_TSdP(frlxa,frlxb);
if isempty(TT),
  error('  Missing relax files for %s\n',datestr(dnmb));
end
[ll,mm,nn]=size(TT);
%  keyboard
%   ----------------  
if f_plt_section>0
  sttl0=sprintf('%4.4i/%2.2i/%2.2i',DV(1:3));
  fld='temp';
  f_lr=1;
  offst=0;
  sub_plot_section(TT,dH,fld,SCT,sttl0,f_lr,offst); 
  btmtxt='remapUV_z2nest.m';
  bottom_text(btmtxt);
  if s_fig>0
    for i=1:3
    ffg=sprintf('%snestSgm0_%s-%4.4i%2.2i%2.2i%2.2i-%i',...
  pthfig,fld,DV(1:4),i);
    fprintf('Saving %s\n',ffg);
    fll=sprintf('-f%i',i);
    print('-dpng',fll,'-r150',ffg);
    end
  end
%  keyboard
end
%   ----------------  

% Convert L. thickness into pressure units  
dP=dH*rg;  % m-> Pa (m*kg/m3*m/s2=kg/m*1/s2=N/m2=Pa)  

% for archm: Use averaged fields sub_average_TSdP->T,S,dP  
% INPUT files with T,S, etc fields
% interpolated into HYCOM grid
if f_avrg==1
  [Tav,Sav,dPav] = sub_average_TSdP(pthrlx0,dnmb0,Nday,hr1,hr2,dhr);
end

hr=0;
finaOLD = frlxa; 
finbOLD = frlxb;

% =====================
%
%     Interpolate U,V into HYCOM layers:
%
% In the mean archm files:
% u_vel=u_btrop+u_bclinic
%
% U and V in the nest files
% are on C-grid:
%     -------------
%     |           |
%     |           |
%     u    dP(i,j)  
%     |           |
%     |           |
%     ----- v ------
%
%  This is particularly important near the land and OBs
%
%  OBs: V are all 2^100 at South, East, North
%       U is 2^100 at East & North only
% =====================
% Read nest UV 
% and Average in time (if selected), z -level
if f_avrg==1
  [Uav,Vav,Eav,Zuv] = sub_avrgUV(pthin0,dnmb,Nday,dP,hr1,hr2,dhr);
end
  
%
% Get NEMO Z-levels:
% These are mid-grid depths
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);
ZMn = ZZN;  % NEMO
clear LONN LATN ZZN

nlrs = length(ZMn);

% Derive ZZ - interface layer depths
ZZn(1,1)=0;
for ik=1:nlrs
  zm1=ZMn(ik);
  dz=abs(2*(zm1-ZZn(ik)));
  ZZn(ik+1,1)=ZZn(ik,1)-dz;
end

%keyboard
% Interpolate into HYCOM hybrid-layers:
fuzlv = sprintf('%suvel_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthintrp,DV(3),DV(2),DV(1));
fvzlv = sprintf('%svvel_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthintrp,DV(3),DV(2),DV(1));

[UU,VV]=sub_interpUV2sgmlrs(fuzlv,fvzlv,ZZn,ZMn,frlxa,frlxb,HH);
%keyboard

f_chck=0;
if f_chck==1
  ipt=5;
  uu=squeeze(UU(ipt,:,:));
  vv=squeeze(VV(ipt,:,:));
  ss=sqrt(uu.^2+vv.^2);

  
  figure(10); clf;
  pcolor(ss); shading flat;
  title(sprintf('HYCOM layer %i',ipt));
  colorbar;
end

%keyboard

% 
% Make C-grid and nan's at OB where needed
if f_avrg==0
  [UU,VV]=sub_putUVonCgrid(UU,VV,HH);
end

%
% Calculate depth-average U,V on C grid
% For the mean archive, u barotropic is used
% to derive u baroclinic (see hycom source: forfun.f - for 3D nest)
% meanar >0: unest=unest - util1, util1=u_btrop
% Looking for "kemix" field in the nest archive file (our case)
%   [otherwise, this is not mean archive and
%    and btrop and bclinic separated]   
% Make sure to write u_brtrop, v_bartrop max/min in *b file
% otherwise (if both are 0s), 
% rd_archive (forfun.f) makes whole field=0
% Ubrtop and Vbrtrop are used at the boundaries (as well as surf. elev)
% NOTE: Barotrop. U,V are not nan's at the OBs (unlike baroclinic 3D)
% Split barotropic and baroclinic parts:
[Ubt,Vbt] = sub_UVbtrop(UU,VV,dH);
Ubt(isnan(Ubt))=0;
Vbt(isnan(Vbt))=0;
if f_avrg==0
  [UU,VV] = sub_Utot2bcl(UU,VV,Ubt,Vbt);
end

% Adjust barotropic components to
% match Yucatan transport if needed
%  if adj_btrop>0;
%    [UU,VV] = sub_adjust_Ubtrop(UU,VV,dH,adj_btrop);    
%  end
  
%
% Add ssh to archv file
fssh = sprintf('%sssh_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
                 pthintrp,DV(3),DV(2),DV(1));
fprintf('Reading %s\n',fssh);
fid=fopen(fssh,'r');
dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
ssh = reshape(dmm,IDM,JDM);
ssh = ssh';
fclose(fid);
I=find(HH>=0);
ssh(I)=nan;

% Calculate Montgomery potential from ssh
fmat_nest=sprintf('%sregression_montg_ssh_IAS003.mat',pthmat);
montg1 = sub_montg1_regr2(ssh,fmat_nest);
  

%keyboard

% ========================================  
%      WRITING OUTPUT files:
%
finaNEW = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthoutp,iyr,iday);
finbNEW = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthoutp,iyr,iday);

% Open files "relax" files (anom+mean old nest) for reading
% and nest files for writing
disp('Openning files ...');
fprintf('relax files: %s\n %s\n',finaOLD,finbOLD);
fprintf('new HYCOM archive: %s\n %s\n',finaNEW,finbNEW);

fida_OLD = fopen(finaOLD,'r');  %
fidb_OLD = fopen(finbOLD,'rt'); %
fida_NEW = fopen(finaNEW,'w');
fidb_NEW = fopen(finbNEW,'wt');

%
% Get hycom model date:
[dS,dE,dC] = get_dtime(dnmb);
dnmb_hycom = dC;

nlevU=0;  % level counter u-vel
nlevV=0;  % level counter v-vel
nlevK=0;  % level counter KE
nlevD=0;  % level counter layer thickness 
nlevT=0;  % level counter layer T 
nlevS=0;  % level counter layer S 
nlevR=0;  % level counter Density anomaly (s.dens - thbase=34)

if ~exist('TDENS','var')
  TDENS=read_targ_dens(finbOLD,'salin','=',4);
end

% Write header:
aa1=sprintf('Daily mean: NEMO 1km + GLORYS (ouside GoM) %2.2i/%2.2i/%4.4i to HYCOM',...
      DV(3),DV(2),DV(1));
aa2='Sigma2*; GDEM3i Jan. init; KPP mixed layer; SeaWiFS mon KPAR; energy-loan ice;';

for nl=1:10
  aa=fgetl(fidb_OLD);
  if nl==1
    aa=aa1;
  end
  if nl==2;
    aa=aa2;
  end
  fprintf(fidb_NEW,[aa,'\n']);

  if nl==8
    dmm=aa(2:8);
    ID=str2num(dmm);
  end
  if nl==9
    dmm=aa(2:8);
    JD=str2num(dmm);
  end

  if exist('JD','var')
    if (ID~=nn | JD~=mm)
      fprintf('ERR: Check dimensions in relax-nest and old nest:\n');
      fprintf('RLX: ID=%i, JD=%i; NEST: ID=%i, JD=%i\n',nn,mm,ID,JD);
      error('*** STOPPING ***');
    end
  end
 
end

%  IJDM=ID*JD;
%  npad=4096-mod(IJDM,4096);
%  toto=ones(npad,1);

%keyboard
% ---------------------------------------------
% Write surface and barotropic part with no change:
% Barotropic fields are used in nest files
% depending on the flags in blkdat.input:
% If barotropic nesting is used then:
% 1.0  'bnstfq' = number of days between baro nesting archive input
% 2  'lbflag' = lateral barotropic bndy flag (0=none, 1=port, 2=input)
% Need: surf elevation, u,v barotropic and montg1
%  
% montg1 is estimated from:
% montg1 = g*zeta - thref*Pb 
%  note that zeta is in m (ssh) and thref=1e-3 m3/kg spec.vol.
%  Pb - barotropic pressure
% ---------------------------------------------
tic;
while ischar(aa)
  aa=fgetl(fidb_OLD);

  if ~ischar(aa), break; end;

  I=strfind(aa,'=');
  bb=aa(I+1:end);
  S=sscanf(bb,'%f');

  if S(3)==0  %# hybr. layer #
    fprintf('<= %s\n',aa);
%keyboard
% Read 1 record from HYCOM relax file:    
    dmm  = fread(fida_OLD,IJDM,'float32','ieee-be');  % read field
    toto = fread(fida_OLD,npad,'float32','ieee-be');  % read npad 
% Save land mask:
%      if ~exist('lmask','var')
%        lmask=find(dmm>hg/10);
%      end
%keyboard
% Barotropic fields need to be changed:
    Imt = strmatch('montg',aa);
    LL  = strmatch('srfhgt',aa);
    II  = strmatch('u_btrop',aa);
    JJ  = strmatch('v_btrop',aa);
    LS1 = strmatch('surflx',aa);
    LS2 = strmatch('salflx',aa);
    LS3 = strmatch('bl_dpth',aa);
    LS4 = strmatch('mix_dpth',aa);

% Barotrop fields not needed in archv
    Imx = strmatch('vmix',aa);
    Smx = strmatch('smix',aa);
    Tmx = strmatch('tmix',aa);
    Thmx= strmatch('thmix',aa);
    Umx = strmatch('umix',aa);
    Vmx = strmatch('vmix',aa);

    if ~isempty(Imx) | ~isempty(Smx) | ~isempty(Tmx) | ...
       ~isempty(Thmx)| ~isempty(Umx) | ~isempty(Vmx)
      fprintf('Skipping field ...\n');
      continue;
    end

% Write Montgomery potential:
    if ~isempty(Imt),
      lmask = LMASK.temp;
      dmm = reshape(montg1',IJDM,1);
      I=find(isnan(dmm));
      dmm(I)=0;
      dmm(lmask)=hg;
      Inh=find(dmm<hg/10);
      minv=min(dmm(Inh));
      maxv=max(dmm(Inh));
      nlev=4;  % flag for correct converting to restart
      td0=34;
% Parse string from old nest *b:
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);
      aa = aa2;
    elseif ~isempty(II),  % write barotropic U
%        lmask = LMASK.uvel;
      dmm = reshape(Ubt',IJDM,1);
      I=find(isnan(dmm));
      dmm(I)=0;
%        dmm(lmask)=0;
      Inh=find(dmm<hg/10);
      minv=min(dmm(Inh));
      maxv=max(dmm(Inh));
      nlev=0;
      td0=0;
% Parse string from old nest *b:
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);
      aa = aa2;
    elseif ~isempty(JJ)
%        lmask = LMASK.vvel;
      dmm = reshape(Vbt',IJDM,1);
      I=find(isnan(dmm));
      dmm(I)=0;
%        dmm(lmask)=0;
      Inh=find(dmm<hg/10);
      minv=min(dmm(Inh));
      maxv=max(dmm(Inh));
      nlev=0;
      td0=0;
% Parse string from old nest *b:
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);
      aa = aa2;
    elseif ~isempty(LL),  % SSH
      lmask = LMASK.thknss;
      SSH=1/ssh2m*ssh;  % SSH m -> srfhgt units
      dmm = reshape(SSH',IJDM,1);
      I=find(isnan(dmm));
      dmm(I)=0;
      dmm(lmask)=hg;
      Inh=find(dmm<hg/10);
      minv=min(dmm(Inh));
      maxv=max(dmm(Inh));
      nlev=0;
      td0=0;
% Parse string from old nest *b:
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);
      aa = aa2;
% keyboard
    elseif ~isempty(LS1) |  ~isempty(LS2) | ~isempty(LS3) | ~isempty(LS4)
      nlev=0;
      td0=0;
      minv = [];
      maxv = [];
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);
      aa = aa2;
    end 

% Write in the new nest file:
    fwrite(fida_NEW,dmm,'float32','ieee-be');
    fwrite(fida_NEW,toto,'float32','ieee-be');

    fprintf('=> %s\n',aa);
    fprintf(fidb_NEW,[aa,'\n']);
   
    if f_avrg==1 
      if ~isempty(JJ)
        Fdmm=sqrt(Ubt.^2+Vbt.^2); % this is not right: U,V are on C-grid
        fldnm='kebtrop';          %but HYCOM does not read this field
        sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);
      end

      if ~isempty(Imx); % add missing fields to *b and *a, mean archive
        Fdmm=zeros(JD,ID);
        fldnm='kemix';
        sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);

        fldnm='covice';
        sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);

        fldnm='thkice';
        sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);

        fldnm='temice';
        sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);
      end;
    end

    continue
    
  end  % if S  - surface and barotropic fields are done
  
  bb=aa(1:I-1);
  wrd = deblank(sscanf(bb,'%s'));
  Imx = strmatch('v-vel.',aa);
  
  switch (wrd)
% -----------------
%  u-vel
% -----------------
   case('u-vel.')
    nlevU=nlevU+1;
    nlev=nlevU;
    
    lmask = LMASK.uvel;
    A = squeeze(UU(nlev,:,:));
    A = A';
    dmm = reshape(A,IJDM,1);
    I=find(isnan(dmm));
%    umm=nanmean(dmm);
    dmm(I)=0;
    dmm(lmask)=hg;
%keyboard
% Keep U for k.e.:
%      uvel = dmm;
    
% Write *a file:           
%          dmm(lmask)=hg;
    fwrite(fida_NEW,dmm,'float32','ieee-be');
    fwrite(fida_NEW,toto,'float32','ieee-be');
    
% Update *b file:      
    Inh=find(dmm<hg/10);
    minv=min(dmm(Inh));
    maxv=max(dmm(Inh));

    ir = strfind(aa,'=');
    SD = sscanf(aa(ir+1:end),'%f');
    td0=TDENS(nlev);

% Parse string from old nest *b:
    aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
    fprintf('<= %s\n',aa);
    fprintf('=> %s\n',aa2);
    fprintf(fidb_NEW,[aa2,'\n']);

%keyboard
% -----------------
%  v-vel
% -----------------
   case('v-vel.')
    nlevV=nlevV+1;
    nlev=nlevV;
    
    lmask=LMASK.vvel;
    A = squeeze(VV(nlev,:,:));
    A = A';
    dmm = reshape(A,IJDM,1);
    I=find(isnan(dmm));
%    umm=nanmean(dmm);
    dmm(I)=0;
    dmm(lmask)=hg;
% Keep V for k.e.:
%      vvel = dmm;
    
% Write *a file:           
%          dmm(lmask)=hg;
    fwrite(fida_NEW,dmm,'float32','ieee-be');
    fwrite(fida_NEW,toto,'float32','ieee-be');
    
% Update *b file:      
    Inh=find(dmm<hg/10);
    minv=min(dmm(Inh));
    maxv=max(dmm(Inh));

    ir = strfind(aa,'=');
    SD = sscanf(aa(ir+1:end),'%f');
    td0=TDENS(nlev);
   
% Parse string from old nest *b:
    aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
    fprintf('<= %s\n',aa);
    fprintf('=> %s\n',aa2);
    fprintf(fidb_NEW,[aa2,'\n']);

    if f_avrg==1 
% ADD missing fields: k.e.
      Fdmm=zeros(JD,ID);
      fldnm='k.e.';
      sub_add_fields2archv(fldnm,Fdmm,fida_NEW,fidb_NEW,aa,toto);
    end

% --------------------------      
   case('thknss')
    nlevD=nlevD+1;                 % level # counter
    nlev=nlevD;
    
    lmask = LMASK.thknss;
    A = squeeze(dP(nlev,:,:));
    A = A';
    dmm = reshape(A,IJDM,1);
    I=find(isnan(dmm));
    umm=nanmean(dmm);
    dmm(I)=umm;
    dmm(lmask)=hg;

    if max(size(dmm))~=IJDM, 
      error('Thickness size is not IJDM');
    end
 
% Write *a file:           
%          dmm(lmask)=hg;
    fwrite(fida_NEW,dmm,'float32','ieee-be');
    fwrite(fida_NEW,toto,'float32','ieee-be');

% Update *b file:      
    Inh=find(dmm<hg/10);
    minv=min(dmm(Inh));
    maxv=max(dmm(Inh));

    ir = strfind(aa,'=');
    SD = sscanf(aa(ir+1:end),'%f');
    td0=TDENS(nlev);
    
% Parse string from old nest *b:
    aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
    fprintf('<= %s\n',aa);
    fprintf('=> %s\n',aa2);
    fprintf(fidb_NEW,[aa2,'\n']);
% Keep str:
    Dstr = aa2;

% ----------------
%    S
% ----------------
  case('salin')
   nlevS=nlevS+1;                 % level # counter
   nlev=nlevS;
   
   lmask = LMASK.temp;
   A = squeeze(SS(nlev,:,:));
   A = A';
   dmm = reshape(A,IJDM,1);
   Slr=dmm;
   I=find(isnan(dmm));
   umm=nanmean(dmm);
   dmm(I)=umm;
   dmm(lmask)=hg;

   if max(size(dmm))~=IJDM, 
     error('Thickness size is not IJDM');
   end

% Write *a file:           
%          dmm(lmask)=hg;
   fwrite(fida_NEW,dmm,'float32','ieee-be');
   fwrite(fida_NEW,toto,'float32','ieee-be');

% Update *b file:      
   Inh=find(dmm<hg/10);
   minv=min(dmm(Inh));
   maxv=max(dmm(Inh));

   ir = strfind(aa,'=');
   SD = sscanf(aa(ir+1:end),'%f');
   td0=TDENS(nlev);

% Parse string from old nest *b:
   aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
   fprintf('<= %s\n',aa);
   fprintf('=> %s\n',aa2);
   fprintf(fidb_NEW,[aa2,'\n']);
% Keep str:
   Dstr = aa2;

% ----------------- -----------------
% Layer density anomaly: sig2 - thbase
% ----------------- -----------------
  case('density')
   if f_avrg==0; continue; end;
   nlevR=nlevR+1;
   nlev=nlevR;

   if nlevR~=nlevT | nlevR~=nlevS
     error('*** ERR: need to get T and S before calculating Rho');
   end
%keyboard
% Calculate sigma0 using hycom function:
   lmask = LMASK.temp;
   Slr(lmask)=nan;
   Tlr(lmask)=nan;
   sgm0_f=sigma0_hycom2218(Slr,Tlr)-thbase;  
   sgm0_f(lmask)=hg;
   dmm=sgm0_f;
% Write *a file:           
   fwrite(fida_NEW,dmm,'float32','ieee-be');
   fwrite(fida_NEW,toto,'float32','ieee-be');

% Update *b file:      
   Inh=find(dmm<hg/10);
   minv=min(dmm(Inh));
   maxv=max(dmm(Inh));

   ir = strfind(aa,'=');
   SD = sscanf(aa(ir+1:end),'%f');
   td0=TDENS(nlev);

% Parse string from old nest *b:
   aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
   fprintf('<= %s\n',aa);
   fprintf('=> %s\n',aa2);
   fprintf(fidb_NEW,[aa2,'\n']);
% Keep String:
   Rstr=aa2;
   
   
  case('temp')
   nlevT=nlevT+1;                 % level # counter
   nlev=nlevT;

   lmask = LMASK.temp;
   A = squeeze(TT(nlev,:,:));
   A = A';
   dmm = reshape(A,IJDM,1);
   Tlr = dmm;
   I=find(isnan(dmm));
   umm=nanmean(dmm);
   dmm(I)=umm;
   dmm(lmask)=hg;
   
   if max(size(dmm))~=IJDM, 
     error('Thickness size is not IJDM');
   end

% Write *a file:           
%          dmm(lmask)=hg;
   fwrite(fida_NEW,dmm,'float32','ieee-be');
   fwrite(fida_NEW,toto,'float32','ieee-be');

% Update *b file:      
   Inh=find(dmm<hg/10);
   minv=min(dmm(Inh));
   maxv=max(dmm(Inh));

   ir = strfind(aa,'=');
   SD = sscanf(aa(ir+1:end),'%f');
   td0=TDENS(nlev);


% Parse string from old nest *b:
   aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv,'hycom_date',dnmb_hycom);

% Write *b file:
   fprintf('<= %s\n',aa);
   fprintf('=> %s\n',aa2);
   fprintf(fidb_NEW,[aa2,'\n']);
% Keep str:
   Dstr = aa2;
   
   
 end;  % switch wrd - field name   

  
  
end;  % while aa - reading *b file all fields

fclose(fida_NEW);
fclose(fidb_NEW);
fclose(fida_OLD);
fclose(fidb_OLD);

fprintf('  ####    1 record: %9.2f sec\n\n',toc);
fprintf(' OUTPUT saved: %s\n',finaNEW);
fprintf(' OUTPUT saved: %s\n\n',finbNEW);


  
