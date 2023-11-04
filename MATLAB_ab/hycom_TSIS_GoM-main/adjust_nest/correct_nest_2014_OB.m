% Nest files prepared from Global HYCOM
% Years 2014-2015
%
% 
% Need to correct nest files: correct_nest_OB.m
% 1) This code corrects created nest files
%    by nest.com in TSIS directory
% /home/ddmitry/codes/HYCOM_TSIS/
%
% 2) Montgomery potential is corrected. It is made corresponding
%    the ssh field. montg1 is approximated by regressing
%    mntg on ssh <--- Caused neg depths, do not do this for now
%
% Only existing nest files are corrected, missing dates - skipped
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

close all
clear

yr1=2014;
yr2=2015;

pthnest0 = '/Net/kronos/ddmitry/hycom/TSIS/';
%pthtmp   = '/nexsan/people/ddmitry/hycom/GOMl0.04/nest/080/old/';
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';
%pthmat = '/nexsan/people/ddmitry/hycom/GOMl0.04/nest/data_mat/';

hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
rg     = 9806;

adj_btrop = 2.5;  % gives ~31.4 Sv during 4 mo, yr 2010
                   % =0 - no btrop adjustmnt
                   % >0- <1 - adjust U,V Btrop to match mean Yuc. transport
%
% Time-varying scale factor
%f_adj = sprintf('%sscale_factor.mat',pthmat);
%ADJ = load(f_adj);
%TMa = ADJ.tvec;
%UBscale = ADJ.sclfac_sm;

fprintf('\n  Correcting nest files: adjust btrop=%3.2f, %i-%i\n\n',...
	adj_btrop,yr1,yr2);

% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
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

%[DX,DY]=sub_dx_dy(elon,alat);


cc=0;
for iy=yr1:yr2
  ic=mod(iy,4);
  d1=1; % 
  d2=365;
  if ic==0, d2=366; end;
  for iday=d1:d2
%  for iday=1:d2
    pthnest022 = sprintf('%s/nest_files_2014_2015/%i/',pthnest0,iy);
    fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthnest022,iy,iday);

    if ~exist(fina,'file'); continue; end;

    cc=cc+1;
    YRPLT(cc,1)=iy;
    YRPLT(cc,2)=iday;
  end
end

% Get land mask
% for C-grid
% for u,v, p,T,S points:
iyr=iy;
LMASK = sub_lmask(YRPLT(1,1),YRPLT(1,2));

iyr0=0;
nrc=size(YRPLT,1);
for ip=1:nrc
  tic;
  iyr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  d1=datenum(iyr,1,1);
  dnmb=d1+iday-1;
  DV=datevec(dnmb);
  
  fprintf('Correcting files from YEAR %i\n',iyr);
  
  pthnest = sprintf('%snest_adjusted/',pthnest0);

% OLD files - with errors in U,V at N and E OB
% NEW files - corrected
  pthtmp  = sprintf('%s/nest_files_2014_2015/%i/',pthnest0,iyr);
  finbOLD = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthtmp,iyr,iday);
  finaOLD = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthtmp,iyr,iday);

  if ~exist(finaOLD,'file')
    continue;
  end
 
  fprintf('OLD nest: %s\n',finaOLD); 
  fprintf('OLD nest: %s\n',finbOLD); 
 
%  % Copy old nests to tmp
%  if iyr~=iyr0
%    scmd=sprintf('csh mv_old.csh %i',iyr);
%    system(scmd);
%    iyr0=iyr;
%  end

% Layers:
%  [F,n,m,l] = read_hycom(finaOLD,finbOLD,'thknss');
%  F=F./rg;  % Layer thickness, m
%  F(F>1e20)=0;
%  dH=F;

% Read in old U,V fields
% correct N and E OB if needed
%  [UU,VV] = sub_correctNE_OB(finaOLD,finbOLD); % corrected brcl U,V

% No need to change baroclinic components
%[F,n,m,l] = read_hycom(finaOLD,finbOLD,'u-vel.');
%UU=F;
%[F,n,m,l] = read_hycom(finaOLD,finbOLD,'v-vel.');
%VV=F;

% Note this nest files
% U, V are not collocated and u,v-vel are baroclinic comp
% without U-btrop added
%clc=0;
%utotal = check_totalU(finaOLD,finbOLD,clc);

% Adjust barotropic components to
% match Yucatan transport if needed
%  iadj=find(TMa==dnmb);
%  if isempty(iadj),
%    error('Cannot find adj_btrop time');
%  end
%  adj_btrop=UBscale(iadj);
%  keyboard
  [F,n,m] = read_hycom(finaOLD,finbOLD,'u_btrop');
  F(F>1e6)=nan;
  Ubt=squeeze(F);

  [F,n,m] = read_hycom(finaOLD,finbOLD,'v_btrop');
  F(F>1e6)=nan;
  Vbt=squeeze(F);

  if adj_btrop>0;
    Ubt=Ubt*adj_btrop;
    Vbt=Vbt*adj_btrop;
  end
  
% Ubrtop and Vbrtrop are used at the boundaries (as well as surf. elev)
% NOTE: Barotrop. U,V are not nan's at the OBs (unlike baroclinic 3D)
%  [Ubt,Vbt] = sub_UVbtrop(UU,VV,dH);  % corrected brtp U,V

% Calculate Montgomery potential from ssh
  montg1 = sub_montg1_regr(finaOLD,finbOLD);
  
%keyboard  
% Creating new files:
% Note in HYCOM TSIS nest files are in 
% format archv.YEAR_DAY_HR.[ab] - see forfun.f
  finaNEW = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthnest,iyr,iday);
  finbNEW = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthnest,iyr,iday);

% Open old files with OB error (anom+mean old nest) for reading
% and files for writing corrected U,V
  disp('Openning files ...');
  fida_OLD = fopen(finaOLD,'r');  % I'm using averaged fields
  fidb_OLD = fopen(finbOLD,'rt'); % use this to copy over *b
  fida_NEW = fopen(finaNEW,'w');
  fidb_NEW = fopen(finbNEW,'wt');

  
  nlevU=0;  % level counter u-vel
  nlevV=0;  % level counter v-vel
  nlevK=0;  % level counter KE
  nlevD=0;  % level counter layer thickness 
  nlevT=0;  % level counter layer T 
  nlevS=0;  % level counter layer S 
  nlevR=0;  % level counter Density anomaly (s.dens - thbase=25)

  if ~exist('TDENS','var')
    TDENS=read_targ_dens(finbOLD,'salin','=',4);
  end
  
% Write header:
%  aa1='Daily nest files: GLBu0.08 anom. + nest climatology GOMl0.04';
%  aa2='Anomalies from: 3hr GLBu0.08 on Z-levels -> 3-day averaged';
  
  for nl=1:10
    aa=fgetl(fidb_OLD);
    if nl==1
      aa=sprintf('TSIS Nest from HYCOM Rnls, Crct Montg,  Ubrtp*%3.2f',adj_btrop);
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
      fprintf('GRID: ID=%i, JD=%i\n',ID,JD);
    end
   
  end

  IJDM=ID*JD;
  npad=4096-mod(IJDM,4096);
  toto=ones(npad,1);

% Surface and barotropic info  
  while ischar(aa)
% Read 1 record from HYCOM nest archv:    
    aa   = fgetl(fidb_OLD);
    if ~ischar(aa), break; end;
    dmm  = fread(fida_OLD,IJDM,'float32','ieee-be');  % read field
    toto = fread(fida_OLD,npad,'float32','ieee-be');  % read npad 


    I=strfind(aa,'=');
    bb=aa(I+1:end);
    S=sscanf(bb,'%f');

%    IM=0;
    IM  = strncmp('montg1',aa,5);
    II  = strncmp('u_btrop',aa,5);
    JJ  = strncmp('v_btrop',aa,5);
    if S(3)==0  | IM  %# hybr. layer = 0 - surface and barotropic field
       aa2=[];
%      fprintf('<= %s\n',aa);
%keyboard
      
      if II,  % write barotropic U
%	keyboard
       fprintf('<= %s\n',aa);
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
        aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);
	aa = aa2;
      elseif JJ
%        lmask = LMASK.vvel;
       fprintf('<= %s\n',aa);
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
        aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);
	aa = aa2;
      elseif IM  % montg potential
       fprintf('<= %s\n',aa);
        lmask = LMASK.temp;
        dmm = reshape(montg1',IJDM,1);
        I=find(isnan(dmm));
        dmm(I)=0;
        dmm(lmask)=hg;
        Inh=find(dmm<hg/10);
        minv=min(dmm(Inh));
        maxv=max(dmm(Inh));
        nlev=S(3);
	td0=S(4);
% Parse string from old nest *b:
        aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);
	aa = aa2;
      end 

% Write in the new nest file:
      fwrite(fida_NEW,dmm,'float32','ieee-be');
      fwrite(fida_NEW,toto,'float32','ieee-be');

      if ~isempty(aa2)
        fprintf('=> %s\n',aa);
      else
        fprintf('<=> %s\n',aa);
      end
      fprintf(fidb_NEW,[aa,'\n']);
      
      continue % read next line from *b / record form *a
      
    end  % if S  - surface and barotropic fields are done
    
    bb=aa(1:I-1);
    wrd = deblank(sscanf(bb,'%s'));

    
    switch (wrd)
% -----------------
%  u-vel
% -----------------
     case('u-vel.NNOT NEEDED')
      nlevU=nlevU+1;
      nlev=nlevU;
      
      lmask = LMASK.uvel;
      A = squeeze(UU(nlev,:,:));
      A = A';
      dmm = reshape(A,IJDM,1);
      I=find(isnan(dmm));
      umm=nanmean(dmm);
      dmm(I)=umm;
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
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);

% Write *b file:
      fprintf('<= %s\n',aa);
      fprintf('=> %s\n',aa2);
      fprintf(fidb_NEW,[aa2,'\n']);

%keyboard
% -----------------
%  v-vel
% -----------------
     case('v-vel.NOT NEEDED')
      nlevV=nlevV+1;
      nlev=nlevV;
      
      lmask=LMASK.vvel;
      A = squeeze(VV(nlev,:,:));
      A = A';
      dmm = reshape(A,IJDM,1);
      I=find(isnan(dmm));
      umm=nanmean(dmm);
      dmm(I)=umm;
      dmm(lmask)=hg;
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
      aa2 = sub_parse_string_B(aa,nlev,td0,minv,maxv);

% Write *b file:
      fprintf('<= %s\n',aa);
      fprintf('=> %s\n',aa2);
      fprintf(fidb_NEW,[aa2,'\n']);

% --------------------------      
     otherwise  % write to new with no changes
      fwrite(fida_NEW,dmm,'float32','ieee-be');
      fwrite(fida_NEW,toto,'float32','ieee-be');
      fprintf('<=> %s\n',aa);
%      fprintf('=> %s\n',aa);
      fprintf(fidb_NEW,[aa,'\n']);
   end;  % switch wrd - field name 	 

    
    
  end;  % while aa - reading *b file all fields

  fclose(fida_NEW);
  fclose(fidb_NEW);
  fclose(fida_OLD);
  fclose(fidb_OLD);
  
  fprintf('  ####    1 record: %9.2f sec\n\n',toc);
  fprintf(' OUTPUT saved: %s\n',finaNEW);
  fprintf(' OUTPUT saved: %s\n\n',finbNEW);


  if ip>1 & iyr~=YRPLT(ip-1,1)
    fprintf('Year %i is done ...\n',iyr-1);
%    fprintf('Cleaning OLD files %s: *.%i*old.[ab] ...\n\n\n',pthtmp,iyr-1);
%    scmd=sprintf('rm %s*.%i*old.[ab]',pthtmp,iyr-1);
%    system(scmd);
  end
  

end;

fprintf('Year %i is done ...\n',iyr);
%fprintf('Cleaning OLD files %s: *.%i*old.[ab] ...\n',pthtmp,iyr);
%scmd=sprintf('rm %s*.%i*old.[ab]',pthtmp,iyr);
%system(scmd);


