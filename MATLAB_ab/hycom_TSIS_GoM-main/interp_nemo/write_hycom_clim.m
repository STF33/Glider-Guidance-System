% Write interpolated NEMO-GLORYS --> IAS HYCOM 
% in climatology files on Z-levels
% see interp3D_nemo_glorys_hycom.m 
% to interpolate data
% 
% HYCOM climatology files are needed to create relax files
% with HYCOM hybrid layers
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


fid = [];

f_write = 0;
fldnm = 'temp';
%fldnm = 'saln';
%fldnm = 'uvel';
%fldnm = 'vvel';

sgm = 'sig2';  % pressure reference
nlrs = 75;     % depth layers in NEMO

dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

pthnemo   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_tmp/';
pthtopo   = '/Net/kronos/ddmitry/hycom/TSIS/';
pthglorys = '/nexsan/people/abozec/TSIS/data/GLORYS/';
pthdata   = '/Net/kronos/ddmitry/hycom/TSIS/datamat/';
pthoutp   = '/Net/kronos/ddmitry/hycom/TSIS/nemo_glorys_interp/';

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

% Load NEMO grid:
% Get NEMO grid and depth levels
f_get_grid=0;
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
fprintf('Loading NEMO grid %s\n',fgrd);
load(fgrd);


% 
% Read 
fin = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
													pthoutp,fldnm,DV(3),DV(2),DV(1));
fprintf('Reading: %s\n',fin);
fid=fopen(fin,'r');

fouta = sprintf('%s%s_zlev_hycom_%s_%2.2i%2.2i%4.4i.a',...
             pthoutp,fldnm,sgm,DV(3),DV(2),DV(1));
foutb = sprintf('%s%s_zlev_hycom_%s_%2.2i%2.2i%4.4i.b',...
             pthoutp,fldnm,sgm,DV(3),DV(2),DV(1));

fida = fopen(fouta,'w');
fidb = fopen(foutb,'wt');


% Write headers in *b file
stl = 'NEMO+GLORYS on IAS HYCOM-TSIS\n';
fprintf(fidb,stl);
switch(fldnm),
 case('temp');
  stl='Potential Temperature';
  bline = 'potential temperature: depth,range =';
 case('saln');
  stl='Salinity';
  bline = '             salinity: depth,range =';
end
fprintf(fidb,stl);
stl=' \n';
fprintf(fidb,stl);
fprintf(fidb,stl);
fprintf(fidb,stl);

stl=sprintf('i/jdm = %i %i\n',IDM,JDM);
fprintf(fidb,stl);

for ik=1:nlrs
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
  F = dmm;

  fmin = min(F);
  fmax = max(F);
  fwrite(fida,F,'float32','ieee-be');
  fwrite(fida,toto,'float32','ieee-be');  

% *b file:
  zz0=ZZN(ik);
  if (ik==1); zz0=0; end;
  aa2 = sub_parse_cline(bline,zz0,fmin,fmax);
  fprintf('===>  %s\n',aa2);
  fprintf(fidb,[aa2,'\n']);

end
fclose(fid);
fclose(fida);
fclose(fidb);










