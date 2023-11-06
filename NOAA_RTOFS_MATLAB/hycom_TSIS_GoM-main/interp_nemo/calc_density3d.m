% Compute potential density
% Needed for creating relax climatology files
% for vertical interpolation onto HYCOM
% hybrid layers
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



% 
% Read 
foutT = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
													pthoutp,'temp',DV(3),DV(2),DV(1));
fprintf('Reading: %s\n',fout);
fidT=fopen(foutT,'r');
foutS = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
													pthoutp,'saln',DV(3),DV(2),DV(1));
fidS=fopen(foutS,'r');

for ik=1:nlrs
  fprintf('Reading %s layer %i\n',foutT,ik);
  dmm=fread(fidT,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fidT,npad,'float32','ieee-be');  % read npad   
  T = reshape(dmm,IDM,JDM);
  T = T';
%
  fprintf('Reading %s layer %i\n',foutS,ik);
  dmm=fread(fidS,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fidS,npad,'float32','ieee-be');  % read npad   
  S = reshape(dmm,IDM,JDM);
  S = S';

  Rho = sw_dens0(S,T);
keyboard


		if f_write==1
		% Fill land
				F = sub_fill_land(Ahycom);
				F = F';
				F = reshape(F,IJDM,1);

				if isempty(fid)
						fout = sprintf('%s%s_nemoglorys2hycom_%2.2i%2.2i%4.4i.dat',...
																			pthoutp,'dens',DV(3),DV(2),DV(1));
						fid  = fopen(fout,'w','ieee-be');
				end

				fwrite(fid,F,'float32','ieee-be');
				fwrite(fid,toto,'float32','ieee-be');
				fprintf('Written files: %s\n',fout);

		end



end
fclose(fidT);
fclose(fidS);







