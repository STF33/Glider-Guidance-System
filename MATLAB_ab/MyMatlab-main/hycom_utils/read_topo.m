function [HH,plon,plat] = read_topo(flda,flga,flgb);
% function [HH,plon,plat] = read_topo(flda,flga,glfb);
% Read binary topo file
% from binary file flda
% Need flgb - to get info about region
% if flga not empty will also read LON/LAT

% Topo/grid files, input:
fprintf('Reading topo %s\n',flda);
% Read in topo from *a file:
f1=fopen(flgb);  % read I,J from regional.grid.b
aa=fgetl(f1);
dmm=aa(2:8);
IDM=str2num(dmm);
aa=fgetl(f1);
dmm=aa(2:8);
JDM=str2num(dmm);
IJDM=IDM*JDM;
fclose(f1);

npad=4096-mod(IJDM,4096);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read lon/lat from regional grid file
if isempty(flga),
  plon=[];
  plat=[];
  return
end

grid_fid=fopen(flga,'r');

[plon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
stat=fseek(grid_fid,4*(npad+IJDM),-1);
if stat<0
  error('Reading grid file ...');
end
[plat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

disp('Reading lat/lon  ...')
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(grid_fid);

% --------------------------
I=find(plon>180);
plon(I)=plon(I)-360;
I=find(plon<-180);
plon(I)=plon(I)+360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read bathymetry from regional.depth.a

depth_fid=fopen(flda,'r');

%fseek(depth_fid,6*4*(npad+IJDM),-1) % this is confusing, should not be needed
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
%y=find(h>1e10);
%h(y)=nan;
HH=reshape(h,IDM,JDM)';

fclose(depth_fid);

I=find(HH>1e20);
HH=-HH;
HH(I)=100;

return
