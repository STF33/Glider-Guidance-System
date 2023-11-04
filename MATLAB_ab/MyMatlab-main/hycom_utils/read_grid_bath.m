function GRD = read_grid_bath(flbase,varargin);
%
% function GRD = read_grid_bath(flbase,varargin);
% reads grid info from flbase.b and flbase.a files
% if topo file (*.nc or *.a) is specified, 
% then topo is also read 
%
% e.g., read_grid_bath(hycom_topo/regional.grid,...
%             hycom_topo/depth_A08_T11.a);
%
% Returns GRD.File_grid = flbase;
%         GRD.File_topo = flda;
%         GRD.PLON = plon;
%         GRD.PLAT = plat;
%         GRD.IDM = IDM;
%         GRD.JDM = JDM;
%         GRD.Topo = HH;
%

flgb=[flbase,'.b'];
flga=[flbase,'.a'];
fprintf('Reading hycom grid from %s\n',flga);
fprintf('Reading hycom grid from %s\n',flgb);
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

fprintf('Grid info: IDM = %i, JDM = %i\n',IDM,JDM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read lon/lat from regional grid file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_fid=fopen(flga,'r');

plon=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,4*(npad+IJDM),-1);
plat=fread(grid_fid,IJDM,'float32','ieee-be');

plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(grid_fid);

% --------------------------
I=find(plon>180);
plon(I)=plon(I)-360;
I=find(plon<-180);
plon(I)=plon(I)+360;

if nargin>0
  flda=varargin{1};

% Define extension of the depth file:
  nll = length(flda);
  bnd = flda(nll-1:end);
  if strcmp(bnd,'nc');
    ext = 'nc';
  elseif strcmp(bnd,'mat');
    ext = 'mat';
  else
    ext = 'a';
  end
  
%  ext = deblank(flda(idt+1:end));
  switch (ext)
    case('a')
 
      fprintf('Reading HYCOM TOPO  from %s\n',flda);
      depth_fid=fopen(flda,'r');
      fprintf('Reading topo from %s\n',flda);

      fseek(depth_fid,6*4*(npad+IJDM),-1);
      [h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
%      y=find(h>1e10);
%      h(y)=nan;
      HH=reshape(h,IDM,JDM)';
      clear h

      fclose(depth_fid);

% Make depths < 0:
      I = find(HH > 1.e20);
      HH = -HH;
      HH(I) = 100.;
    case('nc');
      fprintf('Reading HYCOM TOPO  from %s\n',flda);
      HH=nc_varget(flda,'Bathymetry');

    case('mat');
      fprintf('Reading HYCOM TOPO  from %s\n',flda);
      HH=load(flda);
  end
%

else
  HH=[];
end


GRD.File_grid = flbase;
GRD.File_topo = flda;
GRD.PLON = plon;
GRD.PLAT = plat;
GRD.IDM = IDM;
GRD.JDM = JDM;
GRD.Topo = HH;

return


