% Interpolate GLORYS into HYCOM-TSIS
% 3D field
% 2D at a time:
% interpolate onto HYCOM grid and interpolate onto a depth level z0
% thus need 2 glorys depth layers
%
% Use bilinear interpolation:
% F= a0+a1*x+a2*y+a3*x*y - basis functions [a0,a1,a2,a3]
% are obtained in calculate_glorys_weights.m
function [AAi, Amsk] = sub_glorys2tsis_3d(dnmb,pthnemo,pthdata,pthglorys,LAT,LON,HH,...
                                   fldnm,zz0,Amsk);

[mm,nn] = size(HH);
[XM,YM]=meshgrid([1:nn],[1:mm]);


% Get GLORYS field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

LONN = nc_varget(flglrs,'longitude');
LATN = nc_varget(flglrs,'latitude');
[LNN,LTN] = meshgrid(LONN,LATN);

ZZN  = nc_varget(flglrs,'depth');
ZZN = -ZZN;

iz1 = max(find(ZZN>=zz0));
if isempty(iz1); error('Couldnot locate depth %6.2f',zz0); end;
iz2 = iz1+1;

% Bottom layer in GLORYS is shallower than 
% in NEMO
if zz0<ZZN(end)
  iz2=iz1;
end

%keyboard
% Sanity check
if ZZN(iz1)<zz0 | (ZZN(iz2)>zz0 & iz2<length(ZZN))
  fprintf('sub_glorys2tsis_3d: Check GLORYS depths \n');
  keyboard;
end


%keyboard

switch(fldnm),
 case('temp');
  fldgl='thetao';
 case('saln');
  fldgl='so';
 case('uvel')
  fldgl='uo';
 case('vvel')
  fldgl='vo';
end

nlon = length(LONN);
nlat = length(LATN);
ndpth= length(ZZN);
zz1  = ZZN(iz1);
zz2  = ZZN(iz2);


AA1 = squeeze(nc_varget(flglrs,fldgl,[0 iz1-1 0 0],[1 1 nlat nlon]));
A01 = AA1;

%switch(fldnm)
% case({'temp','saln'})
%		if isempty(Amsk); 
%				AA1 = sub_fill_land(AA1);
%				Amsk=AA1; 
%		else
%				I=find(isnan(AA1));
%				AA1(I)=Amsk(I);
%				Amsk=AA1;
%		end;
% case({'uvel','vvel'})
%  I=find(isnan(AA1));
%  AA1(I)=0;
%end
%switch(fldnm)
% case({'temp','saln'})
  if isempty(Amsk);
    AA1 = sub_fill_land(AA1);
    Amsk=AA1;
  else
    I=find(isnan(AA1));
    AA1 = sub_fill_land(AA1);
    amx = max(max(AA1));
    if amx<-100
      AA1=Amsk;
    else
      AA1(I)=0.5*(Amsk(I)+AA1(I));
    end
    Amsk=AA1;
  end;
% case({'uvel','vvel'})
%  I=find(isnan(AA));
%  AA(I)=0;
%end


if iz2~=iz1
		AA2 = squeeze(nc_varget(flglrs,fldgl,[0 iz2-1 0 0],[1 1 nlat nlon]));
		A02 = AA2;

%  switch(fldnm)
%   case({'temp','saln'})
%				I=find(isnan(AA2));
%				AA2(I)=Amsk(I);
%				Amsk = AA2;
%   case({'uvel','vvel'})
%    I=find(isnan(AA2));
%    AA2(I)=0;
%  end
  if isempty(Amsk);
    AA2 = sub_fill_land(AA2);
    Amsk=AA2;
  else
    I=find(isnan(AA2));
    AA2 = sub_fill_land(AA2);
    amx = max(max(AA2));
    if amx<-100
      AA2=Amsk;
    else
      AA2(I)=0.5*(Amsk(I)+AA2(I));
    end
    Amsk=AA2;
  end;


%		AA2 = sub_fill_land(AA2); % horizontal fill
  AA = sub_vrtintrp_glorys(AA1,AA2,zz1,zz2,zz0);
else
  AA = AA1;
end

fprintf('GLORYS %s interp to %6.2f, min/max: %7.4f/%7.4f\n',fldnm,zz0,...
        min(min(AA)),max(max(AA)));


% Horizontal interpolation:
mm = size(LATN,1);
nn = size(LONN,1);
% 
% Interpolation points, polynomial coeff. 
fglrint = sprintf('%sGLRS_TO_HYCOM_TSIS_interp_pnts.mat',pthdata);
fprintf('Loading %s\n',fglrint);
GLR = load(fglrint);

AAi = sub_interp_glorys2hycom(LNN,LTN,LON,LAT,HH,GLR,AA,zz0);
%AAi(HH>zz0)=nan;


% ------------------------
% Check & plot:
%
f_chck=0;
if f_chck==1;
% Get NEMO grid:
  fprintf('Checking fields, plotting ...\n');
		f_get_grid=0;
		fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
		if f_get_grid==1
				[LONN,LATN,ZZN] = sub_get_NEMO_grid(dnmb);
				LONN = double(LONN);
				LATN = double(LATN);
				ZZN = double(ZZN);
				fprintf('Saving grid %s\n',fgrd);
				save(fgrd,'LONN','LATN','ZZN');
		else
				fprintf('Loading NEMO grid %s\n',fgrd);
				load(fgrd);
		end


%
% Plot Original GLORYS and interpolated onto HYCOM-TSIS IAS
% HYCOM fields are also vertically interpolated
  nf=1;
  stl1 = sprintf('%s z=%6.2f GLORYS',fldnm,zz1);
  stl2 = sprintf('z=%6.2f IAS HYCOM-TSIS',zz0);
  sub_plot2flds(nf,LNN,LTN,A01,LON,LAT,AAi,LONN,LATN,stl1,stl2,fldnm,zz0);
%  sub_plot2flds(nf,LNN,LTN,AA1,LON,LAT,AAi,LONN,LATN,stl1,stl2,fldnm,iz2);
 
  btx = 'sub_glorys2tsis_3d.m';
  bottom_text(btx,'pwd',1);


  nf=2;
  stl1 = sprintf('%s z=%6.2f GLORYS',fldnm,zz2);
  sub_plot2flds(nf,LNN,LTN,A02,LON,LAT,AAi,LONN,LATN,stl1,stl2,fldnm,zz0);

  bottom_text(btx,'pwd',1);

keyboard 

end

return

