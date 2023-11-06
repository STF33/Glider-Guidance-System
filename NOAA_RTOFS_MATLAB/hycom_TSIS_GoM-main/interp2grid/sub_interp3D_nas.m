function Tzx=sub_interp3D_nas(T,IFL,XX,YY,HHs,ZMnas,ZM,nlon,nlat,ILAND,pfld);
% Interpolate 3D field (T)
% onto NAS grid
% 1st step: fill nans in 3D
% 2nd: interpolate into z levels
% 3rd: interpolate onto horiz grid for all levels
% 4th: put land/bottom masks
huge=1e30;

[msb,nsb]=size(XX);
knas=length(ZMnas);
lnas=length(ZMnas);
[mnas,nnas]=size(nlon);

% fill nans in the surface layer  
t1=squeeze(T(1,:,:));
[t1,IFL]=sub_fill_land(t1,IFL);
T(1,:,:)=t1;

% Fill nans in the subsurface layers (bottom):
T = sub_fill_3z(T);

% Interpolate into z levels
% as a 2D vert sections
Tzi=[];
for ii=1:nsb
  if mod(ii,200)==0
    fprintf('%s interp to z, done %5.2f%% ...\n',pfld,ii/nsb*100);
  end

  y1=YY(:,ii);
  ZMh=squeeze(ZM(:,:,ii));
  F=squeeze(T(:,:,ii));
  Hs = HHs(:,ii);
  Hs(Hs>=0)=-1;

  Fi=sub_interp2z_2D(Hs,F,ZMnas,ZMh);

  f_sct=0;
  if f_sct==1
    sub_check_sect(F,Fi,ZMh,ZMnas,Hs);
  end

  Tzi(:,:,ii)=Fi;
end

% Interpolate onto NAS grid
fprintf('%s Interpolating onto XY grid\n',pfld);
[a1,a2,a3]=size(Tzi);
Tzx=[];
for kk=1:a1
  F=squeeze(Tzi(kk,:,:));
  Fi=interp2(XX,YY,F,nlon,nlat);
  Tzx(kk,:,:)=Fi;
end

% Put land back
for kk=1:knas
  il=ILAND(kk).I;
  Tzx(kk,il)=huge;
end
%




return