% interpolate coarser GLORYS into HYCOM
% Interpolate GLORYS into HYCOM-TSIS
% Outside NEMO domain
% Use bilinear interpolation:
% F= a0+a1*x+a2*y+a3*x*y - basis functions [a0,a1,a2,a3]
% are obtained in calculate_glorys_weights.m
%
function Ai = sub_interp_glorys2hycom(LONN,LATN,LON,LAT,HH,GLR,AA,zz0)
fprintf('Interpolating GLORYS ---> HYCOM ...\n');

if isempty(zz0)
  zz0=0;
end

Ai = HH*nan;
Ihycom = GLR.IndxHYCOM;

tic;
nI=length(Ihycom);
for ii=1:nI
  if mod(ii,80000)==0,
    fprintf(' %6.2f%% processed, %9.6f min ...\n',ii/nI*100,toc/60);
  end

  i0=Ihycom(ii);
%  if HH(i0)>zz0, continue; end; % bottom/land
% Quadr. Points:
  Jglr = GLR.IGLRS(ii,:);
  Phi  = squeeze(GLR.PHIGL(ii,:,:));
  F    = AA(Jglr);
  xx0  = LON(i0);  % HYCOM lon
  yy0  = LAT(i0);  % HYCOM lat
  XY = [1; xx0; yy0; xx0*yy0];
  Fi = F*Phi*XY;

  Ai(i0) = Fi;

  fchck=0;
  if fchck==1
    [jh,ih] = ind2sub(size(HH),i0);

    figure(1); clf;
    hold on;
    plot(LON(i0),LAT(i0),'r*');
    plot(LONN(Jglr),LATN(Jglr),'b.');
    plot(LON(jh-1,ih),LAT(jh-1,ih),'k+');
    plot(LON(jh+1,ih),LAT(jh+1,ih),'k+');
    plot(LON(jh,ih+1),LAT(jh,ih+1),'k+');
    plot(LON(jh,ih-1),LAT(jh,ih-1),'k+');
keyboard
  end

%  keyboard
end

return
