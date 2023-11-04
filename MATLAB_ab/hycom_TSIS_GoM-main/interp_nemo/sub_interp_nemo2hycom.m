% High resolution NEMO into coarser resolution HYCOM-TSIS:
% Interpolate by "binning" - average surrounding NEMO points
% into center of HYCOM
% AA - NEMO field
% Ai - HYCOM interpolated
%
function Ai = sub_interp_nemo2hycom(LONN,LATN,LON,LAT,HH,DX,DY,NMI,AA);

fprintf('Interpolating NEMO --> HYCOM ...\n');
[ah1,ah2]=size(HH);
[an1,an2]=size(LONN);
ln1=min(min(LONN));
ln2=max(max(LONN));
lt1=min(min(LATN));
lt2=max(max(LATN));


Ai = HH*nan;

Ihycom = NMI.IndxHYCOM;


tic;
nI=length(Ihycom);
for ii=1:nI
  if mod(ii,20000)==0,
    fprintf(' %6.2f%% processed, %9.6f min ...\n',ii/nI*100,toc/60);
  end

  i0=Ihycom(ii);
  dmm=NMI.INEMO(ii,:);
  iok = find(dmm>0);
  Jnemo=dmm(iok);
  Ai(i0)=nanmean(AA(Jnemo));

  fchck=0;
  if fchck==1
    [jh,ih] = ind2sub(size(HH),i0);

    figure(1); clf;
    hold on;
    plot(LON(i0),LAT(i0),'r*');
    plot(LONN(Jnemo),LATN(Jnemo),'b.');
    plot(LON(jh-1,ih),LAT(jh-1,ih),'k+');
    plot(LON(jh+1,ih),LAT(jh+1,ih),'k+');
    plot(LON(jh,ih+1),LAT(jh,ih+1),'k+');
    plot(LON(jh,ih-1),LAT(jh,ih-1),'k+');

  end

%  keyboard
end

return

