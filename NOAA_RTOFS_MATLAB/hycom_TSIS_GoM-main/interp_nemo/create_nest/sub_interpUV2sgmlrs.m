% Interpolate U and V velocity components on z-levels
% into HYCOM hybrid layers 
% from z-level files (U,V are in HYCOM horiz grid, and some z-levels)
%
% In the archive archm file, Layer U is u_total = u_bacl+u_brtp
%  in the archv U = u_bacl
% In Z-level, U=u_total - so need to split u_btrip and u_bacl for archv 
%
% In the Z-level relax-format files, U,V are located in the
% p-point, need to put htem back in the U,V points for archv
% for archm - no need to relocate
%
% 1) interpoalte u_total from nest z-levels 
% into target densities
% 2) relocate to U,V points
% 3) split into barotrop and b/clinic components (for archv)
%
function [Usgm,Vsgm] = sub_interpUV2sgmlrs(fuzlv,fvzlv,Zlev,ZMlev,fina,finb,HH);

fprintf('  ======    Interpolating U,V Z-levels -> HYCOM ...\n');

rg = 9806;
hg=2^100;
onemm=rg/1000; % HYCOM threshold for almost mass-less layer
[JD,ID]=size(HH);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
IDM=ID;
JDM=JD;

% Find layer mid-depths:
dz = diff(abs(Zlev));
nZl = length(Zlev);
nlrs = length(ZMlev);
%ZMlev = zeros(nZl-1,1);
%for kk=1:nZl-1
%  ZMlev(kk,1)=Zlev(kk)-0.5*dz(kk);
%end;
Zlev0=Zlev;
DZ0=dz;

%
% Read Z-level U and V:
fprintf('Reading %s\n',fuzlv);
fid=fopen(fuzlv,'r');
for ik=1:nlrs
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
  F = reshape(dmm,IDM,JDM);
  F = F';
  Uav(ik,:,:) = F;
end
fclose(fid);

fprintf('Reading %s\n',fvzlv);
fid=fopen(fvzlv,'r');
for ik=1:nlrs
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
  F = reshape(dmm,IDM,JDM);
  F = F';
  Vav(ik,:,:) = F;
end
fclose(fid);



% Read layers:
[dP,n,m,l] = read_hycom(fina,finb,'thknss');
dP(dP>=0.1*hg) = nan;  % Pressure, Pa

% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
dH=dP/rg;
[ZM,ZZ] = sub_zz_zm(fina, finb,HH,'thknss',dH);
%keyboard
%  -------------------------------------------------- 
%  For each interface in new t.dens. find
%  nearest upper and bottom interfaces in old t.dens.
%
%   Old v. layers             New v. layers
%   -------------
%       Z(kold) *           ----------   
%                            Zn(knew)   *
%   -------------
%
%     Z(kold+1) *           ----------   
%                            Zn(knew+1) *
%   -------------
%  Velocity is calculated as follows:
% Unew(k) = (dHold(k)*Uold(K)+ dHold(k+1)*Uold(k+1)+ ...)/thickness of 
%                                                         hybrid sgm layer  
% this approach conserves mass flux 
% Three cases are possible:
% (1) the upper interf. of new targ. layer is within the old layer K, 
%    the bottom interf. of new t. dens is
% within the old layer k+1
%
% (2) whole new t. layer is within one old t. layer (in this case U new = U old)
% 
% (3) new. t. layer extends over several old t. layers
% --------------------------------------------------- 

IOcean=find(HH<0);
Usgm = zeros(l,m,n)*nan;
Vsgm = zeros(l,m,n)*nan;
%keyboard
for io=1:length(IOcean);
%  fprintf('io=%i\n',io);
  if mod(io,5000)==0;
    fprintf(' === U,V z -> Sigma2:  Done %8.2f%%\n',io/length(IOcean)*100);
  end
  
  
  II=IOcean(io);
  [jj,ii]=ind2sub(size(HH),II);
  Zsgm = squeeze(ZZ(:,jj,ii));
  Uz0  = squeeze(Uav(:,jj,ii));
  Vz0  = squeeze(Vav(:,jj,ii));
  nS   = length(Zsgm);
  
% NEMO z-interpolated fields are only to ~5900m
% stretch it to the bottom in deep ocean
  hb=HH(jj,ii);
  Zlev=Zlev0;
  if abs(hb)>abs(Zlev0(end));
    Zlev(end)=hb-0.01;
  end;
  
  
% Check for spurious U,V:
%  IS=find(abs(Uz0>2));
%  JS=find(abs(Vz0>2));
%  if ~isempty(IS)
%    dmm=Uz0;
%    dmm(IS)=nan;
%    Uz0(IS)=sign(Uz0(IS)).*nanmean(dmm);
%  end
%  if ~isempty(JS)
%    dmm=Vz0;
%    dmm(JS)=nan;
%    Vz0(JS)=sign(Vz0(JS)).*nanmean(dmm);
%  end
  
  Isf=min(find(Uz0==0)); % land
  if Isf==1;  % U,V point on land
    continue;
  end
% Interpolate into grid mid-depths:
% Given U,V are at mid-grid depths - no need to interp to ZM
% Fill Nan's with something  
% Ib = min(find(isnan(Uz0)));
% if ~isempty(Ib)
%   Uz0(Ib:end)=9.99e-5;  % fill NaN bottom values
%   Vz0(Ib:end)=9.99e-5;  % fill NaN bottom values
% end
% 
% Uz = interp1(Zlev,Uz0,ZMlev,'pchip');
%  Uz = interp1(Zlev,Uz0,ZMlev,'linear');
% Vz = interp1(Zlev,Vz0,ZMlev,'pchip');
  Uz = Uz0;
  Vz = Vz0;
  
% U,V velocity components in HYCOM hybrid layers (sigma0 dens.)  
  Us = [];
  Vs = [];
  for k=1:nS-1
    dHs=abs(Zsgm(k+1)-Zsgm(k)); % thicknes, sigma hybrid layer
    s1=Zsgm(k);
    s2=Zsgm(k+1);
    
% Find corresponding Z surfaces
% which are top surfaces of the 
% z-layers, that contain HYCOM surfaces
    k1=max(find(Zlev>=s1));  % top HYCOM surface
    k2=max(find(Zlev>s2)); % btm HYCOM surface
    
    z1=Zlev(k1);
    z2=Zlev(k2);
    
% How many complete z-levels inside the hybrid layer:
    nZ = k2-k1-1;
    
    chck=0;
    if chck==1
      figure(10); clf;
%      xx=ones(nZ,1);
      hold on;
      for kz=1:nZl
        plot([0.98 1.02],[Zlev(kz) Zlev(kz)],'b--');
      end
      xx=ones(nS,1);
      for kz=1:nS
        plot([0.98 1.02],[Zsgm(kz) Zsgm(kz)],'r--');
      end
      plot([0.9 1.1],[z1 z1],'b-','linewidth',1.5);
      plot([0.9 1.1],[z2 z2],'b-','linewidth',1.5);
      plot([0.9 1.1],[s1 s1],'m--','linewidth',1.5); % HYCOM sigma2 interface
      plot([0.9 1.1],[s2 s2],'m--','linewidth',1.5);
    end
%
% This is simplified version of U,V interpolation
% More accurate approach: for each grid cell u=U(z) - a linear
% function between the interfaces
% Then the transport between several grid cells (where bottom and top Z cells 
% are only partial - where HYCOM interfaces are) should be
% intgral(zb:Zinterface(k1)) of U(z) + Integral of full cells of U(z) + 
% partial integral(Zinterface(k2):ztop) of U(z)
%    if io==113981; keyboard; end;
    if nZ<0 % HYCOM layer is within a Z-layer    
%      dZtot = abs(z2-z1);
%      dStot = abs(s2-s1);
%      Uztrp = Uz(k1)*dZtot;
%      Vztrp = Vz(k1)*dZtot;
      Us(k,1) = Uz(k1);
      Vs(k,1) = Vz(k1);
    else % upper isop. surf is in layer k, bottm - layer k+1
% Delta thickness of partial layers      
      dZS1=abs(s1-Zlev(k1+1));
      dZS2=abs(s2-Zlev(k2));
% Full layer thickness:
%      dZ1=abs(Zlev(k1)-Zlev(k1+1));
%      dZ2=abs(Zlev(k2-1)-Zlevk(k2));
% Partial transport
% Use mean vel. within z-layers, i.e. interpolated
% into middle of the grid cell
%     ?? if k>1  this seems wrong
%        pU1=Uz(k1)*dZS1; 
%        pV1=Vz(k1)*dZS1;
%    ??  else % surface - no partial transport above surface
%	pU1=0;
%	pV1=0;
%     ?? end`
%keyboard
      pU1=Uz(k1)*dZS1; 
      pV1=Vz(k1)*dZS1;
      pU2=Uz(k2)*dZS2;
      pV2=Vz(k2)*dZS2;
      
      dZfull=0;
      Ufull=0;
      Vfull=0;
      cc=0;
      for kf=k1+1:k2-1
	cc=cc+1;
	dZ1=abs(Zlev(kf)-Zlev(kf+1));
	dZfull=dZfull+dZ1;
	Ufull=Ufull+Uz(kf)*dZ1;
	Vfull=Vfull+Vz(kf)*dZ1;
      end
% Check total depth:
      dZtot=dZS1+dZfull+dZS2;
      err=abs(abs(dZtot/(s2-s1))-1);
      if err>1e-2
	error('*** Check total depths over %i nZ = %8.2m\n',nZ,dZtot);
      end
% Total transport in HYCOM sigma layer:
      UTs=pU1+Ufull+pU2;
      VTs=pV1+Vfull+pV2;
      Us(k,1)=UTs./dZtot;
      Vs(k,1)=VTs./dZtot;
      
    end;  % if 
%      keyboard
  end;    % v. layers
  
  Usgm(:,jj,ii)=Us;
  Vsgm(:,jj,ii)=Vs;
  
%  if hb<-3000; keyboard; end;

% Check transport:
%  Izr = min(find(Uz == 0));
%  if ~isempty(Izr),
%    ibtm=Izr-1;
%  else
%    ibtm=nZl-1;
%  end
  dHz = -1*abs(cumsum(DZ0));
  ibtm = min(find(dHz<=hb));
  if isempty(ibtm), ibtm=nZl-1; end;
  
  dH=squeeze(dP(:,jj,ii))/rg;
  hb0=sum(dH);
  if abs(1-abs(hb0/hb))>0.1,
    fprintf('Check depths in dP and Z arrays: bottom %8.2f  %8.2f\n',...
	    abs(hb0),abs(hb));
    keyboard;
  end
  
  UZ1=sum(Uz(1:ibtm-1).*abs(DZ0(1:ibtm-1)));
  dUz=Uz(ibtm)*abs(sum(DZ0(1:ibtm-1))-abs(hb)); % near-bottom transp:
  UZtot=UZ1+dUz;
  UStot=sum(Us.*abs(dH));
  RUtrt=abs(1-UStot/UZtot);

  VZ1=sum(Vz(1:ibtm-1).*abs(DZ0(1:ibtm-1)));
  dVz=Vz(ibtm)*abs(sum(DZ0(1:ibtm-1))-abs(hb)); % near-bottom transp:
  VZtot=VZ1+dVz;
  VStot=sum(Vs.*abs(dH));
  RVtrt=abs(1-VStot/VZtot);
  
  if (RUtrt>1 & abs(UZtot>0.01) & abs(UStot>0.01)) | ...
     (RVtrt>1 & abs(VZtot>0.01) & abs(VStot>0.01))
%    error('Check tranport not conserved');
    fprintf('CHECK: sub_itnerpUV2sgm: RUtrt=%5.1f RVtrt=%5.1f\n',...
	    RUtrt,RVtrt);
    fprintf('sub_interpUV2sgm:   Check tranport not conserved\n');
    keyboard
  end
  
  
end;  % ocean points;

%keyboard

f_chck=0;
if f_chck>0
%  io=8000;
  II=IOcean(io);
  [jj,ii]=ind2sub(size(HH),II);
%  jj = m-1;
%  ii = 465;
  Zsgm = squeeze(ZZ(:,jj,ii));
  Uz0  = squeeze(Uav(:,jj,ii)); % U on Z
  Vz0  = squeeze(Vav(:,jj,ii)); % V on Z

  figure(10); clf;
  
  axes('Position',[0.08 0.06 0.35 0.85]);
  hold on;
%  plot(Uz0,ZMlev,'g-'); % Original U on z, cell interfaces
  % Plot U in Z-level interpolated
  % into mid-depth grid layers
  for k=1:nZl-1
    u=Uz(k);
    plot([u u],[Zlev0(k) Zlev0(k+1)],'r-');
    if k>1
      um1=Uz(k-1);
      plot([um1 u],[Zlev0(k) Zlev0(k)],'r-');
    end;
  end
  % Plot U interpolated into HYCOM
  for k=1:nS-1
    u=Usgm(k,jj,ii);
    plot([u u],[Zsgm(k) Zsgm(k+1)],'b-');
    if k>1; 
      um1=Usgm(k-1,jj,ii); 
      plot([um1 u],[Zsgm(k) Zsgm(k)],'b-');
    end;
  end
  hb=HH(jj,ii);
  plot([min(Uz) max(Uz)],[hb hb],'k--');
  title('U: Zgrid red, HYCOM blue' );
  set(gca,'ylim',[1.05*hb 0])
  
  axes('Position',[0.52 0.06 0.35 0.85]);
  hold on
%  plot(Vz0,ZMlev,'g-');
  for k=1:nZl-1
    v=Vz(k);
    plot([v v],[Zlev0(k) Zlev0(k+1)],'r-');
    if k>1
      vm1=Vz(k-1);
      plot([vm1 v],[Zlev0(k) Zlev0(k)],'r-');
    end;
  end
  for k=1:nS-1
    v=Vsgm(k,jj,ii);
    plot([v v],[Zsgm(k) Zsgm(k+1)],'b-');
    if k>1; 
      vm1=Vsgm(k-1,jj,ii);
      plot([vm1 v],[Zsgm(k) Zsgm(k)],'b-');
    end;
  end
  hb=HH(jj,ii);
  plot([min(Vz) max(Vz)],[hb hb],'k--');
  title('V, Zgrid: red, HYCOM: blue');
  set(gca,'ylim',[1.05*hb 0])

  btx = 'sub_interpUV2sgm2.m';
  bottom_text(btx,'pwd',1);

keyboard  
end

    

return
