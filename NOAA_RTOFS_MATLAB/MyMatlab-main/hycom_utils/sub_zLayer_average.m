function Sav = sub_zLayer_average(HH,ZZ,F,zbtm,zz1,zz2);
% Average scalar field
% from HYCOM layers over
% specified fixed-depth lelvels: zz1 and zz2
% 
%  --------------------  ZZtop(kk-1) HYCOM top intrf
%
%
%
%  .......................  zz1 = -300 m 
%   dZ(1)
%  --------------------  ZZbtm(kk-1)=ZZtop(kk)  HYCOM btm intrf
%
%    dZ(2)
%
%  --------------------  ZZbtm(kk)=ZZtop(kk+1)  HYCOM btm intrf
%
%  dZ(3)
%  .......................  zz2 = -500 m 
%
%  --------------------  ZZbtm(kk+1)=ZZtop(kk+2)  HYCOM btm intrf
%
%
%  hZ=sum(dZ(i))
%
% Also possible when zz1 and zz2 are within one HYCOM Layer!
%
% Input: HH - HYCOM topo
% zz1,zz2 - fixed depth layers within which field is averaged
% ZZ - 3D array of HYCOM interface depths
% F  - scalar field
% zbtm - threshold to designate the whole-depth 
%         averaging if zz2 is deeper, note sign! (both <0)
%
Sav = [];
fprintf('::: Layer averging %6.2f - %6.2f\n',zz1,zz2);
[llz,mm,nn]=size(ZZ);
ll=llz-1;
% 
% Find the layers that include the specified
% depth interval
% kZt - model layer with top intrf above or at zz1
%                   and btm intrf below zz1
% kZb - model layer with top intrf above zz2
%                   and btm intrf at or bloew zz2
%
%
%  --------------------  ZZ(kZt)=ZZtop HYCOM top intrf
%
%
%
%  .......................  zz1 = -300 m 
%   dZ(1)
%  --------------------  ZZ(kZt+1)==ZZbtm  HYCOM btm intrf
%
kZt  = HH*0;
kZb  = HH*0;
%cZk  = HH*0; % count v. layers
for kk=1:ll
  ZZtop = squeeze(ZZ(kk,:,:));
  ZZbtm = squeeze(ZZ(kk+1,:,:));
%
% Whole depth:
  if zz2<zbtm % bottom, whole water column layer
    I=find(~isnan(ZZbtm)); % find bottom
    kZb(I)=kk;
    kZt(I)=1; % surface
    continue
  end

% Layers:  
  I=find(HH<zz2 & ZZbtm<zz1);
  if isempty(I), continue; end % too shallow

% Find end of integration 
% for points that are within the zz1-zz2 for given kk
  I=find(HH<zz2 & ZZtop>zz2 & kZt>0);
  if ~isempty(I),   
    kZb(I)=kk;
%    cZk(I)=cZk(I)+1;
  end
  
% Find start for points that haven't started yet 
  I=find(HH<zz2 & ZZtop>=zz1 & ZZbtm<zz1 & kZt==0);
  if ~isempty(I),   
    kZb(I)=kk;
    kZt(I)=kk; % keep 1st layer
%    cZk(I)=cZk(I)+1;
  end
% Too deep:
  if max(max(ZZtop))<zz2, break; end;
end
%keyboard

kZb(kZb==0)=nan; % bottom layer within the depth interval
%      kZt=kZb-cZk+1;   % top layer
%kZt=kZb-cZk;   % need to have the top layer - above zz1
kZt(kZt==0)=nan;  % surface layer
Sav = HH*0; % S mean
lmx = max(max(kZb));
lmn = min(min(kZt));
hZ  = HH*0;
smm = HH*0;
for kk=lmn:lmx
  I=find(kZb>=kk & kZt<=kk);
  if isempty(I), continue; end;
%  SS = HH*nan;
%  SS(I) = squeeze(F(kk,I)); % S at level kk, only pnts that fall within zz1 zz2
  SS    = squeeze(F(kk,:,:)); % 
  ZZtop = squeeze(ZZ(kk,:,:)); % top intrf of model layer
  ZZbtm = squeeze(ZZ(kk+1,:,:)); % btm intrf model layer

% Add dZ from incomplete top layer if needed
% only for HYCOM layers that are partally 
% within the Z-layer
  Iz=find(ZZtop>zz1 & ZZbtm>zz2 & kZb>=kk & kZt<=kk);
  if ~isempty(Iz);
    dZ = abs(ZZbtm-zz1);
    smm(Iz)=smm(Iz)+SS(Iz).*dZ(Iz);
    hZ(Iz)=hZ(Iz)+dZ(Iz);
  end
% Now add dZ from incomplete bottom layer if needed
% to match zz2
  Jz=find(ZZbtm<zz2 & ZZtop<zz1 & kZb>=kk & kZt<=kk);
  if ~isempty(Jz)
    dZ = abs(ZZtop-zz2);
    smm(Jz)=smm(Jz)+SS(Jz).*dZ(Jz);
    hZ(Jz)=hZ(Jz)+dZ(Jz);
  end
% Grid cells that are completely within the depth level zz1 zz2	  
  IL=find(ZZtop<=zz1 & ZZbtm>=zz2 & kZb>=kk & kZt<=kk);
  if ~isempty(IL)
    dZ=abs(ZZbtm-ZZtop);
    hZ(IL)=hZ(IL)+dZ(IL);
    smm(IL) = smm(IL)+SS(IL).*dZ(IL); % 
  end
% Case when 1 Layer completely includes zz1 zz2
  IC=find(ZZtop>=zz1 & ZZbtm<=zz2 & kZb>=kk & kZt<=kk);
  if ~isempty(IC)
    dz=abs(zz2-zz1);
    hZ(IC)=hZ(IC)+dz;
    smm(IC) = smm(IC)+SS(IC)*dz; % 
  end

end; % for kk
I0 = find(hZ>0);
Sav(I0) = smm(I0)./hZ(I0);
Sav(Sav==0)=nan;

%keyboard


return