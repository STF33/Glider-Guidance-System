function Sav = sub_binav2z_1D(Hbtm,ZZ,ZZf,F,f_chck);
% Bin-Average scalar field 1D array 
% from HYCOM layers to fixed z level
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
% Input: Hbtm - HYCOM topo, local depth
% ZZf - fixed depth layers within which field is averaged
% ZZ - 1D array of HYCOM interface depths
% F  - scalar 1D array
% f_chck - draws profiles, stops the subroutine at the end
%
%fprintf('sub_binav2z_1D: bin averaging\n');

llz = length(ZZ);
ll  = llz-1;
lzf = length(ZZf);
Sav = zeros(lzf-1,1)+1e20;

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

for kzf=1:lzf-1 % Fixed depth layers
  zz1  = ZZf(kzf);   % top interface fixed Z lev
  zz2  = ZZf(kzf+1); % bottom interface fixed Z lev
  if (abs(Hbtm)<abs(zz1)), break; end; % below bottom
%cZk  = HH*0; % count v. layers
%
  kZt  = 0;
  kZb  = 0;
  for kk=1:ll  % find bounding Z layers fixed depths
    ZZtop = ZZ(kk);   % top intrf HYCOM
    ZZbtm = ZZ(kk+1); % btm intrf HYCOM

  % Too deep:
    if ZZtop<zz2, break; end;

  % Layers:  
    if (ZZbtm>zz1); continue; end % HYCOM layer too shallow

% Find end of integration - HYCOM bottom intrf
    if (Hbtm<zz2 & ZZtop>zz2 & kZt>0)
      kZb=kk;
    end

% Start top index when top& bottom HYCOM intrf above & below Z top intrf
    if (Hbtm<zz2 & ZZtop>=zz1 & ZZbtm<zz1 & kZt==0);
      kZt=kk; % keep 1st layer
      kZb=kk; % indexing is kept by the top index 
    end
  end
  
  if kZb==0, break; end; % bottom 
%  fprintf('** Zlr %i, zz1-zz2=%6.3f:%6.3f\n',kzf,zz1,zz2);
%  fprintf('** HYCOM: kZt=%i, kZb=%i,%6.3f:%6.3f\n',kZt,kZb,ZZ(kZt),ZZ(kZb));
%  keyboard

  %kZb(kZb==0)=nan; % bottom layer within the depth interval
  %kZt(kZt==0)=nan;  % surface layer
%  Sav = 0; % S mean
  hZ  = 0;
  smm = 0;
  for kk=kZt:kZb  % HYCOM layers that contain fixed Z layer
  %  I=find(kZb>=kk & kZt<=kk);
  %  if isempty(I), continue; end;
  %  SS = HH*nan;
  %  SS(I) = squeeze(F(kk,I)); % S at level kk, only pnts that fall within zz1 zz2
    SS    = F(kk); % 
    ZZtop = ZZ(kk); % top intrf of model layer
    ZZbtm = ZZ(kk+1); % btm intrf model layer

  % Add dZ from incomplete top layer if needed
  % only for HYCOM layers that are partially 
  % within the Z-layer
    if (ZZtop>zz1 & ZZbtm>zz2);
      dZ  = abs(ZZbtm-zz1);
      smm = smm+SS.*dZ;
      hZ  = hZ+dZ;
    end
  % Now add dZ from incomplete bottom layer if needed
  % to match zz2
    if (ZZbtm<zz2 & ZZtop<zz1);
      dZ  = abs(ZZtop-zz2);
      smm = smm+SS*dZ;
      hZ  = hZ+dZ;
    end
  % Grid cells that are completely within the depth level zz1 zz2	  
    if (ZZtop<=zz1 & ZZbtm>=zz2);
      dZ  = abs(ZZbtm-ZZtop);
      hZ  = hZ+dZ;
      smm = smm+SS.*dZ; % 
    end
  % Case when 1 Layer completely includes zz1 zz2
    if (ZZtop>=zz1 & ZZbtm<=zz2)
      dz  = abs(zz2-zz1);
      hZ  = hZ+dz;
      smm = smm+SS*dz; % 
    end

  end; % for kk - HYCOM layers
  Sav(kzf,1) = smm./hZ;
end % kzf - fixed z layers

Sav(Sav>1e10)=nan;
%fprintf('Averaging completed\n');

%f_chck=0;
if f_chck==1
  figure(11); clf;
  hold on;
  
  z0=ZZ(1);
  f0=F(1);
  for kk=1:llz-1
    z1=ZZ(kk);
    z2=ZZ(kk+1);
    f1=F(kk);
    plot([f1 f1],[z1 z2],'r');
    plot([f1 f0],[z1 z0],'r');
    z0=z2;
    f0=f1;
%    fprintf('kk=%i,z0=%6.2f,f0=%2.6f\n',kk,z0,f0);
  end

  z0=ZZf(1);
  f0=Sav(1);
  for kkz=1:lzf;
    z1=ZZf(kkz);
    z2=ZZf(kkz+1);
    f1=Sav(kkz);
    if isnan(f1), break; end;
    plot([f1 f1],[z1 z2],'b');
    plot([f1 f0],[z1 z0],'b');
    z0=z2;
    f0=f1;
  end    
  set(gca,'xgrid','on',...
	  'ygrid','on');
%  keyboard
  
end
%keyboard


return