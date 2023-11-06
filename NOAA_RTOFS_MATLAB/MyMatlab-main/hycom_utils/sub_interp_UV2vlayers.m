    function [Uout,Vout] = sub_interp_UV2vlayers(Uin,Vin,dPin,dPout,HH);
% Interpolate U and V velocity components on HYCOM hybrid layers
% from N layers to M layers
% Input: all huge values -> nans
%   Uin,Vin - are Layer U/V,  can be u_total = u_bacl+u_brtp
%   or u_bacl
%   Zlev - array of fixed depths
%   dP - layer thicknesses, Pa or m
%    
%   HH - bathymetry corresponding to dP
%
%  In this method, layer thicknesses are converted into depths of
% layer interfaces and treated as fixed "z-levels" for any given point
%  
[lin,mm,nn]=size(dPin);
[lout,mout,nout]=size(dPout);
fprintf('\n Interpolating U,V %i layers ---> %i layers ...\n',lin,lout);

if mout~=mm | nout~=nn
  fprintf('Input grid i,j   old: %i, %i, new: %i, %i\n',nn,mm,nout,mout);
  error('Dimensions of the New and old Input grids do not agree');
end


rg = 9806;
hg=0.1*2^100;
onemm=rg/1000; % HYCOM threshold for alsmot mass-less layer
%Tv= 72; % topography version

%pthtopo = '/Net/ocean/ddmitry/HYCOM/GoM/topo/GOMl0.04/';
%ftopo = sprintf('%s/depth_GOMl0.04_%2.2i.nc',pthtopo,Tv); % 
%HH  = nc_varget(ftopo,'Bathymetry');
%LON = nc_varget(ftopo,'Longitude');
%LAT = nc_varget(ftopo,'Latitude');
[JD,ID]=size(HH);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

% Find layer mid-depths:
%dz = diff(abs(Zlev));
%nZl = length(Zlev);
%ZMlev = zeros(nZl-1,1);
%for kk=1:nZl-1
%  ZMlev(kk,1)=Zlev(kk)-0.5*dz(kk);
%end;
%Zlev0=Zlev;
%DZ0=dz;

dPin(dPin>=0.1*hg) = nan; 
dPout(dPout>=0.1*hg) = nan;  % Pressure, Pa

% COnvert Pa -> m if needed
md=nanmax(nanmax(nanmax(abs(dPin))));
if md>2*rg,
  dPin=dPin./rg;  % Pa -> m
end;
md=nanmax(nanmax(nanmax(abs(dPout))));
if md>2*rg,
  dPout=dPout./rg;  % Pa -> m
end;

% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
fprintf('Calculating interface depths ...\n');
[ZZin,ZMin] = sub_thck2dpth(dPin);
[ZZout,ZMout] = sub_thck2dpth(dPout); % 

% Check:
I=find(abs(HH)>2000,1);
h0=HH(I);
hin=max(abs(ZZin(:,I)));
hout=max(abs(ZZout(:,I)));

if abs(1-abs(h0/hin))>0.01
  fprintf('Check depth HH= %10.2f, ZZin = %10.2f\n',h0,hin);
  error('Check dPin: interface depths do not sumup to HH');
end

if abs(1-abs(h0/hout))>0.01
  fprintf('Check depth HH= %10.2f, ZZout = %10.2f\n',h0,hout);
  error('Check dPin: interface depths do not sumup to HH');
end


%  -------------------------------------------------- 
%  For each interface in new t.dens. find
%  nearest upper and bottom interfaces in old t.dens.
%     Zin
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
Uout = zeros(lout,mout,nout)*nan;
Vout = zeros(lout,mout,nout)*nan;

for io=1:length(IOcean);
  
  if mod(io,5000)==0;
    fprintf(' === U,V z -> Sigma0:  Done %8.2f%%\n',io/length(IOcean)*100);
  end
  
  
  II=IOcean(io);
  [jj,ii]=ind2sub(size(HH),II);
  Zout = squeeze(ZZout(:,jj,ii));
  Uin0 = squeeze(Uin(:,jj,ii));
  Vin0 = squeeze(Vin(:,jj,ii));
  Zin  = squeeze(ZZin(:,jj,ii));
  nS   = length(Zout);
  nZl  = length(Zin);
  hb   = HH(jj,ii);
  dpin = squeeze(dPin(:,jj,ii));
  dpout= squeeze(dPout(:,jj,ii));
  
% Check for spurious U,V:
  IS=find(abs(Uin0>3));
  JS=find(abs(Vin0>3));
  if ~isempty(IS)
    dmm=Uin0;
    dmm(IS)=nan;
    Uin0(IS)=sign(Uin0(IS)).*nanmean(dmm);
  end
  if ~isempty(JS)
    dmm=Vin0;
    dmm(JS)=nan;
    Vin0(JS)=sign(Vin0(JS)).*nanmean(dmm);
  end
  
  Isf=min(find(Uin0==0)); % land
  if Isf==1;  % U,V point on land
    continue;
  end
% Interpolate into grid mid-depths:
% Fill Nan's with something  
%  Ib = min(find(isnan(Uin0)));
%  if ~isempty(Ib)
%    Uin0(Ib:end)=9.99e-5;  % fill NaN bottom values
%    Vin0(Ib:end)=9.99e-5;  % fill NaN bottom values
%  end
  
%  Uin = interp1(ZZout,Uin0,ZMin,'cubic');
%  Vin = interp1(ZZout,Vin0,ZMin,'cubic');

% U,V velocity components in HYCOM hybrid layers (sigma0 dens.)  
  Us = [];
  Vs = [];
  for k=1:nS-1
    if dpout(k)<1e-3, continue; end; % 0-thickness layer
    s1=Zout(k);
    s2=Zout(k+1);
    
% Find corresponding old "Z surfaces"
% which are top surfaces of the 
% old hybrid-layers, that contain new HYCOM hybdrid surfaces
    k1=max(find(Zin>=s1));  % top HYCOM surface
    k2=max(find(Zin>s2)); % btm HYCOM surface

    if k1>k2
      error('Finding corresponding k1 and k2 for new hybrid layer: k1>k2');
    end
    
    
% How many complete z-levels of old grid inside the new hybrid layer:
    nZ = k2-k1-1;
   
% ------------------------------------    
    chck=0;
    if chck==1
%      keyboard
      z1=Zin(k1);
      z2=Zin(k2);
    
      figure(10); clf;
%      xx=ones(nZ,1);
      hold on;
      for kz=1:nZl
        plot([0.98 1.02],[Zin(kz) Zin(kz)],'b-');
      end
      xx=ones(nS,1);
      for kz=1:nS
        plot([0.98 1.02],[Zout(kz) Zout(kz)],'r-');
      end
      plot([0.9 1.1],[z1 z1],'c--','linewidth',1.5); % old HYCOM layer
      plot([0.9 1.1],[z2 z2],'c--','linewidth',1.5);
      plot([0.9 1.1],[s1 s1],'m:','linewidth',1.5); % new HYCOM layer
      plot([0.9 1.1],[s2 s2],'m:','linewidth',1.5);
    end
% ------------------------------------    

%
% This is simplified version of U,V interpolation
% More accurate approach: for each grid cell u=U(z) - a linear
% function between the interfaces
% Then the transport between several grid cells (where bottom and top Z cells 
% are only partial - where HYCOM interfaces are) should be
% intgral(zb:Zinterface(k1)) of U(z) + Integral of full cells of U(z) + 
% partial integral(Zinterface(k2):ztop) of U(z)
    
    if nZ<0 % HYCOM new layer is within 1 old layer    
            % same U and V
      Us(k,1) = Uin0(k1);
      Vs(k,1) = Vin0(k1);
    else % upper isop. surf is in layer k, bottm - layer k2
% Delta thickness of partial layers      
      dZS1=abs(s1-Zin(k1+1));
      dZS2=abs(s2-Zin(k2));
% Full layer thickness:
%      dZ1=abs(Zin(k1)-Zin(k1+1));
%      dZ2=abs(Zin(k2-1)-Zink(k2));
% Partial transport
      if k>1
        pU1=Uin0(k1)*dZS1; 
        pV1=Vin0(k1)*dZS1;
      else % surface, no partial transport above surface
        pU1=0;
	pV1=0;
      end
      
      pU2=Uin0(k2)*dZS2;
      pV2=Vin0(k2)*dZS2;
      
      dZfull=0;
      Ufull=0;
      Vfull=0;
      cc=0;
      for kf=k1+1:k2-1
	cc=cc+1;
	dZ1=abs(Zin(kf)-Zin(kf+1));
	dZfull=dZfull+dZ1;
	Ufull=Ufull+Uin0(kf)*dZ1;
	Vfull=Vfull+Vin0(kf)*dZ1;
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
  
  Uout(:,jj,ii)=Us;
  Vout(:,jj,ii)=Vs;
  
%  if hb<-3000; keyboard; end;

% Check transport:
%  Izr = min(find(Uin == 0));
%  if ~isempty(Izr),
%    ibtm=Izr-1;
%  else
%    ibtm=nZl-1;
%  end
%keyboard
  Uin0(isnan(Uin0))=0;
  UZtot=Uin0'*dpin;
  Us(isnan(Us))=0;
  UStot=Us'*dpout;
  RUtrt=abs(1-UStot/UZtot);

  Vin0(isnan(Vin0))=0;
  VZtot=Vin0'*dpin;
  Vs(isnan(Vs))=0;
  VStot=Vs'*dpin;
  RVtrt=abs(1-VStot/VZtot);
  
  if RUtrt>1 | RVtrt>1
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
  Zout  = squeeze(ZZout(:,jj,ii));
  Uin0  = squeeze(Uin(:,jj,ii)); % U original old hbrd lrs
  Vin0  = squeeze(Vin(:,jj,ii)); % V original old hybrid lrs
  Zin   = squeeze(ZZin(:,jj,ii));
  
  figure(10); clf;
  
  axes('Position',[0.08 0.06 0.35 0.85]);
  hold on;
  plot(Uin0,Zin,'g-'); % Original U on z, cell interfaces
  % Plot U interpolated into new HYCOM v. grid
  for k=1:nS-1
    u=Uout(k,jj,ii);
    plot([u u],[Zout(k) Zout(k+1)]);
  end
  hb=HH(jj,ii);
  plot([min(Uin) max(Uin)],[hb hb],'k--');
  title('U, OrigZ(g), MidZ(r), HYCOM(b)' );
  set(gca,'ylim',[1.05*hb 0])
  
  axes('Position',[0.52 0.06 0.35 0.85]);
  hold on
  plot(Vin0,Zin,'g-');
  for k=1:nZl-1
    v=Vin(k);
    plot([v v],[Zin(k) ZZin(k+1)],'r-');
  end
  for k=1:nS-1
    v=Vout(k,jj,ii);
    plot([v v],[Zout(k) Zout(k+1)]);
  end
  hb=HH(jj,ii);
  plot([min(Vin) max(Vin)],[hb hb],'k--');
  title('V, HYCOM: blue');
  set(gca,'ylim',[1.05*hb 0])
end

    

return