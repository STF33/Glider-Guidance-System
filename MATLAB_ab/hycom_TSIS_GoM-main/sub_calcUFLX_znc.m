% Calculate UFLX for a "box" along the sides
% U/V are collocated, interpolated on Z, in netcdf
function [UFLX,TRP] = sub_calcUFLX_znc(fina,UFLX,NB,TRP,cc,alat,elon,dd);


SSH=squeeze(nc_varget(fina,'surf_el'));

for kb=1:NB
	OB = UFLX(kb).OB;
	js1= UFLX(kb).IJ(1,2);
	js2= UFLX(kb).IJ(end,2);
	is1= UFLX(kb).IJ(1,1);
	is2= UFLX(kb).IJ(end,1);

	clear ssh
	if strcmp(OB,'S') | strcmp(OB,'N')
		j0=UFLX(kb).IJ(1,2);
		II=UFLX(kb).IJ(:,1);

		for ip=1:length(II)
			i0=II(ip);
			ssh(ip)=0.5*(SSH(j0,i0)+SSH(j0-1,i0));
		end;

	elseif strcmp(OB,'E') | strcmp(OB,'W');
		i0=UFLX(kb).IJ(1,1);
		II=UFLX(kb).IJ(:,2);
		for ip=1:length(II)
			j0=II(ip);
			ssh(ip)=0.5*(SSH(j0,i0)+SSH(j0,i0-1));
		end

	end
	UFLX(kb).SSH=ssh;
end;

%  =================
% Create array of interface depths from layer thicknesses:
%  =================
ZZ=nc_varget(fina,'depth');
nl=length(ZZ);
dZ(1,1)=0.5*(ZZ(2)-ZZ(1));
for kk=2:nl-1
	dZ(kk,1)=0.5*abs(ZZ(kk+1)-ZZ(kk-1));
end
dZ(nl)=0.5*abs(ZZ(nl)-ZZ(nl-1));

for ks=1:NB
	UFLX(ks).ZZ=ZZ;
end


% Note that in the raw binary, U, V are on u,v grid
% h is on p-grid
% !!! NOTE: in the instanteneous output files archv u-vel, v-vel 
% are BAROCLINIC (not total) components
% need to add u barotropic !!!
%  In nc z-interpolated U,V H collocated
% U = total U
%keyboard

% Get velocities along the sides of the box
for kb=1:NB
	OB = UFLX(kb).OB;
	js1=UFLX(kb).IJ(1,2);
	js2=UFLX(kb).IJ(end,2);
	is1=UFLX(kb).IJ(1,1);
	is2=UFLX(kb).IJ(end,1);

	dj=js2-js1+1;
	di=is2-is1+1;
	U = squeeze(nc_varget(fina,'water_u',...
			 [0 0 js1-1 is1-1],[1 nl dj di]));
	V = squeeze(nc_varget(fina,'water_v',...
			 [0 0 js1-1 is1-1],[1 nl dj di]));

	% Barotropic U & V:
	dmm=U';
	dmm(isnan(dmm))=0;
	I=dmm*0+1;
	I(dmm==0)=0;
	Hbz=I*dZ; % Bottom depth in z cells
	Ubrt=dmm*dZ./Hbz;

	dmm=V';
	dmm(isnan(dmm))=0;
	I=dmm*0+1;
	I(dmm==0)=0;
	Vbrt=dmm*dZ./Hbz;

	UFLX(kb).sumDZ=Hbz; % integrated Z cells to bottom 

	if strcmp(OB,'S') | strcmp(OB,'N')
		UFLX(kb).V=V;
		UFLX(kb).UBRT=Vbrt';
	elseif strcmp(OB,'E') | strcmp(OB,'W')
		UFLX(kb).V=U;
		UFLX(kb).UBRT=Ubrt';
	end

end


% For plotting purposes add extra layer to show the bottom:
for kb=1:NB
	V=UFLX(kb).V;  % keep 
	[a1,a2]=size(V);
	for ii=1:a2
		V(a1+1,ii)=V(a1,ii);
	end;

	if ~isfield(UFLX,'Vav') | isempty(UFLX(kb).Vav)
		UFLX(kb).Vav=V*0;
	end
end


% Calculate flux:
Ttot=0;
for kb=1:NB
	js1 = UFLX(kb).IJ(1,2);
	js2 = UFLX(kb).IJ(end,2);
	is1 = UFLX(kb).IJ(1,1);
	is2 = UFLX(kb).IJ(end,1);

	ni1 = is2-is1+1;
	nj1 = js2-js1+1;

	Hs  = UFLX(kb).Hs;
%      thk = UFLX(kb).thk;
	V   = UFLX(kb).V;
	SSH = UFLX(kb).SSH;
	Vav = UFLX(kb).Vav;
	Ubt = UFLX(kb).UBRT;
%      sthk= UFLX(kb).sumthk;
%    SSH=SSH*0;  % ignore SSH

% Note that positive influx is "IN" the volume
% so at the Northern & Eastern bndry, posit. flux is opposite to +V/U
	ob=UFLX(kb).OB;
	if strcmp(ob,'N')>0 | strcmp(ob,'E')>0
		sgn=-1;
	else
		sgn=1;
	end

	VF=V*nan;
	D1=Hs*0;
	ni=max([ni1,nj1]);

	Dst=UFLX(kb).Dst;
  if isempty(Dst)
    Dst=zeros(ni,1);
  end
%keyboard
% ===========================
% FLUX CALCULATION:
% ===========================
	clear AdZ
	for ii=1:ni
		i0=UFLX(kb).IJ(ii,1);
		j0=UFLX(kb).IJ(ii,2);
		if ii<ni
			i1=UFLX(kb).IJ(ii+1,1);
			j1=UFLX(kb).IJ(ii+1,2);
		else
			i1=UFLX(kb).IJ(ii-1,1);
			j1=UFLX(kb).IJ(ii-1,2);

		end

		if Dst(ni)==0;
			dx=distance_spheric_coord(alat(j0),elon(i0),...
				 alat(j1),elon(i1));
			Dst(ii)=dx;
		end

		if ii>1 & ii<ni
			dx=Dst(ii);
		else
			dx=0.5*Dst(ii);
		end

		dnn=0;
		vv=V(:,ii);
		dz=dZ;
		dz(1)=dz(1)+SSH(ii);
		VF(:,ii)=sgn*vv.*dz.*dx;
		in=find(isnan(vv));
		dz(in)=nan;
		AdZ(ii)=nansum(dz.*dx); % area of w colmn * dx
	end
	TRbtp=sgn*(AdZ.*UFLX(kb).UBRT);
	Tbtp=nansum(TRbtp);
%
% Average U - layer thickness is fixed on Z grid
	for ii=1:ni
		Vav(1:nl,ii)=Vav(1:nl,ii)+V(:,ii);  % U
	end


	UFLX(kb).Vav=Vav;
	UFLX(kb).dZ = dZ;
%      UFLX(kb).THKav=sthk;  % averaged layer thickness
	UFLX(kb).counter=cc;  % # of fields averaged
  UFLX(kb).Dst=Dst;

	TRP(kb).VvelAv = V;    % U profile
	TRP(kb).TrLayers(cc,:)=nansum(VF,2);
	Ttot=nansum(TRP(kb).TrLayers(cc,:));
	TRP(kb).DepthPrf = Hs;
	TRP(kb).Time(cc) = dd;
	TRP(kb).TBtp(cc) = Tbtp; % transport from Barotropic U only
%keyboard
  if kb==1, fprintf('         Hourly  %s\n',datestr(dd)); end
	fprintf('====> KB=%i, NET Transport T=%5.2f Sv\n',kb, Ttot*1e-6);
	fprintf('====> Max Uav=%5.2f \n',...
	 nanmax(nanmax(abs(V))));
end;



return 
