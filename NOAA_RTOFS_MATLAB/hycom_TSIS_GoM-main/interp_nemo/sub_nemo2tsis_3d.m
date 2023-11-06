%  Interpolating 3D fields
%
% INterpolate 1km NEMO onto HYCOM-TSIS
% Outside NEMO domain - use GLORYS
%
% Weights (or points) for interpolation are derived and saved in 
% calculate_weights.m
%
function [AAi, Amsk] = sub_nemo2tsis_3d(dnmb,pthnemo,pthtopo,pthdata,...
                                LAT,LON,HH,DX,DY,fldnm,iz0,NMI,...
                                LONN,LATN,ZZN,Amsk)

DV = datevec(dnmb);

[ah1,ah2]=size(HH);
%[DX,DY]=sub_dx_dy(LON,LAT);
%[XM,YM]=meshgrid([1:ah2],[1:ah1]);

% Get NEMO grid:
fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
if isempty(LONN)
  fprintf('Loading NEMO grid %s\n',fgrd);
  load(fgrd);
end 

%
% NEMO interpolated weights/indices:
if isempty(NMI)
  fnemo_indx = sprintf('%sNEMO_TO_HYCOM_TSIS_interp_pnts01.mat',pthdata);
  fprintf('Loading %s\n',fnemo_indx);
  NMI = load(fnemo_indx);
end


switch(fldnm)
 case('temp')
  fldnemo = 'toce';
 case('saln')
  fldnemo = 'soce';
 case('uvel')
  fldnemo = 'uoce';
 case('vvel')
  fldnemo = 'voce';
end

INN =[];
dmean = 0;  % do not demean SSH
switch(fldnm)
 case({'temp','saln'})
  AA = sub_get_NEMO_TS(dnmb,fldnemo,iz0);
 case({'uvel','vvel'})
  AA = sub_get_NEMO_UV(dnmb,fldnemo,iz0);
end
AA0 = AA;

%keyboard
%
% 
%AA = sub_fill_land(AA);
% for vel, land = 0
%switch(fldnm)
% case({'temp','saln'})
		if isempty(Amsk);
				AA = sub_fill_land(AA);
				Amsk=AA;
		else
				I=find(isnan(AA));
    AA = sub_fill_land(AA);
    amx = max(max(AA));
    if amx<-100
      AA=Amsk;
    else
  				AA(I)=0.5*(Amsk(I)+AA(I));
    end
				Amsk=AA;
		end;
% case({'uvel','vvel'})
%		I=find(isnan(AA));
%  AA(I)=0;
%end

AAi = sub_interp_nemo2hycom(LONN,LATN,LON,LAT,HH,DX,DY,NMI,AA);

f_chck=0;
if f_chck==1;
		switch(fldnm)
			case({'temp','saln'})
				c1=0;
				c2=30;
			case({'uvel','vvel'});
				c1=0;
				c2=1;
		end

  if strncmp(fldnm,'saln',4)
    c1=30;
    c2=37;
  end

  figure(1); clf;
  axes('Position',[0.07 0.2 0.4 0.7]);
  pcolor(AA0); shading flat;
  caxis([c1 c2]);
  colorbar('SouthOutside');
  axis('equal');
  set(gca,'xtick',[],...
          'ytick',[]);
  title(sprintf('%s NEMO',fldnm));

  axes('Position',[0.53 0.2 0.4 0.7]);
  pcolor(AAi); shading flat;
  caxis([c1 c2]);
  colorbar('SouthOutside');
  axis('equal');
  set(gca,'xlim',[1 700],...
          'ylim',[200 860]);
  set(gca,'xtick',[],...
          'ytick',[]);
%          'xlim',[1 ah2],...
%          'ylim',[1 ah1]);
  title('HYCOM-TSIS');
keyboard 
end 

return

