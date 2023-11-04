function TrC = sub_avrgUtsis(pthmat,pthtopo,esim,ys,ye);
% Average U vertical section
% in Yucatan Channel
% Derived from HYCOM Global and GoM reanalysis fields
%
%ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
%alat = nc_varget(ftopo,'mplat');
%elon = nc_varget(ftopo,'mplon');

TB=[];
TM=[];
TH=[];
YR=[];  % Time in terms of year days
ZZ=[];
Dst=[];
icc=0;
for iyr=ys:ye
  fmat = sprintf('%s%s_VolTrt_Yucatan_%4.4i-%4.4i.mat',...
		 pthmat,esim,iyr,iyr);
  if ~exist(fmat,'file'), % create empty arrays
    fprintf('Does not exist: %s\n',fmat);
    t1=datenum(iyr,1,1);
    t2=datenum(iyr,12,31);
    Time=[t1:t2];
    Tb=Time*nan;
    T=Time*nan;
    Time=Time(:);
    Tb=Tb(:);
    T=T(:);
  else
    fprintf('Loading %s\n',fmat);
    load(fmat);
  
    icc=icc+1;
		if icc==1
			for kb=1:4
				ZZ  = -UFLX(kb).ZZ;
% Add ghost layer for plotting
        ZZ(end+1)=ZZ(end)-0.01;
				TrC(kb).ZZ  = ZZ;
				TrC(kb).Hs  = UFLX(kb).Hs;
        TrC(kb).Dst = UFLX(kb).Dst;
        TrC(kb).Vav = TRP(kb).VvelAv;
      end
    else
      for kb=1:4
        TrC(kb).Vav = TrC(kb).Vav+TRP(kb).VvelAv;
      end
		end
	 
  end
end

for kb=1:4
  TrC(kb).Vav = TrC(kb).Vav/icc;
end

% Extend V over the bottom for plotting
for kb=1:4
  Vav = TrC(kb).Vav;
  [nz,nx] = size(Vav);
  Vav(isnan(Vav)) = 0;
  TrC(kb).Vav=Vav;
end

return
