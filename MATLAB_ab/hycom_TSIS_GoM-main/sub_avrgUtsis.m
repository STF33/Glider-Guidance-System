function TrC = sub_avrgUtsis(pthmat,pthtopo,esim,ys,ye);
% Average U vertical section
% in Yucatan Channel
% For HYCOM native vertical grid
% with layer thicknesses
ftopo = sprintf('%sias_gridinfo.nc',pthtopo);
alat = nc_varget(ftopo,'mplat');
elon = nc_varget(ftopo,'mplon');

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
				TrC(kb).ZZ = UFLX(kb).ZZ_intrp;
				TrC(kb).Hs = UFLX(kb).Hs;
        TrC(kb).Vav = TRP(kb).VvelAv;

        if isfield(UFLX,'Dst') & ~isempty(UFLX(kb).Dst)
          TrC(kb).Dst = UFLX(kb).Dst;
        else
					ni = length(UFLX(kb).IJ);
					for ii=1:ni
						i0=UFLX(kb).IJ(ii,1);
						j0=UFLX(kb).IJ(ii,2);
						dx=distance_spheric_coord(alat(j0,i0),...
               elon(j0,i0),alat(j0,i0+1),elon(j0,i0+1));
						TrC(kb).Dst(ii)=dx;
					end
        end

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

return
