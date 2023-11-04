function TrC = sub_combine_transp(pthmat,esim,ys,ye);
% Combines time series of Volume Transports
% in Yucatan Channel
% extracted for individual years
TB=[];
TM=[];
TH=[];
YR=[];  % Time in terms of year days

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
  
    aa=TRP(1).TrLayers;
    T1=nansum(aa,2); % layer averaged
    Tb1=TRP(1).TBtp;
    aa=TRP(2).TrLayers;
    T2=nansum(aa,2); % layer averaged
    Tb2=TRP(2).TBtp;
    T=T1+T2;
    Tb=Tb1+Tb2;

    Time=TRP(1).Time;
    Time=Time(:);
    Tb=Tb(:);
    TH=TH(:);
%keyboard
  end
 
 
% Missing Jan 1, some years:
  dj1=datenum(iyr,1,1);
  i1=find(Time==dj1);
  if isempty(i1),
    Time=[dj1;Time];
    TB=[NaN;TB];
    T=[NaN;T];
  end

  TM=[TM;Time];
  TB=[TB;Tb];
  TH=[TH;T];

  DV=datevec(Time);
  ndy = datenum(iyr,12,31)-dj1+1;
  yrdd = DV(1,1)+(Time-Time(1))/ndy;
  YR=[YR;yrdd];
end

%ndy=365; % no leap years!
%DV=datevec(TM);
%dj1=datenum(DV(1,1),1,1);
%d1=TM(1)-dj1+1;
%dj2=datenum(DV(end,1),1,1);
%d2=TM(end)-dj2+1;
%dstrt=DV(1,1)+d1/ndy;
%dend=DV(end,1)+d2/ndy;
%YR=[dstrt:1/ndy:dend];
YR=YR(:);


TrC.Time=TM;
TrC.Barotrop=TB;
TrC.Total=TH;
TrC.YR=YR;
%keyboard
return
