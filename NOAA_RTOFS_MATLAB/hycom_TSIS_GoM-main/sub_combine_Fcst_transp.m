function TrC = sub_combine_Fcst_transp(fmat);
% Combines time series of Volume Transports
% in Yucatan Channel
% extracted for individual years
TB=[];
TM=[];
TH=[];

for iyr=ys:ye
%  fmat = sprintf('%s%s_VolTrt_Yucatan_%4.4i-%4.4i.mat',...
%		 pthmat,esim,iyr,iyr);
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
end

ndy=365; % no leap years!
DV=datevec(TM);
dj1=datenum(DV(1,1),1,1);
d1=TM(1)-dj1+1;
dj2=datenum(DV(end,1),1,1);
d2=TM(end)-dj2+1;
dstrt=DV(1,1)+d1/ndy;
dend=DV(end,1)+d2/ndy;
YR=[dstrt:1/ndy:dend];
YR=YR(:);

TrC.Time=TM;
TrC.Barotrop=TB;
TrC.Total=TH;
TrC.YR=YR;
return