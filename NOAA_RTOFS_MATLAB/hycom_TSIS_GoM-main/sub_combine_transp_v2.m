function TrC = sub_combine_transp(pthmat,esim,ys,ye);
% Combines time series of Volume Transports
% in Yucatan Channel
% extracted for individual years
TB=[];
TM=[];
TH=[];

for iyr=ys:ye
  fmat = sprintf('%s%s_VolTrt_Yucatan_%4.4i-%4.4i.mat',...
		 pthmat,esim,iyr,iyr);
  if ~exist(fmat,'file'), % create empty arrays
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
end

%ndy=365; % no leap years!
nTM = length(TM);
for ll = 1:nTM
  DV=datevec(TM);
  iyr = DV(1,1);
  dj1=datenum(iyr,1,1);
  ndays_year = datenum(iyr,12,31)-dj1+1;
  d1=TM(ll)-dj1+1;
  YR(ll,1) = iyr+d1/ndays_year;
end
YR=YR(:);

TrC.Time=TM;
TrC.Barotrop=TB;
TrC.Total=TH;
TrC.YR=YR;
%keyboard
return
