% subsample TSIS domain to 
% GoM domain similar to NEMO
% and izobath Z0 if Z0<0
%
function  [XX,YY,HHs,INH,Iocn,xt1,xt2,yt1,yt2,its1,its2,jts1,jts2] = ...
            sub_subsample_HYCOM2GoM(HH,LON,LAT,Z0);

[m,n]=size(HH);

% GoM region HYCOM:
% TSIS domain
GOM=[366   489
   476   531
   583   560
   576   646
   508   827
   336   848
   204   829
    64   798
    19   746
    16   662
    12   578
    25   455
    71   382
   165   356
   281   400];

[XM,YM]=meshgrid([1:n],[1:m]);
INH = inpolygon(XM,YM,GOM(:,1),GOM(:,2));
Iocn = find(INH==1 & HH<Z0);
clear XM YM

% 
% Subsample to a smaller domain:
xnas1 = min(LON(Iocn));
xnas2 = max(LON(Iocn));
ynas1 = min(LAT(Iocn));
ynas2 = max(LAT(Iocn));
[J,I]=find(LON<xnas1);
its1=max(I);
[J,I]=find(LON>xnas2);
its2=min(I);
[J,I]=find(LAT<ynas1);
jts1=max(J);
[J,I]=find(LAT>ynas2);
jts2=min(J);

xt1=LON(jts1,its1);
yt1=LAT(jts1,its1);
xt2=LON(jts2,its2);
yt2=LAT(jts2,its2);

XX=LON(jts1:jts2,its1:its2);
YY=LAT(jts1:jts2,its1:its2);
HHs=HH(jts1:jts2,its1:its2);


return
