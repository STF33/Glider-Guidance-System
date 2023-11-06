function [UFLX,TRP] = sub_UFLX_TRP(BB);
%   For not collocated U&V !!!
%
% Preapre structured arrays for 
% Transport calculation within
% the box if NB=4
%
% Southern OB:
% Note that to make contour closed (for not collocated U,V)
% V point on Southern OB needs to end 1 grid point 
% before first U point on the Eastern OB
% also, V points are 1 grid cell down
esim  = BB.esim;
iSW   = BB.iSW;
jSW   = BB.jSW;
iSE   = BB.iSE;
jSE   = BB.jSE;
iNE   = BB.iNE;
jNE   = BB.jNE;
iNW   = BB.iNW;
jNW   = BB.jNW;
HH    = BB.HH;
NB    = BB.NB; % # of bndrs where flux is calculated, if NB<4, OBs>NB are not counted
f_plt = BB.fplt;

is1=[iSW:iSE-1]';
js1=ones(size(is1))*jSW;  % Make sure v is nan at the 1st point along Southern OB
Hs1=squeeze(HH(js1(1),is1(1):is1(end)));
Hs2=squeeze(HH(js1(1)-1,is1(1):is1(end)));
UFLX(1).esim=esim;
UFLX(1).OB='S';
UFLX(1).IJ=[is1,js1];
UFLX(1).Hs=0.5*(Hs1+Hs2);


% Eastern OB:
% First grid point on the EB with U directed into the contoured volume
% is +1 grid cell from the northern OB
js1=[jSE:jNE-1]';
is1=ones(size(js1))*iSE;
%is1=ones(size(js1))*535;
%is1=ones(size(js1))*300;
Hs1=squeeze(HH(js1(1):js1(end),is1(1)));
Hs2=squeeze(HH(js1(1):js1(end),is1(1)-1));
UFLX(2).esim=esim;
UFLX(2).OB='E';
UFLX(2).IJ=[is1,js1];
UFLX(2).Hs=0.5*(Hs1+Hs2);

% Northern OB:
% 1 point less in Eastern Dir.
% than U point on the Eastern OB
is1=[iSW:iNE-1]';
js1=ones(size(is1))*jNE;
Hs1=squeeze(HH(js1(1),is1(1):is1(end)));
Hs2=squeeze(HH(js1(1)-1,is1(1):is1(end)));
UFLX(3).esim=esim;
UFLX(3).OB='N';
UFLX(3).IJ=[is1,js1];
UFLX(3).Hs=0.5*(Hs1+Hs2);


% Western OB:
% First grid point on the West OB with U directed into the contoured volume
% is +1 grid cell from the northern OB
js1=[jSW:jNW-1]';
is1=ones(size(js1))*iSW;
Hs1=squeeze(HH(js1(1):js1(end),is1(1)));
Hs2=squeeze(HH(js1(1):js1(end),is1(1)-1));
UFLX(4).esim=esim;
UFLX(4).OB='W';
UFLX(4).IJ=[is1,js1];
UFLX(4).Hs=0.5*(Hs1+Hs2);


for kp=1:NB
  ob=UFLX(kp).OB;
  TRP(kp).OB=ob;
  TRP(kp).TDav=[];
  TRP(kp).VvelAv=[];
  TRP(kp).DepthPrf=[];
  TRP(kp).Time=[];
end;

% Plot x-section:
if f_plt>0
  figure(5); clf;
  hold on;
  contour(HH,[0 0],'k','linewidth',2);
  contour(HH,[-1000 -1000],'Color',[0.5 0.5 0.5]);

  plot([iSW iSE],[jSW jSE],'r','linewidth',1);
  plot([iSE iSE],[jSE jNE],'r','linewidth',1);
  plot([iSW iSE],[jNE jNE],'r','linewidth',1);
  plot([iSW iSW],[jSW jNW],'r','linewidth',1);
  drawnow
end



return