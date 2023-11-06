function [UUc,VVc] = sub_correctNE_OB(fina,finb);
% Corrects U and V and N and E OBs
%
%
hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
rg     = 9806;

[F,n,m,l] = read_hycom(fina,finb,'thknss');
F=F./rg;
F(F>1e20)=0;
dP=F;


% Fix Eastern OB:
[F,n,m,l] = read_hycom(fina,finb,'u-vel.'); % mean archive - total U
Iu=find(F>1e20);
F(Iu)=nan;

dp1=squeeze(dP(:,:,n-1));
dp2=squeeze(dP(:,:,n-2));
u1=squeeze(F(:,:,n-1));
u2=squeeze(F(:,:,n-2));

i02=find(dp2==0);
u2(i02)=0;
i0=find(dp1==0);
dp1(i0)=1e-6;
Tr2=u2.*dp2;
u1c=Tr2./dp1;
u1c(i0)=0;
dp1(i0)=0;

% Fix spurious U's due to
% very thin layers that are next
% to thick layers => transport increases
% distribute dlta(U) - U change from row 2 to row 1 
% over all layers
%dU0=1.5; % max difference btw U's in row 2 and row 1
%dU=u1c-u2;
%[J,I]=find(abs(dU)>dU0);
%ni=length(I);
%for kk=1:ni
%  ii=I(kk);
%  jj=J(kk);
%  ucrct=sign(dU(jj,ii))*dU0+u2(jj,ii); % correction
%  du=ucrct-u1c(jj,ii);  % change
%  u1c(jj,ii)=ucrct;
%  addU=-du*dp1(jj,ii)./sum(dp1(:,ii)); % compensate change, distribute over depth
%  u1c(:,ii)=u1c(:,ii)+addU;
%end

% Check very high U:
Umx=2;
[J,I]=find(abs(u1c)>Umx);
ni=length(I);
for kk=1:ni
  ii=I(kk);
  jj=J(kk);
  du=sign(u1c(jj,ii))*Umx-u1c(jj,ii);
  u1c(jj,ii)=u1c(jj,ii)+du;
  addU=-du*dp1(jj,ii)./sum(dp1(:,ii)); % compensate change, distribute over depth
  u1c(:,ii)=u1c(:,ii)+addU;
end


% Check transport, it should be 0
% if not, then there is mismatch
% in the # of layers in row1 and row2
% it shouldn't be too much
idP=sum(dp1); % eliminate land points, where mismatch can be
iL=find(idP==0);
iTr2=sum(u2.*dp2);
iTr1=sum(u1c.*dp1); % corrected U
dTr=iTr1-iTr2;
dTr(iL)=0;
IX=find(abs(dTr)>0.1);
np=length(IX);

if max(abs(dTr))>50,
  fprintf('Max Transp error=%9.4f\n',max(abs(dTr)));
  fprintf('sub_correctNE_OB: E Trt not conserved...\n');
  keyboard;
end

for kk=1:np
  ii=IX(kk);
  addU=dTr(ii)./sum(dp1(:,ii));
  u1c(:,ii)=u1c(:,ii)-addU;
end
u1c(i0)=0;
iTr1=sum(u1c.*dp1);
dTr=iTr1-iTr2;
dTr(iL)=0;% eliminate land points, where mismatch can be

if max(abs(dTr))>0.1
  fprintf('Check transport: E OB\n');
  keyboard
end

F(:,:,n-1)=u1c;
UUc=F;

f_chck=0;
if f_chck==1
  i=259;
  zz1=-cumsum(dp1(:,i));
  zz2=-cumsum(dp2(:,i));
  figure(1); clf;
  hold on;
  plot(u1c(:,i),zz1,'r.-'); % corrected
  plot(u1(:,i),zz1,'b.-');
  plot(u2(:,i),zz2,'g.-');
  legend('Crcted','Orig','Row2');
end



% =================
% Fix Northern OB:
% =================
[F,n,m,l] = read_hycom(fina,finb,'v-vel.'); % mean archive - total U
Iv=find(F>1e20);
F(Iv)=nan;

dp1=squeeze(dP(:,m-1,:));
dp2=squeeze(dP(:,m-2,:));
u1=squeeze(F(:,m-1,:));
u2=squeeze(F(:,m-2,:));

i02=find(dp2==0);
u2(i02)=0;
i0=find(dp1==0);
dp1(i0)=1e-6;
Tr2=u2.*dp2;
u1c=Tr2./dp1;
u1c(i0)=0;
dp1(i0)=0;

% Fix spurious U's due to
% very thin layers:
% distribute dlta(U) >1m/s over all layers
dU0=1.5;
dU=u1c-u2;
[J,I]=find(abs(dU)>dU0);
ni=length(I);
for kk=1:ni
  ii=I(kk);
  jj=J(kk);
  u1c(jj,ii)=dU0;
  dud1m=(abs(dU(jj,ii))-dU0)*dp1(jj,ii)./sum(dp1(:,ii));
  addU=sign(dU(jj,ii))*dud1m;
  u1c(:,ii)=u1c(:,ii)+addU;
end

f_chck=0;
if f_chck==1
  i=510;
  zz1=-cumsum(dp1(:,i));
  zz2=-cumsum(dp2(:,i));
  figure(1); clf;
  hold on;
  plot(u1c(:,i),zz1,'r.-'); % corrected
  plot(u1(:,i),zz1,'b.-');
  plot(u2(:,i),zz2,'g.-');
  legend('Crcted','Orig','Row2');
end

% Check transport, it should be 0
% if not, then there is mismatch
% in the # of layers in row1 and row2
% it shouldn't be too much
idP=sum(dp1); % eliminate land points, where mismatch can be
iL=find(idP==0);
iTr2=sum(u2.*dp2);
iTr1=sum(u1c.*dp1);
dTr=iTr1-iTr2;
dTr(iL)=0;
IX=find(abs(dTr)>0.1);
np=length(IX);

if max(abs(dTr))>100,
  fprintf('Max Transp error=%9.4f\n',max(abs(dTr)));
  fprintf('sub_correctNE_OB: N Trt not conserved...\n');
  keyboard;
end

for kk=1:np
  ii=IX(kk);
  addU=dTr(ii)./sum(dp1(:,ii));
  u1c(:,ii)=u1c(:,ii)-addU;
end
u1c(i0)=0;
iTr1=sum(u1c.*dp1);
dTr=iTr1-iTr2;
dTr(iL)=0;% eliminate land points, where mismatch can be

if max(abs(dTr))>0.1,
  fprintf(' Check N transport ...\n');
  keyboard;
end



F(:,m-1,:)=u1c;
VVc=F;


%keyboard

return