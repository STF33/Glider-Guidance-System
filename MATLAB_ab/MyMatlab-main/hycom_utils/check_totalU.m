function utotal = check_totalU(fina,finb,clc);
% ------------
%  utotal = check_totalU(fina,finb,clc);
%
% Check if u-vel. is
% utotal = ubarocl + u barotrop (depth-averaged)
% the total baroclinic velocity or baroclinic anomalies:
% in instanteneous archive files
% u-vel. = utot - ubtrop
% clc - collocated (1) or not (0) U,V, dH
% in mean fields - collocated
% in archm u_vel=utot;
% ------------
fprintf('Checking if utotal = ubrcl+ubrtrop in %s\n',fina);

rg = 9806;

[F,n,m] = read_hycom(fina,finb,'u_btrop');
F(F>1e6)=nan;
Ub=squeeze(F);

%[F,n,m] = read_hycom(fina,finb,'v_btrop');
%F(F>1e6)=nan;
%Vb=squeeze(F);

[F,n,m,l] = read_hycom(fina,finb,'u-vel.');
F(F>1e6)=nan;
U=squeeze(F);

%[F,n,m,l] = read_hycom(fina,finb,'v-vel.');
%F(F>1e6)=nan;
%V=squeeze(F);

[F,n,m,l] = read_hycom(fina,finb,'thknss');
F=F/rg;
F(F>1e10)=nan;
dH=squeeze(F); % note this is not U,V thickness, center grid
	       % there are thck-u & thck-v 
ih=find(dH<1e-1);
U(ih)=nan;
V(ih)=nan;

HH = -abs(squeeze(nansum(dH,1)));
I=find(HH<-2000);
i0=I(10);
[jj,ii] = ind2sub([m,n],i0);


% collocate U and H
uu=squeeze(U(:,jj,ii));
%dh=squeeze(dH(:,jj,ii));
if clc==0
  dh1=dH(:,jj,ii-1);
  dh2=dH(:,jj,ii);
  dh=0.5*(dh1+dh2);
else
  dh = dH(:,jj,ii);
end
uua=nansum(uu.*dh)./nansum(dh);
ub=Ub(jj,ii);
% uua and ub should be close:
if abs(uua-ub)/abs(ub)<0.01;
  fprintf('\n');
  fprintf('Depth-intgrated anom U=%6.4f, brtrp U=%6.4f\n',uua,ub);
  fprintf('u_vel. and v_vel. are total baroclinic+barotrop U\n');
  utotal=1;
else
  fprintf('\n');
  fprintf('Depth-integrated anom U=%6.4f, brtrp U=%6.4f\n',uua,ub);
  fprintf('u_vel. and v_vel. are barocl. anomlaies U\n');
  utotal=0;
end


return
