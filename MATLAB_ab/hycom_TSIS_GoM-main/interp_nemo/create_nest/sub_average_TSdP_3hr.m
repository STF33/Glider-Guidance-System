       function [Tav,Sav,dPav] = sub_average_TSdP(pthrlx,dnmb);
% 
% average T,S, dP fields on HYCOM grid
% from 3-hourly fields -> daily mean
% to avoid the problem of wierd isopycnals
%
DV=datevec(dnmb);
iyr=DV(1);
dj1=datenum(iyr,1,1);
iday=dnmb-dj1+1;

rg = 9806;
hg = 2^100;
hgg= 1e20;

Tav=[];
Sav=[];
dPav=[];

cc=0;
for hr=0:3:21
  fina = sprintf('%srelax.%4.4i%3.3i%2.2i.a',pthrlx,iyr,iday,hr);
  finb = sprintf('%srelax.%4.4i%3.3i%2.2i.b',pthrlx,iyr,iday,hr);
  
  fprintf('    sub_aveage_TSdP: averaging %2.2ihr\n',hr);
  
  if ~exist(fina,'file'), 
    fprintf('Missing: %s\n',fina);
    continue; 
  end;
  cc=cc+1;

  
  [F,n,m,l] = read_hycom(fina,finb,'temp');
  F(F>hgg)=nan;
  T=F;
  
  [F,n,m,l] = read_hycom(fina,finb,'salin');
  F(F>hgg)=nan;
  S=F;

  [F,n,m,l] = read_hycom(fina,finb,'thknss');
  F(F>hgg)=nan;
  F=F./rg;
  dP=F;

  if cc==1
    Tav  = zeros(l,m,n);
    Sav  = zeros(l,m,n);
    dPav = zeros(l,m,n);
  end
  
  Tav=Tav+dP.*T;
  Sav=Sav+dP.*S;
  dPav=dPav+dP;
end;
Tav=Tav./dPav;
Sav=Sav./dPav;
dPav=dPav./cc;

%keyboard

f_chck=0;
if f_chck==1
  SCT=[];
  sub_plot_section(T,dP,'temp',SCT);
end


% convert layer thickness back to pressure
%dPav=dPav*rg;

return
