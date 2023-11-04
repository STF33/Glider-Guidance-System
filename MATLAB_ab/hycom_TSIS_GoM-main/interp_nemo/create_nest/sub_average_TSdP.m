function [Tav,Sav,dPav] = sub_average_TSdP(pthrlx0,dnmb0,Nday,hr1,hr2,dhr);
% 
% average T,S, dP fields on HYCOM grid
% Do Nday averaging
% from 3-hourly fields -> Nday- mean
% to avoid the problem of wierd isopycnals
% Nday = 0 - no averaging

bug - does not work when no averaging
Nday = 00000

%Nday = 3; % has to be an odd number!
dNd  = (Nday-1)/2;
if Nday==0,
  fprintf('Relax: No time averaging of relax fields\n');
  dNd = 0;
end

DV=datevec(dnmb0);
iyr=DV(1);
dj1=datenum(iyr,1,1);
iday=dnmb0-dj1+1;

rg = 9806;
hg = 2^100;
hgg= 1e20;

Tav=[];
Sav=[];
dPav=[];

%keyboard    

cc=0;
for dnmb=dnmb0-dNd:dnmb0+dNd
  DV=datevec(dnmb);
  iyr=DV(1);
  dj1=datenum(iyr,1,1);
  iday=dnmb-dj1+1;
  pthrlx  = sprintf('%s%4.4i/',pthrlx0,iyr);
%  pthrlx  = pthrlx0;

  for hr=hr1:dhr:hr2
    fina = sprintf('%srelax.%4.4i%3.3i%2.2i.a',pthrlx,iyr,iday,hr);
    finb = sprintf('%srelax.%4.4i%3.3i%2.2i.b',pthrlx,iyr,iday,hr);

    fprintf('    sub_aveage_TSdP: reading %2.2ihr\n',hr);

    if ~exist(fina,'file'), 
      fprintf('!!!    Missing Relax: %s\n',fina);
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

    Tav = Tav  + dP.*T;
    Sav = Sav  + dP.*S;
    dPav= dPav + dP;

    if Nday == 0, return; end; % no averaging, get records
    
  end;
end;
Tav = Tav./dPav;
Sav = Sav./dPav;
dPav= dPav./cc;

%keyboard

f_chck=0;
if f_chck==1
  SCT=[];
  f_lr=1;
  sttl='check'
  sub_plot_section(T,dP,'temp',SCT,sttl,f_lr);
end


% convert layer thickness back to pressure
%dPav=dPav*rg;

return
