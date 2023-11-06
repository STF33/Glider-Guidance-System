    function [Uav,Vav,Eav,ZZ] = sub_avrgUV(pthin,dnmb0,Nday,dP,hr1,hr2,dhr);
% TIme average U,V, SSH fields fields
% from nest files on z-levels
% THe output will be interpolated
% into HYCOM layers in final nest files

dNd  = (Nday-1)/2;
if Nday==0,
  fprintf('UV: No time averaging of relax fields');
  dNd = 0;
end

DV=datevec(dnmb0);
iyr=DV(1);
dj1=datenum(iyr,1,1);
iday=dnmb0-dj1+1;

rg = 9806;
hg = 2^100;
hgg= 1e20;

Uav=[];
Vav=[];
Eav = [];
ZZ = [];

cc=0;
for dnmb=dnmb0-dNd:dnmb0+dNd
  DV    = datevec(dnmb);
  iyr   = DV(1);
  dj1   = datenum(iyr,1,1);
  iday  = dnmb-dj1+1;
%  pthin = sprintf('%s%4.4i/',pthin0,iyr);

  for hr=hr1:dhr:hr2
    fmatUV = sprintf('%seuv_nest_%4.4i%3.3i%2.2i.mat',...
		   pthin,iyr,iday,hr);

    if ~exist(fmatUV,'file'), 
      fprintf('UV:  Missing: %s\n',fmatUV);
      Uav = -1;
      Vav = -1;
      Eav = -1;
      continue; 
    end;

   fprintf('    sub_avrgUV: Reading %4.4i-%3.3i-%2.2ihr\n',iyr,iday,hr);

    cc=cc+1;
    load(fmatUV);
    
    if cc==1
      [l,m,n]=size(u_nest);
      Uav  = zeros(l,m,n);
      Vav  = zeros(l,m,n);
      Eav  = zeros(m,n);
      ZZ   = -double(depth);
    end

    Uav=Uav+double(u_nest);
    Vav=Vav+double(v_nest);
    Eav=Eav+double(e_nest);
    
  end;
end;
Uav=Uav/cc;
Vav=Vav/cc;
Eav=Eav/cc;
if Uav == -1, return; end; % no files

% Correct Uav, Vav at the OB if needed
% in the nest from 022 experiment, the last row/column u
% were half of what it should be - error in creating nest
% files
u1=squeeze(Vav(1,m-1,:));
u2=squeeze(Vav(1,m-2,:));
i0=find(u1==max(u1));
rr=abs(u1(i0)/u2(i0));
if rr<0.8,
  fprintf(' =====   Nest files: OB row needs correction ===\n');
  fprintf('row1 u = %8.3f, row2  u = %8.3f\n',u1(i0),u2(i0));

  [UUc,VVc] = sub_correctUV_OB(Uav,Vav,dP);
%keyboard
  Uav = UUc;
  Vav = VVc;
else
  fprintf(' ==  Nest files: OB rows no correction needed\n');
end


f_chck=0;
if f_chck==1
  dp=diff(abs(ZZ));
  dp(1)=dp(1)-0.5;
  dp=[0.5;dp];
  
  [dm1,dPr,dm2]=meshgrid([1:m],dp,[1:n]);
  
  SCT=[];
  f_lr=0;
  sttl=sprintf('Time-avrg %i d',Nday);
  sub_plot_section(Uav,dPr,'U',SCT,sttl,f_lr);
end



return