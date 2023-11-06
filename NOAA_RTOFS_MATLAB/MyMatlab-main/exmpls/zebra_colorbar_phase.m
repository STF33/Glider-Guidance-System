addpath /Net/Movies0/ddmitry/MyMatlab

close all
clear

h=ones(50,50);

for j=1:50
for i=1:50
  h(i,j)=(i)*360/50*sqrt(j);
end
end
cf=max(max(h))/360;
h=h/cf;

c1=0;
c2=360;
nint=5;
cct=[c1:(c2-c1)/nint:c2];
cmp=['M';'B';'G';'Y';'O';'R'];
%cmp=[ ];
nsint2=4;

jjs = clrmp_zebra (cct,cmp,nsint2);
pcolor(h);
shading flat
colormap(jet(length(cct)));
caxis([c1 c2]);
colormap(jjs)
% =======================
% Colorbar
% =======================
      hght=[];
      lngth=[];
      mint=1;
      fsz=12;
      bxc='k';
      posc=[0.3 0.025 0.6 0.05];
%	    nsint2=mint;
      [az,axc] = pcolorbar_zebra_horiz(jjs,cct,hght,lngth,mint,fsz,...
                             bxc,posc,nsint2);












