function CMP = sub_conc_clrmpA(c1,c2);
% 
% colormap for tracer conc.
% nat-log scale
n=40;
A=ones(n,3);
cl1=colormap_gray(n);
for ik=1:n
  cl2(ik,:)=[0,0,1];
  cl3(ik,:)=[0,1,1];
  cl4(ik,:)=[0,0.8,0];
  cl5(ik,:)=[0.6,1,0];
  cl6(ik,:)=[1,0.5,0];
  cl7(ik,:)=[0.8,0,0];
end

cmp=[cl1;cl2;cl3;cl4;cl5;cl6;cl7];
cmp=smooth_colormap(cmp,25,2);
cmp(1,:)=[1 1 1];

nint=length(cmp);
%c1=-6;
%c2=1;
cnt=(c1:(c2-c1)/nint:c2);


CMP.colormap=cmp;
CMP.intervals=cnt;
CMP.c1=c1;
CMP.c2=c2;

return
