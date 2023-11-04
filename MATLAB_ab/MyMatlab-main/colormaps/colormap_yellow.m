  function cmp = colormap_yellow(nint);
% yellow->ornage
cmp=[];
cl1=[1 1 0.9];
cl2=[1 1 0];

n1=round(nint/2);
CL=mix_2colors(cl1,cl2,n1);

cl1=CL(end,:);
cl2=[1 0.6 0];
n2=nint-n1+1;
CL2=mix_2colors(cl1,cl2,n2);

cmp=[CL;CL2(2:end,:)];

return;

