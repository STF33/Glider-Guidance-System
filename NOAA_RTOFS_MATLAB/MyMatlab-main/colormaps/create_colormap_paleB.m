  function CMP = create_colormap_paleB(nint,c1,c2);
% CMP = create_colormap_paleB(nint,c1,c2);
% Create colormap: pale light-blue-blue 
%  example: 
%     c1=-0.3;
%     c2=0.6;
%  CMP=create_colormap1(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[0.9 0.6 0.9;
      0.86  0.86  1;
      1 0.8  0.8;
      1 0.9 0.95;
      0.8 1 0.8;
      0.92 1 0.92;
      1  1 0.6;
      0.8 1 1;
      1 0.4 1;
      0.9 0.7 0.7];
nsb=size(sclr,1)-1;
ni=ceil(nint/nsb);



cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;
  
nint=length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

CMP.colormap=cmp;
CMP.intervals=cnt;

return
