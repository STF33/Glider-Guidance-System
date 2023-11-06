  function CMP = colormap_WG(nint,c1,c2);
% CMP = colormap_WG(nint,c1,c2);
% Create colormap: white - green
%  example: 
%     c1=-0.3;
%     c2=0.6;
%  CMP=create_colormap1(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[1, 1., 1;
      .8, 1, .9;  % 
      .7, 1, 0.8;
      .5 0.8 0.6;
      0.3 0.5 0.4];
nsb=size(sclr,1)-1;
ni=nint/nsb;



cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;
  
cmp=smooth_colormap(cmp,round(0.05*nint));

CMP.colormap=cmp;
CMP.intervals=cnt;
