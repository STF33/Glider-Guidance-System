  function CMP = create_colormap_bath(nint,c1,c2);
% CMP = create_colormap_bath(nint,c1,c2)
% Creates colormap for bathymetry: lightgreen - dark blue
% 
%  example: 
% nint=200;
%     c1=0;
%     c2=20;
%  CMP=create_colormap_bath(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[0.95 1 0.95;
      0.8  1 0.85;
      0.1  0.9 0.7;
      0.   0.7 0.6;
      0.   0.6 0.8;
      0    0.4 0.8;
      0.   0.2  0.5];

%sclr = sclr/255;
nsb=size(sclr,1)-1;
ni=round(nint/nsb);

cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;
  
CMP.colormap=cmp;
CMP.intervals=cnt;
