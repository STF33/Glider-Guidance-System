  function CMP = create_colormap1(nint,c1,c2);
% CMP = create_colormap1(nint,c1,c2);
% Create colormap: white - dark-gray - l. purple - d. purple -  ...
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
sclr=[1, 1, 1;
      .5, .5, .5;
      1,  0.4, 1;
      0.4, 0, 0.4;  % 
      0, 0.9, 1;
      0, 0, 0.8;
      1, 1, .5;
      0, .7, 0;
      .8, 1, .8;
      .7, 0, 0;
      1, .6, .6];

nsb=size(sclr,1)-1;
ni=nint/nsb;

if ni~=round(ni)
  fprintf('ERROR: nint=%i, #colors=%i \n',nint,ni);
  fprintf('ERROR: nint/ncolors should give integer\n');
  error('Stopping colormap2');
end


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
