  function CMP = create_colormap(nint,c1,c2,sclr);
% CMP = create_colormap(nint,c1,c2)
% Create colormap using specified "target" colors in sclr
% Best if nint is divisible by # of target colors
%  example: sclr = [0,0,0; 0,1,0; 0,0,1; 1,0,0; 1,1,1];
%           5 target colors
%           c1=0; c2=10;
%           nint = 10  -
% nint = 240; % should be divisible by 12 - # of colors
%     c1=-0.3;
%     c2=0.6;
%  CMP=create_colormap2_3(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

nsb=size(sclr,1)-1;
ni=nint/nsb;

if ni~=round(ni)
  fprintf('  COLORMAP:  # of nint=%i, #colors=%i \n',nint,ni);
  fprintf('  COLORMAP:  nint/ncolors does not  give integer\n');
  fprintf('  COLORMAP: adjusting nint = %i\n',round(ni)*nsb);
  nint = round(ni)*nsb;
  ni = nint/nsb;
end

%keyboard
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
 
return
 