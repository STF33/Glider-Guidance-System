  function CMP = create_colormap2_1(nint,c1,c2);
% CMP = create_colormap2_3(nint,c1,c2);
% Create colormap: white-grey for the interval that is 
% out of interest, then purple - blue - yellow - red for main interval
% Best if nint is divisible by 12
%  example: 
% nint = 240; % should be divisible by 12 - # of colors
%     c1=-0.3;
%     c2=0.6;
%  CMP=create_colormap2_3(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
sclr=[1 1 1;
      0.3 0.3 0.3;
      0.5 0 0.5;
      1, 0, 1;
      0, 0, 0.4;
      0, 0, 1;  % 
      0, 1, 1;
      .2, 1, .7;
      1, 1, .5;
      1, 0.5, 0;
      1, 0.5, 0.5;
      1, 0, 0;
      0.3, 0, 0];

nsb=size(sclr,1)-1;
ni=nint/nsb;

if ni~=round(ni)
  fprintf('   create_colormap2_3: nint=%i, #colors=%i \n',nint,ni);
  fprintf('           nint/ncolors should give integer\n');
  fprintf('           Using  nint = %i\n',round(ni)*nsb);
  nint = round(ni)*nsb;
  ni=nint/nsb;
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
 
return
 