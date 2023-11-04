  function CMP = create_colormap5(nint,c1,c2);
% CMP = create_colormap5(nint,c1,c2);
% similar to colormap used in Fieg et al., Oc.Dyn., 2010
% figure shows T in Fram Strait section
% Create colormap: purple - green - yell - red
% 
%  example: 
% nint=200;
%     c1=0;
%     c2=20;
%  CMP=create_colormap5(nint,c1,c2);
%  cmp=CMP.colormap;
%  cnt=CMP.intervals;


%nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

% Create colormap by mixing 2 colors
% within 1 interval
%
sclr=[76  0  153;
      200 146 235;
      120 50 255;
      51 51  255;  % 
      0 102  51;
      51 255 51;
      255 255 0;
      255 178 102;
      255 128 0;
      255 51  51;
      153 0  0];
%sclr=[76  0  153;
%      153 51 255
%      204 153 255;
%      51 51  255;  % 
%      0 102  51;
%      51 255 51;
%      255 255 0;
%      255 178 102;
%      255 128 0;
%      255 51  51;
%      153 0  0];
sclr = sclr/255;
nsb=size(sclr,1)-1;
ni=round(nint/nsb);
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
