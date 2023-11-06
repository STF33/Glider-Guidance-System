  function CMP = colormap_sclr1(nint,c1,c2);
% Create colormap: darkblue - turqouise white - yellow - brown

cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

%
sclr=[    0.2081                 0.1663                    0.5292
       0.00595714285714286       0.408614285714286         0.882842857142857
        0.0640571428571428       0.556985714285714         0.823957142857143
        0.0589714285714286       0.683757142857143         0.725385714285714
         0.6                     1                          0.8
	  1                      1                         1
                    0.9763                  0.9831                    0.0538
         0.999042857142857       0.765314285714286         0.216414285714286
         1                       0.6                       0.
         1                       0.4                       0.
         0.5                     0.2                       0.1];
    
nsb=size(sclr,1)-1;
ni=ceil(nint/nsb);
%keyboard

cmp=[];
for ik=1:nsb;
  cl1=sclr(ik,:);
  cl2=sclr(ik+1,:);
  clrM=mix_2colors(cl1,cl2,ni);
  icc=size(cmp,1);
  cmp(icc+1:icc+ni,:)=clrM;
end;

nsm = round(ni);
cmp = smooth_colormap(cmp,nsm);

CMP.colormap=cmp;
CMP.intervals=cnt;

return