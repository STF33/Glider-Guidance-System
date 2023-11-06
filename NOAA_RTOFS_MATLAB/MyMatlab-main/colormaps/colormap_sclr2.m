  function CMP = colormap_sclr2(nint,c1,c2);
% Create colormap: dark purple - blue - turqouise - yellow - red - brown

% Create colormap by mixing 2 colors
% within 1 interval
%
% Dmitry Dukhovskoy, COAPS FSU, 
% Jan 2014
% Modified 2016
%
sclr=[    0.4                    0                         0.4
	  0.9                    0                         0.9
	  0.5                    0                         1
          0.2081                 0.1663                    0.5292
       0.00595714285714286       0.408614285714286         0.882842857142857
        0.0640571428571428       0.556985714285714         0.823957142857143
        0.0589714285714286       0.683757142857143         0.725385714285714
         0.395257142857143                  0.7459         0.524442857142857
         0.752485714285714                  0.7384         0.376814285714286
         0.999042857142857       0.765314285714286         0.216414285714286
                    0.9763                  0.9831                    0.0538
         0.9                      0.5                       0.
	 0.8                      0.3                       0
	 0.7                     0                         0
	 0.6                     0.4                       0.3
         0.5                     0.1                       0.1];
    
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
  
cmp = smooth_colormap(cmp,ni);
nint = length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals


CMP.colormap=cmp;
CMP.intervals=cnt;

return