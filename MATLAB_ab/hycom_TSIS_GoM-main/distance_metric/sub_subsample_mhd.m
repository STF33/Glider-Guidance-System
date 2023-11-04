function mhd=sub_subsample_mhd(MHD,dlt0);
% Exclude mhd values exceeding threshold value
%
nn=length(MHD);

m0=MHD(1);
if m0==0
  ip=min(find(MHD<0.2 & MHD>0));
  m0=MHD(ip);
  MHD(1)=m0;
end
  
if m0>0.2
  ip=min(find(MHD<0.2));
  m0=MHD(ip);
end

cc=0;
IX=[];
mhd=MHD;
for im=1:nn
  dm=abs(MHD(im)-m0)/m0;
  dM(im)=dm;
%  if dm<=dlt0 & MHD(im)<1.5  % large MHD - LCE shedding
  if dm<=dlt0   % large MHD - LCE shedding
    m0=MHD(im);
    cc=cc+1;  
    IX(cc)=im;
  else
    mhd(im)=nan;
  end
end

%plot(MHD); 
%hold on;
%plot(IX,MHD(IX),'r.');


return