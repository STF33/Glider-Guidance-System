function jmp=sub_mhd_jump(MHD,dlt0,mhd0);
% Find first occurence of large 
% mhd score
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
 
mhd=MHD;
for im=1:nn
  dm=abs(MHD(im)-m0)/m0;
  dM(im)=dm;
  jmp=im;
%  fprintf('im=%i, dm=%6.5f\n',im,dm);
%  if dm<=dlt0 & MHD(im)<mhd0   % large MHD - LCE shedding
  if MHD(im)<mhd0   % large MHD - LCE shedding
    m0=MHD(im);
  else
    break;
  end
end



return