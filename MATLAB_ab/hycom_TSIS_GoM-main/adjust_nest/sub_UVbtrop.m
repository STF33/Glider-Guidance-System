function [Ubt,Vbt] = sub_UVbtrop(UU,VV,dH);
% U and V  and dH (layer thicknesses)
% are on the C grid
% Make sure to write u_brtrop, v_bartrop max/min in *b file
% otherwise (if both are 0s), 
% rd_archive (forfun.f) makes whole field=0
% Ubrtop and Vbrtrop are used at the boundaries (as well as surf. elev)
% NOTE: Barotrop. U,V are not nan's at the OBs (unlike baroclinic 3D)

[kk,mm,nn]=size(dH);

clear dHu
for k=1:kk
  dh1=squeeze(dH(k,:,1:nn-1));
  dh2=squeeze(dH(k,:,2:nn));
  dh=0.5*(dh1+dh2);
  dHu(k,:,:)=[zeros(mm,1),dh];
  dh1=squeeze(dH(k,1:mm-1,:));
  dh2=squeeze(dH(k,2:mm,:));
  dh=0.5*(dh1+dh2);
  dHv(k,:,:)=[zeros(1,nn);dh];
end
dmm=UU.*dHu;
Hu=squeeze(sum(dHu,1));
Hu(Hu==0)=nan;
Ubt=squeeze(sum(dmm,1))./Hu;

dmm=VV.*dHv;
Hv=squeeze(sum(dHv,1));
Hv(Hv==0)=nan;
Vbt=squeeze(sum(dmm,1))./Hv;

% 
% Ubarotropic & Vbtrop are not 0 at the OBs
Ubt(1,:)=Ubt(2,:); %S bndry
Vbt(1,:)=Vbt(2,:);
Ubt(:,end)=Ubt(:,end-1); %E bndry
Vbt(:,end)=Vbt(:,end-1);
Ubt(end,:)=Ubt(end-1,:); %E bndry
Vbt(end,:)=Vbt(end-1,:);



%keyboard
return