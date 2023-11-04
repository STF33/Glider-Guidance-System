function [U,V] = sub_adjust_Ubtrop(UU,VV,dH,adj_btrop);
% Adjust U barotropic in the nest files
% to match Yucatan flow
% After adjusting the barotropic field recalculate
% new total U,V ("u-vel.", "v-vel.") using 
% baroclinic u,v (need to calculate those first)
% In the nest files:
% u_vel=u_btrop+u_bclinic
%
% So: (1) adjust u_btrop  -> u_btrop_new
%     (2) adjust u total:
%        get b/clinic u,v: u_bclinic=u_vel-u_btrop_old
%        u_vel(new)=u_bclinic+u_btrop_new
%
%
[kk,mm,nn]=size(dH);
fprintf('Barotropic U adjustment: %6.2f\n',adj_btrop);
 
% Calculate U, V barotropic:
[Ubt,Vbt] = sub_UVbtrop(UU,VV,dH);
Ubadj=Ubt*adj_btrop;
Vbadj=Vbt*adj_btrop;

% Calc. anomalies:
U=UU*0;
V=VV*0;
for k=1:kk
  u=squeeze(UU(k,:,:));
  v=squeeze(VV(k,:,:));
  ubcl = u-Ubt;
  vbcl = v-Vbt;
  utot = ubcl+Ubadj;
  vtot = vbcl+Vbadj;
  U(k,:,:)=utot;
  V(k,:,:)=vtot;
end

return