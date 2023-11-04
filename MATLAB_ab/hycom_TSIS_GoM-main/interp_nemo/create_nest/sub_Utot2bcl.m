% Subtract Ubtrop from Utot to create Ubcl
% for archv
%
function [Ubcl,Vbcl] = sub_Utot2bcl(UU,VV,Ubt,Vbt);
fprintf('Removing barotopic: U/V total ---> U/V baroclinic\n');

[kk,mm,nn]=size(UU);

Ubcl = UU*0;
Vbcl = VV*0;
for k=1:kk
  a=squeeze(UU(k,:,:));
  b=squeeze(VV(k,:,:));
%  Iu=find(isnan(a));
%  Iv=find(isnan(b));
  a = a-Ubt;
  b = b-Vbt;
%  a(Iu) = nan;
%  b(Iv) = nan;
  Ubcl(k,:,:)=a;
  Vbcl(k,:,:)=b;
end

%keyboard
return
