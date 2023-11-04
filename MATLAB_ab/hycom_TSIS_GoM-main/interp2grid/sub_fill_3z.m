function F = sub_fill_3z(F);
%
% 3D array
% Fill all nans below the 1st upper layer
% with values for vertical interpolation
% 1 layer must not have nans
%
fprintf('Filling nans 3z\n');
a=squeeze(F(1,:,:));
I1=find(isnan(a));
if ~isempty(I1)
  fprintf('sub_fill_3z: 1st lyaer cannot have nans\n');
  error(' quitting ');
end

[kk,mm,nn]=size(F);
for ik=2:kk
  In=find(isnan(F(ik,:,:)));
  F(ik,In)=F(ik-1,In);
end


return