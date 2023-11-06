function Ssm = sub_fill_bottom_nans(Ssm);
% 2D array
% for plotting and interpolation
% fill bottom nans with closest values
% Fill nans
I=find(isnan(Ssm));
if isempty(I), return; end;

fprintf('Filling bottom nans ...\n');

[mm,nn]=size(Ssm);
if isnan(Ssm(1,1)),
  i1=min(find(~isnan(Ssm(1,:))));
  for kk=1:mm
    Ssm(kk,1:i1-1)=Ssm(kk,i1);
  end
end
if isnan(Ssm(end,1)),
  i1=max(find(~isnan(Ssm(1,:))));
  for kk=1:mm
    Ssm(kk,i1+1:end)=Ssm(kk,i1);
  end
end

s0=nanmean(nanmean(Ssm));
for ii=1:nn
  i1=min(find(isnan(Ssm(:,ii))));
  if ~isempty(i1) & i1>1
    Ssm(i1:end,ii)=Ssm(i1-1,ii);
  elseif ~isempty(i1) & i1==1 % iland 
    Ssm(:,ii)=Ssm(:,ii-1);
  end
end

return