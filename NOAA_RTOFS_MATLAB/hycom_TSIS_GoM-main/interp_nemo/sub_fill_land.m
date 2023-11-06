% Fill land points - Nan in NEMO fields
% with closest ocean values
% for interpolation
function AA = sub_fill_land(A0);
fprintf('Filling land ...\n');

[a1,a2]=size(A0);

I=find(isnan(A0(1,:)));
A0(1,I)=A0(2,I);
I=find(isnan(A0(a1,:)));
A0(a1,I)=A0(a1-1,I);
AA = A0;

In=find(~isnan(A0));
if isempty(In); % all nans
  AA=zeros(a1,a2)-999.;	
  return;
end

for kk=1:a1
  dmm = A0(kk,:);
%
% First check if the whole row is nan
% Typically this is at the boundaries
% or deep layers
  Igd = find(~isnan(dmm));
  if isempty(Igd)
    dmm(1:end)=nanmean(nanmean(A0));
  end;

  In=find(isnan(dmm));
  Igd = find(~isnan(dmm)); 

  if isempty(In)  % no land
    AA(kk,:)=dmm;
    continue;
  end

  In=dmm*0;
  In(isnan(dmm))=100;

  dIn=abs(diff(In));
  Jn = find(dIn>1);
%keyboard

  if In(1)>0, 
    Jn=[0,Jn];
  end
  if In(end)>0
    Jn=[Jn,a2];
  end
%
% Nan segments:
%  Jn=[0,Jn];
% Also in case the row ends with NaN:
%  if isnan(dmm(end)) & Jn(end)~=length(In)
%    Jn=[Jn,length(In)];
%  end

  for jj=1:2:length(Jn)-1
    j1=Jn(jj)+1;  % start of land
    j2=Jn(jj+1);  % end of land
    imm = [j1:j2];
    if j1>1,
%      dmm(j1)=dmm(j1-1);
      d1=dmm(j1-1);
    else             % left wall
%      fprintf('* kk=%i\n',kk);
%      dmm(j1)=dmm(j2+1);
      d1=dmm(j2+1);
    end
    if j2<a2
      d2=dmm(j2+1);
    else    % right wall
      d2=d1;
    end
 
    if length(imm)>1    
% 1D interpolation:
      idx = [j1:j2];
      phi1 = j2/(j2-j1)-1/(j2-j1)*idx;
      phi2 = -j1/(j2-j1)+1/(j2-j1)*idx;
      ll1 = d1*phi1+d2*phi2;      
%      dmm(j1+1:j2)=dmm(j1);
      dmm(j1:j2)=ll1;
    else
      dmm(j1)=0.5*(d1+d2);
    end

  end
%keyboard
  AA(kk,:) = dmm;

end

%keyboard

return
