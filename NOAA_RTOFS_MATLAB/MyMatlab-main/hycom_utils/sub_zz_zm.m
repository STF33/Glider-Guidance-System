function [ZM,ZZ] = sub_zz_zm(fina, finb,HH,varargin);
% [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',[0,1]); 
% derive layer depths: zz - interface depths
%   zm - mid cell depths
% f_btm=1 - makes ZZ=nan below bottom
% if 'thknss',dH is provided as optional
% input field, then dH is not 
% read in from fina files again

% Dmitry Dukhovskoy, COAPS FSU, Nov. 2017
%              May 2018 - adding f_btm option
%
nV = length(varargin);
f_btm=0;
%keyboard
for k=1:nV
  aa=varargin{k};
  aa=lower(aa);
  if strncmp(aa,'f_btm',2);
    f_btm=varargin{k+1}; % f_btm=1 - make ZZ=nan below bottom
  elseif strncmp(aa,'thknss',6);
    F=varargin{k+1};
  end;  
end

fprintf('Creating ZZ, ZM from HYCOM dH, f_btm=%i...\n',f_btm);

rg = 9806;

if ~exist('F','var')
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e10)=nan;
  F(F<0.01)=0;
  F=F/rg;
end
[l,m,n]=size(F);
I = find(HH<0);

fprintf('Deriving ZZ, ZM \n');
clear ZZ
ZZ = F*0;
ZZ(l+1,:,:)=ZZ(l,:,:);

for kk=2:l+1
  ZZ(kk,:,:)=ZZ(kk-1,:,:)-F(kk-1,:,:);
  if f_btm==1
    dmm=squeeze(F(kk-1,:,:));
    Ib=find(dmm<1e-3);
    if ~isempty(Ib),
      ZZ(kk,Ib)=nan;
    end
  end
end

% Depths of the middle of the cells:
ZM = F*0;
for kk=1:l
  ZM(kk,:,:)=0.5*(ZZ(kk+1,:,:)+ZZ(kk,:,:));
end

%keyboard

return
