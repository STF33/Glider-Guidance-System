% Create land masks
% using old nest files:
function LMASK=sub_lmask(fina,finb);
fprintf('Land mask indices ...\n');

fld='temp';
[F,n,m,l] = read_hycom(fina,finb,fld);
ID=n;
JD=m;
IJDM=ID*JD;
a=squeeze(F(1,:,:));
a=reshape(a',IJDM,1);
lmask=find(a>1e20);
LMASK.temp=lmask;

[F,n,m,l] = read_hycom(fina,finb,'thknss');
ID=n;
JD=m;
IJDM=ID*JD;
a=squeeze(F(1,:,:));
a=reshape(a',IJDM,1);
lmask=find(a>1e20);
LMASK.thknss=lmask;

[F,n,m,l] = read_hycom(fina,finb,'u-vel.');
ID=n;
JD=m;
IJDM=ID*JD;
a=squeeze(F(1,:,:));
a=reshape(a',IJDM,1);
lmask=find(a>1e20);
if isempty(lmask),
  lmask=find(a==0);  % 0-masking for u,v is used
end
LMASK.uvel=lmask;

%keyboard

[F,n,m,l] = read_hycom(fina,finb,'v-vel.');
ID=n;
JD=m;
IJDM=ID*JD;
a=squeeze(F(1,:,:));
a=reshape(a',IJDM,1);
lmask=find(a>1e20);
if isempty(lmask),
  lmask=find(a==0);  % 0-masking for u,v is used
end
LMASK.vvel=lmask;



return
