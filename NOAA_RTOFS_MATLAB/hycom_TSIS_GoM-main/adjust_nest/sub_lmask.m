function LMASK=sub_lmask(iyr,iday);
% Create land masks
% using old nest files:
fprintf('Land mask indices ...\n');

pthnest022 = '/Net/kronos/ddmitry/hycom/TSIS/nest_files/';
AB='a';
fina=sprintf('%sarchv.%4.4i_%3.3i_00.a',pthnest022,iyr,iday);
finb=sprintf('%sarchv.%4.4i_%3.3i_00.b',pthnest022,iyr,iday);
fprintf('%s\n',fina);

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
LMASK.uvel=lmask;

[F,n,m,l] = read_hycom(fina,finb,'v-vel.');
ID=n;
JD=m;
IJDM=ID*JD;
a=squeeze(F(1,:,:));
a=reshape(a',IJDM,1);
lmask=find(a>1e20);
LMASK.vvel=lmask;



return
