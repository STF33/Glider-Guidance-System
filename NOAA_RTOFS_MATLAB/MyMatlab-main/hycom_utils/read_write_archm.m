% Test for checking reading/writing archm HYCOM output
% Note that Alan uses direct access files
% they require all records be the same length and this record 
% length needs to be specified in OPEN fortran statment
% record length is machine-dependent, use inquire to 
% figure this out (see example in fortran/hycom/ codes)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

f_read=2; %=1 - read HYCOM by skipping record length
          %=2 - read HYCOM consecutievly reading all records
expt='060';
yr=2004;
iday=8;
pthbin = sprintf(...
    '/Net/tholia/ddmitry/hycom/ARCc0.08/060/%4.4i/bin_outp/',yr);
pthout = sprintf(...
    '/Net/tholia/ddmitry/hycom/ARCc0.08/060/%4.4i/',yr);
fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
fouta= sprintf('%s%s_archm.%4.4i_%3.3i_12test.a',pthout,expt,yr,iday);

fprintf('Input: %s\n',fina);
fprintf('Output: %s\n',fouta);


% Check *b file
% count # of all fields in the *a
f1 = fopen(finb,'r');  % read I,J from *.b
for nl=1:7
  aa=fgetl(f1);
%  disp(aa);
end

aa=fgetl(f1);
%dmm=aa(2:8);
%ID=str2num(dmm);
 [ID, c1,c2,c3,c4,c5,c6] = strread(aa,'%d%s%s%s%s%s%s');
aa=fgetl(f1);
%dmm=aa(2:8);
%JD=str2num(dmm);
 [JD, c1,c2,c3,c4,c5,c6] = strread(aa,'%d%s%s%s%s%s%s');

disp(['Grid I=',num2str(ID),' J=',num2str(JD)]);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

aa=fgetl(f1);

a=1;
cntr=0;
nrec=0;
FLOC=[];
while (a),  
  aa=fgetl(f1);
  fprintf('%s\n',aa);
  if ~ischar(aa); break; end
  cntr=cntr+1;
  nrec=nrec+1;
  FLOC(nrec)=cntr;
end;
fclose(f1);
l=length(FLOC);

fid1=fopen(fina,'r');
fout1=fopen(fouta,'w');
if f_read==1
  fprintf('Read HYCOM by skipping record length ...\n');
  F=[];
  for ii=1:l
    fprintf('Writing l=%i out of %i\n',ii,l);
    frewind(fid1);
    k0=FLOC(ii)-1;
    stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
    dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
    dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 
    fwrite(fout1,dmm,'float32','ieee-be');
    fwrite(fout1,dm1,'float32','ieee-be'); % padding
  end;
elseif f_read==2
  fprintf('Read HYCOM by sequential access ...\n');
  frewind(fid1);
  ii=0;
  while ~feof(fid1)
    dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
    dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 
    if isempty(dm1),
      fprintf('EOF fid1\n');
      break
    end
    ii=ii+1;
    fprintf('Writing l=%i out of %i\n',ii,l);
%    keyboard
    fwrite(fout1,dmm,'float32','ieee-be');
    fwrite(fout1,dm1,'float32','ieee-be'); % padding

  end;
end


fclose(fid1);
fclose(fout1);


