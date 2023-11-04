% Prepare *b clim file for T or S *a 
% in HYCOM format 
% in climatology files on Z-levels
% see interp3D_nemo_glorys_hycom.m 
% to interpolate data
% 
% HYCOM climatology files are needed to create relax files
% with HYCOM hybrid layers
%
function sub_create_climb(pthoutp,fldnm,flnma,dnmb,IDM,JDM,ZZN);

fid = [];

sgm = 'sig2';  % pressure reference
nlrs = 75;     % depth layers in NEMO

%dnmb = datenum(2011,5,1);  % interpolation date
DV = datevec(dnmb);

IJDM = IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
% 
% Read 
fina = sprintf('%s%s',pthoutp,flnma);
fprintf('Reading: %s\n',fina);
fid=fopen(fina,'r');

flnmb = flnma(1:end-2);
foutb = sprintf('%s%s.b',pthoutp,flnmb);

%fida = fopen(fouta,'w');
fidb = fopen(foutb,'wt');

fprintf('Creating climatology %s\n',fidb);

% Write headers in *b file
stl = 'NEMO+GLORYS on IAS HYCOM-TSIS\n';
fprintf(fidb,stl);
switch(fldnm),
 case('temp');
  stl='Potential Temperature';
  bline = 'potential temperature: depth,range =';
 case('saln');
  stl='Salinity';
  bline = '             salinity: depth,range =';
end
fprintf(fidb,stl);
stl=' \n';
fprintf(fidb,stl);
fprintf(fidb,stl);
fprintf(fidb,stl);

stl=sprintf('i/jdm = %i %i\n',IDM,JDM);
fprintf(fidb,stl);

for ik=1:nlrs
  fprintf('Reading layer %i\n',ik);
  dmm=fread(fid,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid,npad,'float32','ieee-be');  % read npad   
  F = dmm;
 
  I=find(abs(F)<1e20);
  fmin = min(F(I));
  fmax = max(F(I));
%  fwrite(fida,F,'float32','ieee-be');
%  fwrite(fida,toto,'float32','ieee-be');  

% *b file:
  zz0=ZZN(ik);
  if (ik==1); zz0=0; end;
  aa2 = sub_parse_cline(bline,zz0,fmin,fmax);
  fprintf('===>  %s\n',aa2);
  fprintf(fidb,[aa2,'\n']);

end
fclose(fid);
%fclose(fida);
fclose(fidb);

fprintf(' Created %s\n',foutb);








