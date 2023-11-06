    function [F,n,m,l] = read_hycom_arche(fina,finb,fld);
%
%  function [F,n,m,l] = read_hycom_arche(fina,finb,fld);
%  reads hycom binary archive files surface fields (model output), 
%
% Read surface fluxes exchange btw ocean/sea ice
% output arche files:
%sst      =   sst, deg. C
%sss      =   sss
%ssu      =   ocean current, m/s
%ssv      =   
%ssh      =   sea surf. height, m
%ssfi     =   ocean heat flux to sea ice (dwnwrd), W/m2
%mlt      =   ocean mixed layer thickn, m
%sic      =   sea ice conc
%sitxdown =   downard sea ice tau x to ocean
%sitydown =   
%siqs     =   solar heat flux through ice to ocean, W/m2
%sifh     =   ice freezing/melting heat flux (-1 from CICE), W/m2
%sifs     =   ice freez/melt salt flux (kg / m2 s)
%sifw     =   ice net water flux, dwnrd (kg /m2 s)
%sit      =   sea ice T, (deg C)
%sih      =   sea ice thickness, m
%siu      =   sea ice x vel., m/s
%siv      =   
%surtx    =  surface wind flux (with no ice) 
%surty    =   
%sflice   =   ice salt flux
%

% --------------------
% Get grid dimensions
%
% Formatted structure of *.b is presumed:
% after idm and jdm lines follows a header line
% then the data
% --------------------
f1 = fopen(finb,'r');  % read I,J from *.b
nl0=0;
while nl0<100
  nl0=nl0+1;
  aa=fgetl(f1);
  ii0=findstr(aa,'idm');
  if ~isempty(ii0), break; end
end
if isempty(ii0);
  fclose(finb);
  error('No idm found: Reading hycom %s\n',finb);
end
nl0=nl0-1;
frewind(f1);

for nl=1:nl0
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

%disp(['Grid I=',num2str(ID),' J=',num2str(JD)]);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

aa=fgetl(f1);

% Find locations of required field:
a=1;
cntr=0;
nrec=0;
nf=length(fld);
FLOC=[];
while (a),  
  aa=fgetl(f1);
  if ~ischar(aa); break; end
  cntr=cntr+1;
  if strncmp(aa,fld,nf);
    FLOC=cntr;
  end
end;

fclose(f1);

if isempty(FLOC),
  error('read_hycom:: could not find %s in %s',fld,finb);
end


n=ID;
m=JD;
ll=length(FLOC);

fprintf('\nReading HYCOM arche, field: %s\n',fld);
fprintf('%s\n',fina);
fprintf('Model dim: n(x)=%i, m(y)=%i \n',n,m);

lr1 = 1;

%keyboard

fid1=fopen(fina,'r');
F=[];
ccL=0;

frewind(fid1);
k0=FLOC-1;
stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 
%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
%keyboard
dmm=reshape(dmm,ID,JD);
ccL=ccL+1;
F(ccL,:,:)=dmm';
%keyboard

if ll==0
  fprintf('!!!!    %s is not found in the output  !!!\n',fld);
  fprintf('!!!!    Check field spelling in *.b file  !!!\n');
end

l = size(F,1);

fclose(fid1);

return








