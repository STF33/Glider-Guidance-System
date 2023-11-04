    function TDENS = read_targ_dens_blkdat(fnm);
%
% TDENS = read_targ_dens_blkdat(fnm,ss,sep,fld);
% Read target density from HYCOM blkdat.input file, fnm
% fnm - full path + blkdat.input file name
% 
ss='''sigma ';
fid=fopen(fnm,'r');
%keyboard
pp=' ';
kl=0;
TDENS=[];
while ischar(pp)
  pp = fgetl(fid);
  if ~ischar(pp); break; end;

  I=strfind(pp,ss);
  if ~isempty(I)
    b=pp(1:I-2);
    S=sscanf(b,'%f');

    kl=kl+1;
    TDENS(kl,1)=S;
    if kl>1000; 
      fprintf('Endless loop Reading blkdat.input ...\n');
      error('Check formating in blkdat');
    end
    
  end
%keyboard
end;

if kl==0, 
  fprintf(' No t. dens in %s\n',fnm);
  error('Could not find target densities ...'); 
end;

fprintf('###  Target densities read, # t. densities  =  %i\n',kl);
fprintf('###    min T.dens=%f, max =%f\n',min(TDENS),max(TDENS));

fclose(fid);



