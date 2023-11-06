    function TDENS = read_targ_dens(fnm,ss,sep,fld);
%
% TDENS = read_targ_dens(fnm,ss,sep,fld);
% Read target density from HYCOM *b file, fnm
% ss is a piece of a string where to start reading targ densities 
% sep - separator character after which numeric values begin in the string
% fld  - field # of the density in the numeric values
%
%e.g.: relax_sal.b:
% PHC 3.0 (NOAA WOA-98 with improved Arctic) Climatology    
%Expt 06.0 nhybrd=32 nsigma=14 ds00= 0.50 dp00= 3.00 dp00x= 120.0 dp00f=1.180  
%Layered averages w.r.t. Sigma-2,  levtop=2                                     
%Salinity                                                                       
%i/jdm = 1600 2520                                                              
% sal: month,layer,dens,range = 01 01 28.100 3.7587984 38.609257   % 
% sal: month,layer,dens,range = 01 02 28.900 3.7587984 38.609257    
%
% ...
%%
% fnm = 'relax_sal.b'
% TDENS=read_targ_dens(fnm,'sal','=',3); 3rd value is density
%

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
    iq=strfind(pp,sep);
    b=pp(iq+1:end);
    S=sscanf(b,'%f');

% In some *b files, several fields are written, and target densities
% loop over and over, need to check this:
    if kl>1
      if S(fld)==TDENS(1); break; end
    end

    kl=kl+1;
    TDENS(kl,1)=S(fld);
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



