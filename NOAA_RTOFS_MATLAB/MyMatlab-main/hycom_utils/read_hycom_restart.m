%  function F = read_hycom_restart(fina,finb,fld,ID,JD);
%  reads hycom binary restart files (model output), 
%  returns specified field 'fld' 
% if fields is saved at several time levels
% on ly the 1st is read in
%  ID,JD - J and I dimensions
%%  Dmitry Dukhovskoy, FSU, 2016
%%  April 2017 - added option to read N tracers, in N layer(s)
%
% if several tracers, options are:
% read tracer N: 'r_tracer',1
% read all tracers by default
% any variable can be read in 1 layer Nl:
%       'r_layer',1
%
% If all tracers are read, recommended to specify 
% only 1 layer to read 
%
% Only first time level is read (there are 2 in the restart)
%
% Examples: [F,n,m,l] = read_hycom_restart(arch1.a,arch1.b,'tracer','r_tracer',2)
%     returns tracer 2 in all v. layers
% 
%   [F,n,m,l] = read_hycom_restart(arch1.a,arch1.b,'tracer','r_layer',1)
%   reads all tracers in layer 1, returns as Ntracers : n : m matrix
%
function [F,n,m,l] = read_hycom_restart(fina,finb,fld,ID,JD,varargin);

%
f1 = fopen(finb,'r');  % read I,J from *.b
for nl=1:2
  aa=fgetl(f1);
%  disp(aa);
end

IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


% Find locations of required field:
a=1;
cntr=0;
nrec=0;
nf=length(fld);
FLOC=[];
tstp0=0;
while (a),  
  aa=fgetl(f1);
  if ~ischar(aa); break; end
  IL=strfind(aa,' ');
  Lstr=IL(1)-1; % str. length in *b
  cntr=cntr+1;
  if Lstr~=nf, continue; end;
  isb=strfind(aa,'=');
  if isempty(isb),
    error('*b file: Check formatting, = symb in file');
  end
  bb=aa(isb+1:end);
  [lr,tstp,bmin,bmax]=strread(bb,'%d%d%f%f');
  if tstp0==0
    tstp0=tstp;
  end
  if tstp~=tstp0
    continue; % next time level
  end
  
  if strncmp(aa,fld,nf);
    nrec=nrec+1;
    FLOC(nrec)=cntr;
  end
end;

fclose(f1);

if isempty(FLOC),
  error('read_hycom:: could not find %s in %s',fld,finb);
end


n=ID;
m=JD;
ll=length(FLOC);
lr1 = 1;
lr2 = ll;
% If fld = tracer and # of tracers >1
% need to distinguish # of layers 
% vs # of tracers*layers
%keyboard
%if strmatch(fld,'tracer')
if ll==1
  nVlev = 1;
  nTR = 1; 
else
  dI=diff(FLOC);
  nTR=length(find(dI>1))+1; % # tracers or field
  nVlev = ll/nTR;   % actual # of v layers
end


fprintf('Reading HYCOM restart, field: %s\n',fld);
fprintf('Output dim: n(x)=%i, m(y)=%i, l(z or time)=%i, \n',n,m,nVlev);
if nTR~=1
  fprintf('!!: Found %i variables %s per layer\n',nTR,fld);
end

% Note that l - can be vertical layers or/and time steps, 1,2, 3


% Find layers to read, if specified
% and tracers, if specified
nV = length(varargin); 
if nV>0
  for k=1:nV
    vfld = varargin{k};
    if strmatch(vfld,'r_tracer')
      Rtrc = varargin{k+1}; % tracer to read
      if nTR<Rtrc
	error('# of saved tracers %i < # tracer to read %i',...
	      nTR, Rtrc);
      end
      dmm = FLOC;
      FLOC= dmm((Rtrc-1)*nVlev+1:Rtrc*nVlev);
      lr2 = length(FLOC);
      ll  = lr2;
    end
    if strmatch(vfld,'r_layer'); % read 1 layer
      rLayer = varargin{k+1};
      lr1 = rLayer;
      lr2 = rLayer;
    end
  end
end



fid1=fopen(fina,'r');
F=[];
ccL=0;
%keyboard
for ii=lr1:lr2
  frewind(fid1);
  k0=FLOC(ii)-1;
  stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
  dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 
%  fprintf('ii=%i, dmm=%ix%i\n',ii,size(dmm));
%keyboard
  dmm=reshape(dmm,ID,JD);
  ccL=ccL+1;
  F(ccL,:,:)=dmm';
end;
%keyboard

if ll==0
  fprintf('!!!!    %s is not found in the output  !!!\n',fld);
  fprintf('!!!!    Check field spelling in *.b file  !!!\n');
end

fclose(fid1);

[l,m,n] = size(F);
F = squeeze(F);

return





