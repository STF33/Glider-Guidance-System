    function [F,n,m,l] = read_hycom(fina,finb,fld,varargin);
%
%  function [F,n,m,l] = read_hycom(fina,finb,fld);
%  reads hycom binary archive files (model output), 
%  returns specified field 'fld'
%  and dimensions of the grid
%
%%  Dmitry Dukhovskoy, FSU, 2010
%  2017: added options for nlayer, n tracers
% if several tracers, options are:
% read tracer N: 'r_tracer',1
% read all tracers by default
% any variable can be read in 1 layer Nl:
%       'r_layer',1
%
% If all tracers are read, recommended to specify 
% only 1 layer to read 
%
% Examples: [F,n,m,ll] = read_hycom(arch1.a,arch1.b,'tracer','r_tracer',2)
%     returns tracer 2 in all v. layers
% 
%   [F,n,m,ll] = read_hycom(arch1.a,arch1.b,'tracer','r_layer',1)
%   reads all tracers in layer 1, returns as Ntracers : n : m matrix
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
  nTR=find(dI>1,1); % # tracers or field 
  nVlev = ll/nTR;   % actual # of v layers
end

fprintf('\nReading HYCOM output, field: %s\n',fld);
fprintf('%s\n',fina);
fprintf('Model dim: n(x)=%i, m(y)=%i, l(z)=%i \n',n,m,nVlev);
if nTR~=1
  fprintf('!!: Found %i variables %s per layer\n',nTR,fld);
end



% Find layers to read, if specified
% and tracers, if specified
lr1=0;
lr2=0;
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
      FLOC= dmm(Rtrc:nTR:end);
      if lr1==0 | lr2==0 % not to overide r_layer flag
        lr2 = length(FLOC);
        ll  = lr2;
      end
      
    end
    if strmatch(vfld,'r_layer'); % read 1 layer
      rLayer = varargin{k+1};
      lr1 = rLayer;
      lr2 = rLayer;
    end
  end
end

% If no layers specified - read all layers
if lr1==0 | lr2==0
  lr1 = 1;
  lr2 = ll;
end  


%keyboard

fid1=fopen(fina,'r');
F=[];
ccL=0;
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

l = size(F,1);

fclose(fid1);

return








