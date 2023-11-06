    function [F,n,m] = read_pang(fina,finb);
%
%  function [F,n,m] = read_pang(fina,finb);
%  fina, finb - regional grids
%
%  Reads angles to rotate U,V from North-East orientation
% back to X (right) - Y (up) orientation on a grid
% or vice versa
%
% See Alan's email:
% The native model velocity variables are x-wards and y-wards (your #2), 
% but as part of the interpolation to fixed vertical levels we rotate them 
% to eastwards and northwards (your #1).  
% These are the same south of 47N, 
% where the grid is recti-linear, 
% but not north of 47N, where the grid is curvi-linear 
% in the Arctic bi-polar patch part of the tri-pole grid.
%
% If you need velocities North of 47N you have two options:
%
%  a) Use the GLBu0.08 grid fields, 
%  which are rectilinear everywhere but do not go north of 80N.
% 
% b) Rotate the velocities back to x-wards y-wards.  
%The array pang in regional.grid can be used to do this.  
%The eastwards,northwards to x-wards,y-wards code is in ALL/force/src/wi.f:
%
%              COSPANG  = COS(PANG(I,J))
%              SINPANG  = SIN(PANG(I,J))
%              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
%              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
%
fld = 'pang';
% --------------------
% Get grid dimensions
% --------------------
f1 = fopen(finb,'r');  % read I,J from *.b
% Read dimensions
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

n=ID;
m=JD;
l=length(FLOC);
fprintf('Reading HYCOM output, field: %s\n',fld);
fprintf('Model dim: n(x)=%i, m(y)=%i, l(z)=%i, \n',n,m,l);


fid1=fopen(fina,'r');
F=[];
for ii=1:l
  frewind(fid1);
  k0=FLOC(ii)-1;
  stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
  dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 

  dmm=reshape(dmm,ID,JD);
  F(ii,:,:)=dmm';
end;

F=squeeze(F);

fclose(fid1);

return








