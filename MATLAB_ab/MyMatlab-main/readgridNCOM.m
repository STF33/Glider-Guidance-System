  function [h,alat,elon,z] = readgridNCOM(fname,n,m,l)
%  function [h,alat,elon,z] = readgridNCOM(fname,n,m,l)
%     readgridNCOM(fname,n,m,l)
% 
%     read grid info
%%%   reads dimensions and flags
%%%   reads grid data alat,elon,h,z
%%%   initializes array space
% 
% Input data: n,m,l - grid dimensions obtained from readflagsNCOM.m
% ---------------------------------------

  NCOMgrid = [fname,'.grd'];
%%%  OPEN uvwtse3d.grd AND READ GRID INFO  %%%

  grid=fopen(NCOMgrid,'r','b');
  Dummy=fread(grid,1,'float');
  h=fread(grid,[n m],'float');
  Dummy=fread(grid,1,'float');
  h=h';                        % orient matrix such that Y is rows, X - columns

  Dummy=fread(grid,1,'float');
  alat=fread(grid,[n m],'float');
  Dummy=fread(grid,1,'float');
  alat=alat';

  Dummy=fread(grid,1,'float');
  elon=fread(grid,[n m],'float');
  Dummy=fread(grid,1,'float');
  elon=elon';

  Dummy=fread(grid,1,'float');
  for ll=1:l;
    zz=fread(grid,[n m],'float');
    z(:,:,ll)=flipud(rot90(zz));
  end;
  Dummy=fread(grid,1,'float');

  fclose(grid);

  clear Dummy ll grid NCOMgrid ans;

    %%%  OPEN .dat FILE FOR LATER  %%%
%    datun=fopen(NCOMdata,'r','b');
