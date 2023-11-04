function VL = read_layers_blkdat(fnm);
% Read z-level layer thicknesses
% from blkdat file
% that explicetely has this information
% There are several types of blkdat.input
% some have min thicknesses specified by layer

fid=fopen(fnm,'r');
VLAYERS=[];
dp0k=[];
ds0k=[];

% Read 5-line header
for k=1:4
  pp = fgetl(fid);
end
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
fprintf('blkdat: HYCOM version %i\n',dmm(1));
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
fprintf('blkdat: experiment %i\n',dmm(1));
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
IDM=dmm(1);
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
JDM=dmm(1);
fprintf('blkdat: IDM, JDM: %i %i\n',IDM);
pp = fgetl(fid);
pp = fgetl(fid);
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
kdm=dmm(1);
fprintf('blkdat: number of layers %i\n',kdm);
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
nhybrid=dmm(1);
fprintf('blkdat: hybrid levels %i\n',nhybrid);
pp = fgetl(fid);
dmm=sscanf(pp,'%i%s');
nsigma=dmm(1);
fprintf('blkdat: sigma levels %i\n',nsigma);
VL.IDM=IDM;
VL.JDM=JDM;
VL.Nvlayers=kdm;
VL.Nsigma_layers=nsigma;

% Read z-lev spacing deep region
for ik=1:kdm
  pp = fgetl(fid);
  pat={'\d?','dp0'};
  rr=regexp(pp,pat);

  if isempty(rr{2})
    fprintf('Missing field: dp0k in blkdat ...\n');
    error('blkdat.input file is not in the right format ...');
  end
  im=rr{1};
  dz=sscanf(pp(1:im(1)+6),'%f');
  dp0k(ik,1)=dz;
end

% read shallow z-level spacing
for ik=1:nsigma
  pp = fgetl(fid);
  pat={'\d?','ds0'};
  rr=regexp(pp,pat);

  if isempty(rr{2})
    fprintf('Missing field: ds0k in blkdat ...\n');
    error('blkdat.input file is not in the right format ...');
  end
  im=rr{1};
  dz=sscanf(pp(1:im(1)+6),'%f');
  ds0k(ik,1)=dz;
end

VL.dp0k=dp0k;
VL.ds0k=ds0k;

fclose(fid);

return