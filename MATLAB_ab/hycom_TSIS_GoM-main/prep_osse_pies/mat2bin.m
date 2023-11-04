% convert mat format files to bin
% ssh fields
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

YR=2011;

% Synthetis OSSE observations
% frmo NEMO 1km free run 2011
pthosse = '/Net/gleam/dmitry/hycom/TSIS/IASx0.03/OSSE/';
pthout  = pthosse;
pthtopo = '/home/ddmitry/codes/HYCOM_TSIS/';


fmat=sprintf('%salongtrack_2011.mat',pthosse);
load(fmat);
A=alongtrack;
ll=length(A.lon);
nbit=ll*8;


TM=A.time;
dv1=datevec(TM(1));
dj1=datenum(dv1(1),1,1);
Time=TM-dj1;


fout=sprintf('%salongtrack_2011.dat',pthout);
fid = fopen(fout,'w','ieee-be');

fwrite(fid,nbit,'int');
fwrite(fid,ll,'int');  % length
fwrite(fid,nbit,'int');

fwrite(fid,nbit,'int');
fwrite(fid,A.lon,'float64');
fwrite(fid,nbit,'int');

fwrite(fid,nbit,'int');
fwrite(fid,A.lat,'float64');
fwrite(fid,nbit,'int');

fwrite(fid,nbit,'int');
fwrite(fid,Time,'float64');
fwrite(fid,nbit,'int');

fwrite(fid,nbit,'int');
fwrite(fid,A.ssh,'float64');
fwrite(fid,nbit,'int');

fclose(fid);

fprintf(' Data are in %s\n',fout);
