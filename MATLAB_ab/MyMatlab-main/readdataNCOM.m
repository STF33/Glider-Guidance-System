 function [u,v,w,t,s,e] = readdataNCOM(datun,n,m,l,uf,vf,wf,tf,sf,ef);

% function [u,v,w,t,s,e] = readdataNCOM(datun,n,m,l,uf,vf,wf,tf,sf,ef);
%%% readdata(fname,datun,n,m,l,numtimes,alat,elon,h,z,timed,u,v,w,t,s,e)
%%%
%%% read data from 1 output record 
%%%   reads the next time record of data
%%% 
%%% variables used in this procedure (first run readflagsNCOM.m):
%%%   datun - holds the unit number for the .dat file
%          the file should be open before calling this subroutine:
%          exmp= [fname,'.dat'];
%          datun = fopen(exmp,'r','b');
%%%   fname - filename for extracted data, no extention
%%%   n,m,l - x,y,z dimensions. 
%%%   numtimes - number of time steps in extracted data set
%%%   alat, elon - model latitude and longitude arrays
%%%   h - model topography
%%%   z - depths of the middle of each grid cell in extracted domain
%%%   timed - model time in days
%%%   u,v,w - z,y,z velocity components
%%%   t,s - temperature and salinity
%%%   e - surface elevation
    
Dummy=fread(datun,1,'float');
timed=fread(datun,1,'float');
Dummy=fread(datun,1,'float');

if uf == 1
  Dummy=fread(datun,1,'float');  
  for ll=1:l
    uu=fread(datun,[n m],'float');  
    u(:,:,ll)=flipud(rot90(uu));
  end
  Dummy=fread(datun,1,'float');
else 
  u=[];
end
if vf == 1   
  Dummy=fread(datun,1,'float');
  for ll=1:l  
    vv=fread(datun,[n m],'float'); 
    v(:,:,ll)=flipud(rot90(vv));
  end
  Dummy=fread(datun,1,'float');
else
   v=[];
end
if wf == 1
  Dummy=fread(datun,1,'float');
  for ll=1:l
    ww=fread(datun,[n m],'float'); 
    w(:,:,ll)=flipud(rot90(ww));
  end
  Dummy=fread(datun,1,'float');
else
  w=[];
end

if tf == 1
  Dummy=fread(datun,1,'float');
  for ll=1:l  
    tt=fread(datun,[n m],'float');  
    t(:,:,ll)=flipud(rot90(tt));
  end
  Dummy=fread(datun,1,'float');
else
  t=[];
end
  
if sf == 1
  Dummy=fread(datun,1,'float');
  for ll=1:l  
    ss=fread(datun,[n m],'float'); 
    s(:,:,ll)=flipud(rot90(ss));
  end
  Dummy=fread(datun,1,'float');
else
  s=[];
end

if ef == 1
  Dummy=fread(datun,1,'float');
  e=fread(datun,[n m],'float');
  Dummy=fread(datun,1,'float');
else
  e=[];
end
 

