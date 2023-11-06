function [n,m,l,numtimes,uf,vf,wf,tf,sf,ef]=readflagsNCOM(fname);
%function [n,m,l,numtimes,uf,vf,wf,tf,sf,ef]=readflagsNCOM(fname);
%%%  Reads parameters and flags from NCOM output ascii descriptor file
%%%  filename - data filename without extension
%%%  n,m,l - dimensions for extracted data in x,y,z
%%%  numtimes - number of extracted time steps
%%%  uf,vf,wf,tf,sf,ef - logical flags for whether u,v,w,t,s,e are saved

  NCOMtext = [fname,'.txt'];
  textid=fopen(NCOMtext,'r');
  %%%  READ DATA PROPERTIES FROM .txt descriptor file %%%

  is=sscanf(fgetl(textid),'%f');
  ie=sscanf(fgetl(textid),'%f');
  n=sscanf(fgetl(textid),'%f');
  js=sscanf(fgetl(textid),'%f');
  je=sscanf(fgetl(textid),'%f');
  m=sscanf(fgetl(textid),'%f');
  ks=sscanf(fgetl(textid),'%f');
  ke=sscanf(fgetl(textid),'%f');
  l=sscanf(fgetl(textid),'%f');
  numtimes=sscanf(fgetl(textid),'%f');
  uf=sscanf(fgetl(textid),'%f');
  vf=sscanf(fgetl(textid),'%f');
  wf=sscanf(fgetl(textid),'%f');
  tf=sscanf(fgetl(textid),'%f');
  sf=sscanf(fgetl(textid),'%f');
  ef=sscanf(fgetl(textid),'%f');

  st=fclose(textid);

  disp(['Output is read from:',NCOMtext]);
  disp(['X dim, n=',int2str(n),' Y dim, m=',int2str(m),...
       ' Z dim, l=',int2str(l)]);
  disp(['Extracted domain:']);
  disp(['X lft bndry,   is=',int2str(is),';  X Right, is=',int2str(ie)]);
  disp(['Y Lower bndry, js=',int2str(js),';  Y Upper, je=',int2str(je)]);
  disp(['Z srf bndry,   ks=',int2str(ks),';  Z btm,   ke=',int2str(ke)]);
  disp(['# of output records extracted: ',int2str(numtimes)]);

  clear st textid NCOMtext

