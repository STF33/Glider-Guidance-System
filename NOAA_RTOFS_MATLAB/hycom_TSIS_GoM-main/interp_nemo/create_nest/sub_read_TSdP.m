% Read in T,S, dH(m) fields on HYCOM grid
% NO  averaging
% from 3-hourly fields
function [TT,SS,dH] = sub_read_TSdP(fina,finb);

rg = 9806;
hg = 2^100;
hgg= 1e20;

TT=[];
SS=[];
dH=[];

fprintf('    sub_read_TSdP: reading %s\n',fina);

if ~exist(fina,'file'), 
  fprintf('!!!    Missing Relax: %s\n',fina);
  return; 
end;


[F,n,m,l] = read_hycom(fina,finb,'temp');
F(F>hgg)=nan;
TT=F;

[F,n,m,l] = read_hycom(fina,finb,'salin');
F(F>hgg)=nan;
SS=F;

[F,n,m,l] = read_hycom(fina,finb,'thknss');
F(F>hgg)=nan;
F=F./rg; % pressure -> thicknesses, m
dH=F;

%keyboard
% Plot relax files 
f_chck=0;
if f_chck==1
%  SCT=[];
  SCT = sub_sections('ob2');
  f_lr=1;
  sttl='Relax';
  offst=10;
  sub_plot_section(TT,dH,'temp',SCT,sttl,f_lr,offst);
  keyboard
end


% convert layer thickness back to pressure
%dPav=dPav*rg;

return
