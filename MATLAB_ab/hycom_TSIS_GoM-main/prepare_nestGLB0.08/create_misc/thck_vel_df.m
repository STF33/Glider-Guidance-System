% Prepare 
% Thickness and velocity diffusivity
% 2D fields
% for IAS 0.03 from GLBb0.08 GOFS3.1
%
% specified in blkdat.input:
% -0.00257 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissip.
%  -0.02    'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissip.
%   0.0    'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffus.
%   0.01   'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffus.
% if >/=0 - constant, no input file needed
%
% No land mask, so no need to change
%forfun.f:c --- veldf2 is diffusion velocity for laplacian  background diffusion
%forfun.f:c --- veldf4 is diffusion velocity for biharmonic background diffusion
%c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
%c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion 
%c ---             (negative to input spacially varying diffusion velocity)   
% diffusivity = thkdf2*dx or thkdf4*dx^3 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'IAS0.03';
E = '030';
ntopo2=09;
TV = sprintf('%2.2iDD',ntopo2);

pthtopo = '/Net/kronos/ddmitry/hycom/TSIS/topo_IAS0.03/';
pthglb = '/nexsan/people/ddmitry/hycom/GLBb0.08_GOFS3.1_expt73.7/';
pthout = '/Net/kronos/ddmitry/hycom/TSIS/relax_41lrs/';


% 
% Read global fields
IG = 4500;
JG = 3298;
fina = sprintf('%sveldf2.a',pthglb);
finb = sprintf('%sveldf2.b',pthglb);

fida = fopen(fina,'r','ieee-be');
IJG = IG*JG;
Vdf2 = fread(fida,IJG,'float32','ieee-be');
Vdf2 = reshape(Vdf2,IG,JG)';
fclose(fida);

% Read glb topo:
fina=sprintf('%sdepth_GLBb0.08_09m11.a',pthglb);
finb=sprintf('%sdepth_GLBb0.08_09m11.b',pthglb);
fida = fopen(fina,'r','ieee-be');
HHg = fread(fida,IJG,'float32','ieee-be');
HHg = reshape(HHg,IG,JG)';
fclose(fida);

I = find(HHg>1e10);
HHg=-1*HHg;
HHg(I)=100;

% GoM:
i1 = 2300;
i2 = 2900;
j1 = 1580;
j2 = 1920;
VD2 = Vdf2(j1:j2,i1:i2);
Hgom = HHg(j1:j2,i1:i2);

fprintf('Mean IAS domain veldf2=%8.6f\n',mean(mean(VD2)));
%
% Keep veldf* constant - otherwise change blkdat and create input 2D files

% veldf4:
fina = sprintf('%sveldf4.a',pthglb);
finb = sprintf('%sveldf4.b',pthglb);
fida = fopen(fina,'r','ieee-be');
IJG = IG*JG;
Vdf4 = fread(fida,IJG,'float32','ieee-be');
Vdf4 = reshape(Vdf4,IG,JG)';
Vdf4(Vdf4>1e20)=nan;
fclose(fida);
VD4 = Vdf4(j1:j2,i1:i2);


figure(1);clf;
hold on;
pcolor(VD4); shading flat
contour(Hgom,[0 0],'k');
caxis([2e-3 2.6e-3]);
colorbar

fprintf('Mean IAS domain veldf4=%8.6f\n',mean(mean(VD4)));
%
% Keep veldf* constant - otherwise change blkdat and create input 2D files





