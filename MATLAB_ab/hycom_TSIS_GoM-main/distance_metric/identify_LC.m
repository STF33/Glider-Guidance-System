      function LCLCE = identify_LC(elon,alat,ssh,Bisol,varargin);
%
%  Identify eddies and LC using Bisol ssh (usually 0.17 m)
%   
%
% For a given SSH field (demeaned !), check if the SSH contour = Bisol
% goes through the Yucatan Strait and Straits of Florida,
% i.e. follows the LC (stays inside the LC)
%
% LCssh = 0 - not the LC ssh contour
%       = 1 LC ssh contour
% xx,yy - lon, lat of the contour
% 
% Input: elon, alat - 2D arrays of geogr. coordinates of the model grid
%        ssh        - demeaned SSH, 2D array
%        Bisol      - ssh contour used for LC/LCE detection 
%        optional: err_stop = 1 - the code stops returning error message if
%                  the LC was not identified
%                  nattmp - # of attempts to find the LC by changing Bisol
%                  maxBisol - max ssh contour for finding LC during
%                   LC search when initial Bisol failed
%                  minBisol - min ssh contour for LC, used if dBisol<0
%                  dBisol - change of Bisol per iteration during 
%                  LC search if initial Bisol failed
%
%  Output:   LCLCE - structured array with LC and LCEs inside
%                    1st contour is LC, other contours - LCEs arranged
%                    by size
%
%
% Dmitry Dukhovskoy, COAPS FSU, 2015-2021
%
% Please use the following reference when use this code: 
% Reference: Dukhovskoy, D.S.,, R.R. Leben, E.P. Chassignet, C. Hall, S.L. Morey, 
%              and R. Nedbor-Gross, 2015. Characterization 
%              of the Uncertainty of Loop Current Metrics using 
%              a Multidecadal Numerical Simulation and Altimeter Observations. 
%              Deep-Sea Res. I, 100, 140-158
%
% =====================================
%
f_plot=0;     %=0 - keep figure visible off
f_stop=1;     %=1 - stop with error if LC no found
nattmp=25;    % # of attempts to find the LC
maxBisol=0.25; % max contour to search LC
minBisol=0.1; % min contour
dBisol=0.01;  % >0 - contour increases when searching LC, <0 - decreases
nV = length(varargin);
if nV>0
  for k=1:nV
    vfld = varargin{k};
    if strmatch(vfld,'err_stop') % stop or not when LC is not found
      f_stop = varargin{k+1}; %
    end
    if strmatch(vfld,'nattmp'); % # of attempts to find the LC
      nattmp = varargin{k+1};
    end
    if strmatch(vfld,'maxBisol'); % max SSH contour to find LC 
      maxBisol = varargin{k+1};
    end
    if strmatch(vfld,'minBisol'); % minSSH contour for LC, if dBisol<0
      minBisol = varargin{k+1};
    end
    if strmatch(vfld,'dBisol'); % step to increase/decr contour to search LC
      dBisol = varargin{k+1};
    end
  end
end



% Define the check sections in the straits:
% Yucatan:
YS=[-86.86, 21.8; -85.2,21.8];  % Yuc. section where to search for initial point
D=distance_spheric_coord(YS(1,2),YS(1,1),alat,elon);
[j1,i1,]=find(D==min(min(D)));
IJY(1,1)=i1;
IJY(1,2)=j1;
D=distance_spheric_coord(YS(2,2),YS(2,1),alat,elon);
[j1,i1,]=find(D==min(min(D)));
IJY(2,1)=i1;
IJY(2,2)=j1;
% Get section:
jy0=IJY(1,2);
iy0=IJY(1,1);
jy1=IJY(2,2);
iy1=IJY(2,1);
if iy0>iy1; dmm=iy1; iy1=iy0; iy0=dmm; end;
LnLtY=[elon(jy0,iy0:iy1)',alat(jy0,iy0:iy1)'];  % Yucatan Chan.

% Str. Florida
FS=[-82.5, 24.2; -82.5,23.];  % Straits of Florida, check section
D=distance_spheric_coord(FS(1,2),FS(1,1),alat,elon);
[j1,i1,]=find(D==min(min(D)));
IJF(1,1)=i1;
IJF(1,2)=j1;
D=distance_spheric_coord(FS(2,2),FS(2,1),alat,elon);
[j1,i1,]=find(D==min(min(D)));
IJF(2,1)=i1;
IJF(2,2)=j1;
% Get section:
jf0=IJF(1,2);
if0=IJF(1,1);
jf1=IJF(2,2);
if1=IJF(2,1);
if jf0>jf1; dmm=jf1; jf1=jf0; jf0=dmm; end;
LnLtF=[elon(jf0:jf1,if0),alat(jf0:jf1,if0)]; % straits of florida
%keyboard

LCssh=0;

if f_plot==0
  ff=figure('Visible','off');
end

%keyboard
% Some SSH fields have coarse res.
% this may result in broken 0.17m
% near Campeche Bank - problem for LC identif.
% may need to increase Bisol - happens very rare
LCfound=logical(0);
nBisol=0;
while ~LCfound
  kk=1;
  ccT=0;
  LCLCE = struct;
  nBisol=nBisol+1;  
  if nBisol>nattmp; break; end;
  if Bisol>maxBisol, break; end;
  if Bisol<minBisol, break; end;
  if nBisol>1,
    Bisol=Bisol+dBisol;
    fprintf(' ======  identify_LC: cannot locate LC \n');
    fprintf(' ======  identify_LC: change LC isol %4.2fm \n',Bisol);
  end

  
  [cc1,cc2]=contour(elon,alat,ssh,[Bisol Bisol],'k');
  % Identify the LC contour
  np=size(cc1,2);
  d0=20000;  % max distance to the check point
%  n0=200;    % min # of points in the LC contour
%  n0=100;    % min # of points in the LC contour - resolution dependent!
% Roughly estimate min nmb of points in the LC contour in
% retracted phase - shortest distance / resol from Yuc to Fl. Str.:
  distYF = distance_spheric_coord(YS(1,2),YS(1,1),23.0,-81.5);
  dX = mean(diff(elon(end,:)));
  dY = mean(diff(alat(:,end)));
  dltX = distance_spheric_coord(YS(1,2),YS(1,1),YS(1,2),YS(1,1)+dX);
  dltY = distance_spheric_coord(YS(1,2),YS(1,1),YS(1,2)+dY,YS(1,1));
  dYF = mean([dltX,dltY]);
  n0 = round(distYF/dYF);

% if corase res, # of points is small
  if n0>length(cc1)
    n0=round(0.2*length(cc1));
  end
%keyboard  
% 

  while kk<np
    nrd=cc1(2,kk);
    iiS=kk+1;
    iiE=kk+nrd;
    xx=cc1(1,iiS:iiE)';
    yy=cc1(2,iiS:iiE)';

    chck=0;

    if nrd>n0
      clear D
      for mm=1:length(LnLtY)
        D(mm)=min(distance_spheric_coord(yy,xx,LnLtY(mm,2),LnLtY(mm,1)));
      end
      if min(D)<d0; % Check Fl. Straits:
% Check Str. FLorida
        clear D
        for mm=1:length(LnLtF)
          D(mm)=min(distance_spheric_coord(yy,xx,LnLtF(mm,2),LnLtF(mm,1)));
        end
        if min(D)<d0; 
          chck=chck+1; 
          LCLCE(1).Label='LC';
          LCLCE(1).xx=xx;
          LCLCE(1).yy=yy;
       	  LCfound=logical(1);
        end
      end;
    end

    if chck~=1 & nrd>20
%  Check if the contour is closed:
      dD=distance_spheric_coord(yy(1),xx(1),yy(2),xx(2));
      Dend=distance_spheric_coord(yy(1),xx(1),yy(end),xx(end));
      if Dend<0.01*dD, 
        ccT=ccT+1;
        LCLCE(ccT+1).Label='LCE';
        LCLCE(ccT+1).xx=xx;
        LCLCE(ccT+1).yy=yy;
      end
    end
%keyboard

    kk=iiE+1;
  end; % while kk<np  
end;  % while ~LCfound
%keyboard
LCLCE(1).Bisol = Bisol;  
if ~LCfound
  if f_stop==1
    error('*** ERR: identify_LC:  LC was not identified, check %5.2f m contour\n',Bisol);
  else
    fprintf('LC was not identified, contour Yuc-FlStr not found\n');
    LCLCE(1).xx=[];
    LCLCE(1).yy=[];
    close(ff);
    return
  end
end

%keyboard

% Small eddies can be inside bigger eddies
% check this:
inp=0;
ixp=[];
nC=length(LCLCE);
if nC>1
% First, rearrange eddies, putting smaller to the end
  for il=1:nC
    LL(il)=length(LCLCE(il).xx);
  end;
  LL(1)=10*max(LL);  % keep LC first
  [LL,IM]=sort(LL,'descend');
  dmm=LCLCE;
  LCLCE=dmm(IM);
%keyboard
  for ke=2:nC
    x0=mean(LCLCE(ke).xx);
    y0=mean(LCLCE(ke).yy);
    for pe=1:ke-1
      IN=inpolygon(x0,y0,LCLCE(pe).xx,LCLCE(pe).yy);
      if IN==1;
       	inp=inp+1;
	       ixp(inp)=ke;
      end
    end
  end
 
% Leave only not-enclosed eddies:  
  if inp>0
    dmm=LCLCE;
    clear LCLCE
    icc=0;
    for ke=1:nC
      I=find(ixp==ke);
      if isempty(I)
	       icc=icc+1;
	       LCLCE(icc)=dmm(ke);
      end
    end
  end
  
end

if chck==1;  LCssh=1; end;
if chck~=1, xx=[]; yy=[]; end;

if f_plot==0
 close(ff);
end;



return







