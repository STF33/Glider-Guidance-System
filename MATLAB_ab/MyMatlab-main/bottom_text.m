function bottom_text(btxt0,varargin);
% function bottom_text(btxt,varargin);
% Places a string btxt at the bottom (default)
% of the figure, e.g. reference to a code name
% Specify      position: 'Position',[xleft ybtm dx dy]);
%              fontsize: 'fntsz', 12
% add current directory: 'pwd',1
%         
x1=0.011;
y1=0.012;
dx=0.8;
dy=0.05;
fntsz=10;
ff=gcf;
fpwd=0;

nV=nargin-1;
for k=1:2:nV
  aa=varargin{k};
  aa=lower(aa);
  if strncmp(aa,'position',3);
    ab=varargin{k+1};
    x1=ab(1);
    y1=ab(2);
    dx=ab(3);
    dy=ab(4);
  elseif strncmp(aa,'fntsz',3);
    fntsz=varargin{k+1};
  elseif strncmp(aa,'fghnd',3);
    ff=varargin{k+1};  % figure handle
  elseif strncmp(aa,'pwd',3);
    fpwd = 1; % automatically define dir and add to text
  end
  
end
  
vzz=get(ff,'Visible');
H = gca;

%[usr,dmm1]  = evalc('system(''whoami'')');
%[drr,dmm1]  = evalc('system(''pwd'')');
%[unm,dmm1]  = evalc('system(''uname -a'')');

if fpwd
  [a,drr]  = system('pwd');

  I=strfind(drr,'/');
  if length(I)>1
    i0=I(end-1);
  else
    i0=1;
  end

  dnm = drr(i0:end-1);
  %keyboard
  btxt = sprintf('ddmitry@yucatan.coaps.fsu.edu: %s/%s',dnm,btxt0);
  %[a,unm]  = system('uname');

else
  btxt = btxt0;
end



axes('Parent',ff,'Position',[x1 y1 dx dy]);
txx = text(0.1,0.12,btxt,'Interpreter','none','FontSize',fntsz);
set(txx,'FontSize',fntsz);
%keyboard

set(gca,'visible','off');

axes(H);
set(ff,'Visible',vzz);

return