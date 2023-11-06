  function fig_text(txt,pos,varargin);
% fig_text(txt,pos,varargin);
% Put a string txt in a figure
% at specified location POS [xlft ylft xlength ylength]
% Options: fontsize, default=9
if nargin>0
  FS=varargin{1};
else
  FS=9;
end
hh=gca;

axes('position',pos);
ttbs=strrep(txt,'_','\_');
text(0.05,0.1,ttbs,'Fontsize',FS);
set(gca,'visible','off');
axes(hh);