	function [hh] = match_colors(jjc,hh,cct);
	
% =========================================================
% [hh] = match_colors(jjc,hh,cct)
% This function corrects some bugs in the Matlab
% contouring in the contourf command: colors of the
% colormap do not match the specified intervals
% 
% Input parameters:
% jjc - colormap handle, i.e. array of colors for N intervals
%				where N=length(jjc)
% hh - handle of contourf command with all the information
%     about contoured data and patches ([cc,hh]=contourf(...))
% cct - contoured intervals
% ========================================================

vv=version;    % Matlab version, new version has different
               % structure of contourf handle

for in = 1:length(hh),
  datint = get(hh(in),'CData');  % find which data corresponds to 
  where = max(find(cct-datint <= 0)); % compare this data with 
				% given interval cct
				% when =0 -> "where" is the 
				% index of corresp. data from interval cct
				% Now, we need to take the color of this data
				% and  assign it to plotted data
  if ~isempty(where),
%        if (jjc(where,:)==[1.0000,1.0000,0]), % dye one yellow patch
%            jjc(where,:)=[.9,.8,0];               
%        end
    set(hh(in),'FaceColor',jjc(where,:));
  end
end
