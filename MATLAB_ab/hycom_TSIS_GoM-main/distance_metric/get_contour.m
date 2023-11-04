      function [XX,YY] = get_contour(elon,alat,ssh,Bisol,n0);
%
% The code is similar to identify_LC.m but it does not
% need a continuous contour from Yuc Ch to Fl Str
% the user specified any contour Bisol andy field
% and 
% min # of points in the contour
% and all contours will be combined into a single array of point
% 
%
% xx,yy - lon, lat of the contour
%

if isempty(n0)
  n0=200;    % min # of points in the LC contour
end
f_plot=0;     %=0 - keep figure visible off
f_stop=1;     %=1 - stop with error if LC no found



if f_plot==0
  ff=figure('Visible','off');
end

%keyboard
% Some SSH fields have coarse res.
% this may result in broken 0.17m
% near Campeche Bank - problem for LC identif.
% may need to increase Bisol - happens very rare
XX=[];
YY=[];

[cc1,cc2]=contour(elon,alat,ssh,[Bisol Bisol],'k');
np=size(cc1,2);
% 
%keyboard
kk=1;
while kk<np
		nrd=cc1(2,kk);
		iiS=kk+1;
		iiE=kk+nrd;
		xx=cc1(1,iiS:iiE)';
		yy=cc1(2,iiS:iiE)';

		if nrd>n0
				XX=[XX;xx];
				YY=[YY;yy];
		end

  kk=iiE+1;
end; % while kk<np  
%keyboard
if isempty(XX),
  error('Could not find contour %5.1f, min/max Field=%8.4f %8.4f\n',...
         Bisol,min(min(ssh)),max(max(ssh)));
end


if f_plot==0
 close(ff);
end;



return







