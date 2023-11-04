	function [az,axc]  = pcolorbar_zebra_horiz (jjs,cct,hght,lngth,...
                                         mint,fsz,bxc,posc, ...
																				 nsint2)
	
%	function [az,axc]=pcolorbar_zebra_horiz (jjs,cct,hght,lngth,...
%                                         mint,fsz,bxc,posc,nsint2)
% The program draws a horizontal colorbar under the plot
% for given array of colors jjs
% number of intervals on the colorbar = length of jjs
% with subintervals within each color
% The plot can be created using contourf or pcolor commands
% Colormap jjs is created using function clrmp_zebra
% Difference from colorbar_rainbow - out of range limits are not plotted
% and no color is assigned to the values less and larger than limits
% of cct
% 
% Input:
%       jjs   - array of colors for colormap
%       cct   - array with "major" intervals
%       hght  - height (thickness) of the colorbar
%			 lngth  - length of the colorbar
%      mint   - interval for labeling the colorbar:
%               mint=1: every interval will be labelled
%      fsz    - fontsize of the colorbar labels
%      bxc    - box color of the colorbar
%      posc  - position of the colorbar: [xcolb ycolb hcolb wcolb], 
%              default - placed just below the
%               current gca
%      nsint2 - number of subintervals within each interval
%
% Output: 
%      az     - gca of the colorbar
%      axc    - axes of the colorbar
% ====================================

% -----------------------------------
% If parameters are not specified
% use the default
% -----------------------------------
  if (isempty(hght));  hght=0.05; end;
	if (isempty(lngth)); lngth=0.8; end;
	if (isempty(mint));  mint=1;    end;
	if (isempty(fsz));   fsz=8;     end;
	if (isempty(bxc));   bxc='k';   end;

% ------------------------------------
% Create an array cct2 with subintervals
% -------------------------------------
dltD=abs(cct(2)-cct(1))/nsint2; % interval in data
nn=length(cct)-1;               % # of major intervals
nn2=0;
	
clear cct2
for kk=1:nn,
  for jj=1:nsint2,
    nn2=nn2+1;
    cct2(nn2)=cct(1)+dltD*(nn2-1); % data array with subintervals
  end;
end;
cct2(nn2+1)=cct(nn+1);
			  
% ============================================
% Drawing a colorbar
% ============================================
	
% get gca and position of the main plot
gca1=gca;
pos=get(gca1,'position');
hplt=get(gcf,'currentaxes');% handles for the plot axes 


mm1=length(cct2);
ddy=1/mm1;
% ----------------------------------
% Set off the main plot upward 
% to place the colorbar 
% ----------------------------------
if (isempty(posc));
  set(gca1,'position',[pos(1),pos(2)+1.2*hght,pos(3),pos(4)]);
  wcolb=hght;  % this is actually "height" as the colorbar under fig.
  hcolb=lngth*pos(3); % 
  xcolb=pos(1)+(1-lngth)/2;
  ycolb=pos(2);
else,
  xcolb=posc(1);
  ycolb=posc(2);
  hcolb=posc(3);
  wcolb=posc(4);
end;
hx = axes('position',[xcolb ycolb hcolb wcolb]);
set(gca,'xlim',[0 1],'ylim',[0 1])
		
%keyboard
		
% Plot the colorbar
for jw=1:mm1
  ycrd(jw,1:5)=[0.,0.4,0.4,0.,0.];
end
xcrd(1,1:5)=[0.,0.,ddy,ddy,0.];
for jw=2:mm1;
  xcrd(jw,1)=xcrd(jw-1,3);
  xcrd(jw,2)=xcrd(jw,1);
  xcrd(jw,3)=xcrd(jw,1)+ddy;
  xcrd(jw,4)=xcrd(jw,3);
  xcrd(jw,5)=xcrd(jw,1);
end;
		
% =============================================
% filling boxes with corresponding patches
% =============================================
 hold on
 %keyboard
 for jp=1:mm1-1
   ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp,:));
   set(ph(jp),'edgecolor','none');
 end
	 	
% =============
% draw boxes
% ==============
%keyboard	
   jp=1;
 		
  for ll=1:nn
    jp1=(ll-1)*nsint2+1;
    jp2=ll*nsint2;
    p1=plot([xcrd(jp1,1) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,2)],'-k');
    p2=plot([xcrd(jp1,1) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,2)],'-k');
   	p3=plot([xcrd(jp2,3) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,1)],'-k');
   	p4=plot([xcrd(jp2,3) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,1)],'-k');
	set([p1,p2,p3,p4],'Color',bxc);
  end


% ===================================
% Labels on the colorbar
% ===================================
 	ytxt=-0.5;  % x position of the text
 	for ll=1:mint:nn+1,
   	jw=(ll-1)*nsint2+1;
	 	if (ll==nn+1); jw=(ll-1)*nsint2; end;
   	xtxt=xcrd(jw,1);
		if (ll==nn+1); xtxt=xcrd(jw,3); end;
   	text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz,...
		     'HorizontalAlignment','Center');
 	end;
		
 	set(gca,'visible','off','box','on');
	
 	az  = gca; % handle of the colorbar axes;
	axc = get(gcf,'currentaxes');  % axes of the colorbar

% ============
% Go back to original axes
% ============
      set(gcf,'currentaxes',hplt);





