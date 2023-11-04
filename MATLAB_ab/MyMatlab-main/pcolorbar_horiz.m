	function [az,axc]  = pcolorbar_horiz (jjs,cct,hght,lngth,...
                                         mint,fsz,bxc,posc, ...
																				 nsint2)
	
%	function [az,axc]  = pcolorbar_horiz (jjs,cct,hght,lngth,...
%                                         mint,fsz,bxc,posc, ...
%																				 nsint2)
% The program draws a horizontal colorbar under the plot
% for given array of colors jjs
% number of intervals on the colorbar = length of jjs
% with subintervals within each color
% The plot can be created using contourf or pcolor commands
% Colormap jjs is created using function clrmp_rainbow
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
	nn=length(cct);
	nn2=0;
	
  for kk=1:nn-1,
	   for jj=1:nsint2,
		   nn2=nn2+1;
		   cct2(nn2+1)=cct(1)+dltD*(nn2-1); % data array with subintervals
		 end;
	end;
			  
%	cct2(1)=cct(2)-dltD; 
	cct2(1)=cct2(1)-dltD; 
	cct2(nn2+2)=cct2(nn2+1)+dltD;


% ============================================
% Drawing a colorbar
% ============================================
	
% get gca and position of the main plot
  gca1=gca;
	pos=get(gca1,'position');
  hplt=get(gcf,'currentaxes');% handles for the plot axes 


  mm1=length(cct2)-2;
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
   for jp=1:mm1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp+1,:));
     set(ph(jp),'edgecolor','none');
   end
	 	
% =============
% draw boxes
% ==============
%keyboard	
   jp=1;
 		
   for ll=1:nn-1
 	 		jp1=(ll-1)*nsint2+1;
	 		jp2=ll*nsint2;
   		p1=plot([xcrd(jp1,1) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,2)],'-k');
   		p2=plot([xcrd(jp1,1) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,2)],'-k');
   		p3=plot([xcrd(jp2,3) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,1)],'-k');
   		p4=plot([xcrd(jp2,3) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,1)],'-k');
			set([p1,p2,p3,p4],'Color',bxc);
 		end
% =============================
% Draw "Arrow heads" at the end 
% of the colorbar
% ============================
  xhead(1,1)=xcrd(1,1);
	yhead(1,1)=ycrd(1,2);
	xhead(1,2)=xcrd(1,1);
	yhead(1,2)=ycrd(1,1);
	xhead(1,3)=xcrd(1,1)-abs(ddy*nsint2);
	yhead(1,3)=ycrd(1,1)+0.5*ycrd(1,2);
	xhead(1,4)=xcrd(1,1);
	yhead(1,4)=ycrd(1,2);

	it=(nn-1)*nsint2;
	xhead(2,1)=xcrd(it,3);
	yhead(2,1)=ycrd(it,2);
	xhead(2,2)=xcrd(it,3);
	yhead(2,2)=ycrd(it,1);
	xhead(2,3)=xcrd(it,3)+abs(ddy*nsint2);
	yhead(2,3)=ycrd(it,1)+0.5*ycrd(it,2);		
	xhead(2,4)=xcrd(it,3);
	yhead(2,4)=ycrd(it,2);
	
  set(gca,'xlim',[xhead(1,3) xhead(2,3)]);
	pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
	pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
  set([pp1, pp2],'edgecolor',bxc)	



% ===================================
% Labels on the colorbar
% ===================================
 	ytxt=-0.5;  % x position of the text
 	for ll=1:mint:nn,
   	jw=(ll-1)*nsint2+1;
	 	if (ll==nn); jw=(ll-1)*nsint2; end;
   	xtxt=xcrd(jw,1);
		if (ll==nn); xtxt=xcrd(jw,3); end;
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





