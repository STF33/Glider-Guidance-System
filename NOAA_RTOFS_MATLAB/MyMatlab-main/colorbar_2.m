	function [gca0,gcf0,az,hplt,hx] = ...
						colorbar_2 (gca0,gcf0,poscb,wdth,hght,mint,fsz,...
						cct,cct2,jjs,
						
						
						not finished
%
% 
% Program draws colorbar after using CONTOURF_2 to plot the data
% NOTE: The program doesn't work if the contouring has been performed
% 			with built-in matlab programs as in the 
%				contourf_2 I change the colormap
% 			such that the out-of-range values are patched with some color
% 			rather than leaving them white as Matlab does, so 
% 			my length(colormap)=length(colormap of matlab)+1
% 			and the specified intervals has one additional interval at the 
% 			lower values which is <= min(min(data matrix))   
% 
% Colorbar is placed either to the right of the plot or beneath it, specified
% by user (POSCB = 'Right','Under') (default - right position).
% Note that WDTH and HGHT are for vertical colorabr, for a 
% horizontal one they correspond as follows: WDTH = height of horiz. bar
% and HGHT = length of the bar
%
% INPUT parameters:
% 		WDTH - specified colorbar widht (default = 0.05 Normal Units)
% 		HGHT - specified colorbar height (default = 0.8) as fraction
%        			of the hight of the main figure
% 		gca0 - gca of the main figure (with contours)
%			gcf0 - gcf of the main figure
%			poscb - colorbar position ('RIGHT','UNDER')
%     mint - step in labeling on colorbar (1 - every tick is labeled, etc) 
%     fsz - font size for colorbar labels
%	    cct - major intervals (contoured with different colors on plot) 
%			jjs - color codes of major intervals 
%			cct2 - array of all contoured subintervals
%
%	Output:
%			gca1 - modified gca for the main figure (axes position)
%			az   - gca for the colorbar axes
%	    cc1, hh1  - handles for the contourf 
%     hplot - axes of the main figure
%	    hz   - axes of the colorbar 
% 		
% ===============================================================
% 

	if (isempty(wdth));  wdth=.05; 			end;
	if (isempty(hght));  hght=.8; 			end;
	if (isempty(poscb)); poscb='right'; end;
	if (isempty(mint));  mint=1; 				end;
	if (isempty(fsz));   fsz=10; 				end;
	

% =======================================
% Get information 
% ======================================
	nn=length(cct);  % # of major contoured intervals
	nsint=abs(jjs(2)-jjs(1));
	mm1=length(cct2)-2; % # of intervals on colorbar, 2 values at the
											% end are out-of-range (> and < than specified range)
 	ddy=1/mm1; % delta Y

  hplt=get(gcf0,'currentaxes');% handles for the plot axes 
	pos=get(gca0,'position');
	
% -----------------------------------------------------------
	if (poscb=='right')|(poscb=='RIGHT'); 
	
    set(gca0,'position',[pos(1),pos(2),pos(3)-.4*wdth,pos(4)]);
    wcolb=wdth;
    hcolb=hght*pos(4);
    xcolb=pos(1)+pos(3);
    ycolb=pos(2)+(pos(4)-hcolb)/2;
  	hx = axes('position',[xcolb ycolb wcolb hcolb]);
 	  set(gca,'xlim',[0 1],'ylim',[0 1])


 		for jw=1:mm1
   		xcrd(jw,1:5)=[0.,0.,0.4,0.4,0.];
 		end

 		ycrd(1,1:5)=[0.,ddy,ddy,0.,0.];
 		for jw=2:mm1
   		ycrd(jw,1)=ycrd(jw-1,2);
   		ycrd(jw,2)=ycrd(jw,1)+ddy;
   		ycrd(jw,3)=ycrd(jw,2);
   		ycrd(jw,4)=ycrd(jw,1);
   		ycrd(jw,5)=ycrd(jw,1);
 		end

% =============================================
% filling boxes with corresponding patches
% =============================================
   	hold on
   	for jp=1:mm1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp+1,:));
     set(ph(jp),'edgecolor','none');
   	end

% draw boxes

%keyboard	

 		for ll=1:nn-1
 	 		jp1=(ll-1)*nsint+1;
	 		jp2=ll*nsint;
   		plot([xcrd(jp,1) xcrd(jp,2)],[ycrd(jp1,1) ycrd(jp2,2)],'-k')
   		plot([xcrd(jp,2) xcrd(jp,3)],[ycrd(jp2,2) ycrd(jp2,3)],'-k')
   		plot([xcrd(jp,3) xcrd(jp,4)],[ycrd(jp2,3) ycrd(jp1,4)],'-k')
   		plot([xcrd(jp,4) xcrd(jp,5)],[ycrd(jp1,4) ycrd(jp1,5)],'-k')
 		end

% =============================
% Draw "Arrow heads" at the end 
% of the colorbar
% ============================
		xhead(1,1)=xcrd(1,1);
		yhead(1,1)=ycrd(1,1);
		xhead(1,2)=xcrd(1,4);
		yhead(1,2)=ycrd(1,1);
		xhead(1,3)=0.5*xcrd(1,3);
		yhead(1,3)=ycrd(1,1)-abs(ddy*nsint);
		xhead(1,4)=xcrd(1,1);
		yhead(1,4)=ycrd(1,1);

		it=(nn-1)*nsint;
		xhead(2,1)=xcrd(it,1);
		yhead(2,1)=ycrd(it,2);
		xhead(2,2)=xcrd(it,4);
		yhead(2,2)=ycrd(it,2);
		xhead(2,3)=0.5*xcrd(it,3);
		yhead(2,3)=ycrd(it,2)+abs(ddy*nsint);
		xhead(2,4)=xcrd(it,1);
		yhead(2,4)=ycrd(it,2);

		pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
		pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
	
		set(gca,'ylim',[yhead(1,3) yhead(2,3)]);

% ===================================
% Labels on the colorbar
% ===================================

 		xtxt=0.5
 		for ll=1:mint:nn,
   		jw=(ll-1)*nsint+1;
	 		if (ll==nn); jw=(ll-1)*nsint; end;
   		ytxt=ycrd(jw,1);
   		text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz);
 		end;

	end; % IF PSCB='RIGHT' loop

% --------------------------------------------------------------------
%
% =========================
% colorbar under the figure
% =========================
	if (poscb=='under')|(poscb=='UNDER'); 
	
    set(gca0,'position',[pos(1),pos(2)+.4*wdth,pos(3),pos(4)-.4*wdth]);
    wcolb=wdth;  % this is actually "height" as the colorbar under fig.
    hcolb=hght*pos(3); % this is catually "length"
    xcolb=pos(1)+(1-hght)/2;
    ycolb=pos(2);
  	hx = axes('position',[xcolb ycolb hcolb wcolb]);
 	  set(gca,'xlim',[0 1],'ylim',[0 1])


 		for jw=1:mm1
   		ycrd(jw,1:5)=[0.,0.4,0.4,0.,0.];
 		end

 		xcrd(1,1:5)=[0.,0.,ddy,ddy,0.];
 		for jw=2:mm1
   		xcrd(jw,1)=xcrd(jw-1,3);
   		xcrd(jw,2)=xcrd(jw,1);
   		xcrd(jw,3)=xcrd(jw,1)+ddy;
   		xcrd(jw,4)=xcrd(jw,3);
   		xcrd(jw,5)=xcrd(jw,1);
 		end

% =============================================
% filling boxes with corresponding patches
% =============================================
   	hold on
   	for jp=1:mm1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp+1,:));
     set(ph(jp),'edgecolor','none');
   	end

% draw boxes

%keyboard	

 		for ll=1:nn-1
 	 		jp1=(ll-1)*nsint+1;
	 		jp2=ll*nsint;
   		plot([xcrd(jp,1) xcrd(jp,2)],[ycrd(jp1,1) ycrd(jp2,2)],'-k')
   		plot([xcrd(jp,2) xcrd(jp,3)],[ycrd(jp2,2) ycrd(jp2,3)],'-k')
   		plot([xcrd(jp,3) xcrd(jp,4)],[ycrd(jp2,3) ycrd(jp1,4)],'-k')
   		plot([xcrd(jp,4) xcrd(jp,5)],[ycrd(jp1,4) ycrd(jp1,5)],'-k')
 		end

% =============================
% Draw "Arrow heads" at the end 
% of the colorbar
% ============================
		xhead(1,1)=xcrd(1,1);
		yhead(1,1)=ycrd(1,1);
		xhead(1,2)=xcrd(1,4);
		yhead(1,2)=ycrd(1,1);
		xhead(1,3)=0.5*xcrd(1,3);
		yhead(1,3)=ycrd(1,1)-abs(ddy*nsint);
		xhead(1,4)=xcrd(1,1);
		yhead(1,4)=ycrd(1,1);

		it=(nn-1)*nsint;
		xhead(2,1)=xcrd(it,1);
		yhead(2,1)=ycrd(it,2);
		xhead(2,2)=xcrd(it,4);
		yhead(2,2)=ycrd(it,2);
		xhead(2,3)=0.5*xcrd(it,3);
		yhead(2,3)=ycrd(it,2)+abs(ddy*nsint);
		xhead(2,4)=xcrd(it,1);
		yhead(2,4)=ycrd(it,2);

		pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
		pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
	
		set(gca,'ylim',[yhead(1,3) yhead(2,3)]);

% ===================================
% Labels on the colorbar
% ===================================

 		xtxt=0.5
 		for ll=1:mint:nn,
   		jw=(ll-1)*nsint+1;
	 		if (ll==nn); jw=(ll-1)*nsint; end;
   		ytxt=ycrd(jw,1);
   		text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz);
 		end;

	end; % IF PSCB='UNDER' loop
% --------------------------------------------------------------------

 set(gca,'visible','off','box','on');
 az= gca; % handle of the colorbar axes;
% ============
% Go back to original axes
% ============
      set(gcf,'currentaxes',hplt);
