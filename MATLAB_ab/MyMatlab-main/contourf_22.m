	function [gca1,az,cc1,hh1,hplt,hx] = ...
									contourf_22 (XX,YY,ZZ,cct,cmp,nsint2,COLB,poscb,...
															wdth,hght,mint,fsz)   
% 
% =================================================================
% [gca1,az,cc1,hh1,hplt,hx] = ...
%  	contourf_22 (XX,YY,ZZ,cct,cmp,nsint2,COLB,poscb,wdth,hght,mint,fsz) 
%
%  Program creates its own colormap and draws 
% the colorbar for specified major intervals (CCT)
% If the data interval contains 0, i.e. there are negative and
% positive intervals, it is better to have 0 as one of the interval limits
% E.g., CCT=(-4:2:4), will give: -4, -2, 0, 2, 4;
% 
% Recommended number of major intervals (CCT) is 7 = 7 rainbow colors
% If number of intervals is greater than 7, the colors will repeat
%
% Each major interval is assigned its own color in default order 
% unless it is specified by user (CMP).
% CCT - specified major intervals
% CMP - array of characters specifying color orders ('G','B',...);
%   B  - Blue
%   M  - Magenta
%   C  - Cyen
%   O - orange
%   Y	 - Yellow
%   G  - Green
%   R  - Red
%   W  - white
%   K  - Black
% Recommended # of CCT not to exceed 7, white color is used only 
% to wite out intervals which are out of interest (for example,
% interval around zero, -2 to 2 can be white).
% Each color interval (except for WHITE) is then subdivided
% into NSINT subintervals, 1<=NSINT<=20. If NSINT=1, then
% no subintervals are drawn. 
% Patches in the figure are redrawn according to each subinterval color.
% 
% Colorbar is placed either to the right of the plot or beneath it, specified
% by user (POSCB = 'Right','Under') (default - right position).
% Note that WDTH and HGHT are for vertical colorabr, for a 
% horizontal one they correspond as follows: WDTH = height of horiz. bar
% and HGHT = length of the bar
%
% Inputs: 
%   WDTH - specified colorbar widht (default = 0.05 Normal Units)
%   HGHT - specified colorbar height (default = 0.95) as fraction
%          of the hight of the main figure
%   XX - Arrays of X coordinates (NxM, or Mx1)
%   YY - Arrays of Y  coordinates (NxM, or Nx1)
%   ZZ - Data matrix (NxM)
%   CCT - specified (major) intervals for contouring
%   CMP - color order specifed by user (e.g., CMP=['R','W','C',...];
%   nsint2 - number of subintervals in each color
%            NOTE: nsint2 must be even number !!!
%   gca - gca of the main figure
%   COLB - flag for colorbar, 'N' or 'n','NO','no' - no colorbar
%   poscb - colorbar position: 'Right','Under', default - right position 
%   mint - step in labeling on colorbar, 
%   fsz - font size for colorbar labels
%
%	Output:
%			gca1 - modified gca for the main figure (axes position)
%			az   - gca for the colorbar axes
%	    cc1, hh1  - handles for the contourf 
%     hplot - axes of the main figure
%	    hx   - axes of the colorbar 
% 		
% ===============================================================

	nn=length(cct);
	
	if (~isempty(cmp));
		if (nn-1~=length(cmp));
			disp('Number of colors in cmp does not match # of intervals');
			pause;
			return;
		end;
	end;

% ----------------------
% Assign default values	
% ----------------------
	if (isempty(nsint2)); nsint2=6; end;
	if (nsint2>1),
	  nsint2=floor(nsint2/2)*2;  % make sure that nsint2 is even
	end;
	if (isempty(wdth)); wdth=.05; end;
	if (isempty(hght)); hght=.95; end;
	if (isempty(poscb)); poscb='right'; end;
	if (isempty(mint)); mint=1; end;
	if (isempty(fsz)); fsz=10; end;
	if (isempty(COLB)); COLB='Yes'; end;
	
	

% -------------------------------
% Create a colormap array CLRS
% if # of intervals > 7, colors
% will repeat 
% ------------------------------
	if (isempty(cmp));
	  cmp=[''];
	  col=['M';'B';'C';'G';'Y';'O';'R'];
    kcheck=0;  % counter, if # of intervals > 7
		for jj=1:nn;  % ???? should be nn-1
      kcheck=kcheck+1;
      if (kcheck>7); kcheck=1; end;
			cmp(jj,1)=col(kcheck);
		end;
	end;

	clrs=zeros(nn,3);
	i1=find(cmp=='M');
	if (~isempty(i1)); 
    for ll=1:length(i1),
        clrs(i1(ll),:)=[1,0,1]; 
    end;
  end;

	i1=find(cmp=='B');
	if (~isempty(i1)); 
    for ll=1:length(i1),
       clrs(i1(ll),:)=[0,0,1]; 
    end;
	end;
	
	i1=find(cmp=='C');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[0,1,1]; 
    end;
	end;
	
	i1=find(cmp=='Y');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[1,1,0]; 
    end;
	end;
	
	i1=find(cmp=='G');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	     clrs(i1(ll),:)=[0,1,0]; 
    end;
	end;

	i1=find(cmp=='O');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[1,.5,0]; 
    end;
	end;
	
	i1=find(cmp=='R');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[1,0,0]; 
    end;
	end;
	
	i1=find(cmp=='K');
	if (~isempty(i1)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[0,0,0]; 
    end;
	end;
	
	iW=find(cmp=='W');
	if (~isempty(iW)); 
    for ll=1:length(i1),
	    clrs(i1(ll),:)=[1,1,1]; 
    end;
	end;
  

% ---------------------------------------------------------
% Creat new arrays of data intervals (cct2) and colormap (jjs)
% with specified subintervals
% ---------------------------------------------------------

%keyboard
% ------------------------
% the color is changed in two steps,
% first changing RGB index from 0.8 to 0:
% e.g. if CCT2(1,:)=1,0,1 then, it goes
% from 1, .8, 1 (light magenta) to 1,0,1
% and then 1,0,1 changes to .2,0,.2
% ----------------- 
%  
  
	nsint=floor(nsint2/2.); % if nsint2=1, nsint=0
	if (nsint>=1),
	  dltC=0.8/nsint; % interval in colormap
	else,
	  dltC=0;
	end;
	
	dltD=abs(cct(2)-cct(1))/nsint2; % interval in data
	
	clear jjs cct2 
	nn2=0;  % First color for the out of the specified CCT range value
% ------------------------------------
	for kk=1:nn-1; % loop for all major intervals
%%%%%		
    mm=0;
		for jj=1:nsint;
		 nn2=nn2+1;
		 mm=mm+1;
		 cct2(nn2+1)=cct(1)+dltD*(nn2-1); % data array with subintervals
		 j1=clrs(kk,1)+0.8-dltC*(mm-1); % color changes from light to dark
		 j2=clrs(kk,2)+0.8-dltC*(mm-1);
		 j3=clrs(kk,3)+0.8-dltC*(mm-1);

%  if there are negative intervals:
%		 if (~isempty(mneg))&(kk<=mneg);
%		 		j1=clrs(kk,1)+dltC*(mm-1); % for negative values, flip up colors
%		 		j2=clrs(kk,2)+dltC*(mm-1); % R,G, or B changes from 0 to 0.8
%		 		j3=clrs(kk,3)+dltC*(mm-1);
%		 end;

		 if (~isempty(iW)) & (kk==iW); 
		 		j1=1; j2=1; j3=1; % White region does not change
		 end; 
		 
		 if (j1>1); j1=1; end;
		 if (j2>1); j2=1; end;
		 if (j3>1); j3=1; end;
		 jjs(nn2+1,:)=[j1,j2,j3];
		end;
%%%%%

	   mm=0;

		for jj=nsint+1:nsint2;
		 nn2=nn2+1;
		 mm=mm+1;
		 cct2(nn2+1)=cct(1)+dltD*(nn2-1);
		 j1=clrs(kk,1)-dltC*(mm-1);
		 j2=clrs(kk,2)-dltC*(mm-1);
		 j3=clrs(kk,3)-dltC*(mm-1);
		 
%		 if (~isempty(mneg))&(kk<=mneg);
%		 		j1=clrs(kk,1)+dltC*(mm-nsint); % for negative values, flip up colors
%		 		j2=clrs(kk,2)+dltC*(mm-nsint);
%		 		j3=clrs(kk,3)+dltC*(mm-nsint);
%		 end;

		 if (~isempty(iW)) & (kk==iW); 
		 		j1=1; j2=1; j3=1; % White region does not change
		 end; 
		 
		 if (j1<0); j1=0; end;
		 if (j2<0); j2=0; end;
		 if (j3<0); j3=0; end;
		 jjs(nn2+1,:)=[j1,j2,j3];
		end;
		

	end; 

%keyboard
	
% ==============================
%explicitly define orange...kk=6

    clear ior
    ior=find(cmp=='O');
		
		if (~isempty(ior)),
      if (nsint>=1),
	      dltC=0.7/(nsint); % interval in colormap
	    else,
	      dltC=0;
	    end;
	
		  for ll=1:length(ior), % if there are several Orange intervals
		
        for mm=1:nsint;
		      j1=1;
			    j3=0.7-dltC*(mm-1);
			    j2=(j1+j3)/2.;	   
		      jjs((ior(ll)-1)*nsint2+mm+1,:)=[j1,j2,j3]; % +1 is offset for outcolr
		    end;
        for mm=1:nsint;
		      j1=1-dltC*(mm-1);
			    j3=0;
			    j2=(j1+j3)/2.;	   
		      jjs((ior(ll)-1)*nsint2+mm+nsint+1,:)=[j1,j2,j3];
		    end;


      end; % for ll=... several orange intervals
    end; % if ~isempty

% =========================
  nsint=nsint2;  % for colorbar 

% Specify out of the range values and colors
	j1=jjs(2,1)+dltC;
	j2=jjs(2,2)+dltC;
	j3=jjs(2,3)+dltC;

	if (j1>1); j1=1; end;
	if (j2>1); j2=1; end;
	if (j3>1); j3=1; end;

%	j1=jjs(2,1)-dltD;
%	j2=jjs(2,2)-dltD;
%	j3=jjs(2,3)-dltD;
		 
	if (j1<0); j1=0; end;
	if (j2<0); j2=0; end;
	if (j3<0); j3=0; end;
	
	jjs(1,:)=[j1,j2,j3];
	 
	j1=jjs(nn2+1,1)-dltC;
	j2=jjs(nn2+1,2)-dltC;
	j3=jjs(nn2+1,3)-dltC;
		 
	if (j1<0); j1=0; end;
	if (j2<0); j2=0; end;
	if (j3<0); j3=0; end;
	
	jjs(nn2+2,:)=[j1,j2,j3];
	cct2(1)=min(min(min(ZZ)), cct(2)-dltD); 
	cct2(nn2+2)=cct2(nn2+1)+dltD;
	
% ----------------------------------------------
% if there are negative values - flip up-down
% the colormap for that values
% ----------------------------------------------
%keyboard

	ineg=find(cct<0);
	mneg=length(ineg);
	
  if (mneg>0),
	  for kkk=1:mneg,
		  i1=(kkk-1)*nsint2+1+1;
			i2=kkk*nsint2+1;
			AA=jjs(i1:i2,:);
			AA=flipud(AA);
			jjs(i1:i2,:)=AA;
		end;
		
% Out of range colors:
		
	  j1=jjs(2,1)-dltC;
	  j2=jjs(2,2)-dltC;
	  j3=jjs(2,3)-dltC;

	  if (j1<0); j1=0; end;
	  if (j2<0); j2=0; end;
	  if (j3<0); j3=0; end;
		
	  jjs(1,:)=[j1,j2,j3];
	end;

%keyboard

  [cc1,hh1] = contourf(XX,YY,ZZ,cct2);
	set(hh1,'edgecolor','none');
% ===============================================
% Make the colors to match the specified 
% intervals
% ==============================================	
	hh1 = match_colors(jjs,hh1,cct2);
	gca1=gca;
  hplt=get(gcf,'currentaxes');% handles for the plot axes 
	pos=get(gca1,'position');
  
		
% ==============================================================
% ==============================================================

%keyboard

% ============================================
% Drawing a colorbar
% ============================================
	if (COLB(1)=='n')|(COLB(1)=='N'); 
%		disp('No colorbar');
		az=[];
		hx=[]; 
		return; 
	end;
	
	mm1=length(cct2)-2; % # of intervals on colorbar
  ddy=1/mm1; % delta Y
	
% ----------------------------------------
	if (poscb=='right')|(poscb=='RIGHT'); 
	
    set(gca1,'position',[pos(1),pos(2),pos(3)-1.15*wdth,pos(4)]);
    wcolb=wdth;
    hcolb=hght*pos(4);
    xcolb=pos(1)+pos(3)-.85*wdth;
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
	end;
% ----------------------------------------
	
	if (poscb=='under')|(poscb=='UNDER'); 
	
    set(gca1,'position',[pos(1),pos(2)+1.2*wdth,pos(3),pos(4)]);
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
 		for jw=2:mm1;
   		xcrd(jw,1)=xcrd(jw-1,3);
   		xcrd(jw,2)=xcrd(jw,1);
   		xcrd(jw,3)=xcrd(jw,1)+ddy;
   		xcrd(jw,4)=xcrd(jw,3);
   		xcrd(jw,5)=xcrd(jw,1);
 		end;
		
  end;
% ----------------------------------------




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
% ----------------------------------------
	if (poscb=='right')|(poscb=='RIGHT'); 
 		for ll=1:nn-1
 	 		jp1=(ll-1)*nsint+1;
	 		jp2=ll*nsint;
   		plot([xcrd(jp,1) xcrd(jp,2)],[ycrd(jp1,1) ycrd(jp2,2)],'-k')
   		plot([xcrd(jp,2) xcrd(jp,3)],[ycrd(jp2,2) ycrd(jp2,3)],'-k')
   		plot([xcrd(jp,3) xcrd(jp,4)],[ycrd(jp2,3) ycrd(jp1,4)],'-k')
   		plot([xcrd(jp,4) xcrd(jp,5)],[ycrd(jp1,4) ycrd(jp1,5)],'-k')
 		end
	end;
% ----------------------------------------

	if (poscb=='under')|(poscb=='UNDER'); 
 		for ll=1:nn-1
 	 		jp1=(ll-1)*nsint+1;
	 		jp2=ll*nsint;
   		plot([xcrd(jp1,1) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,2)],'-k')
   		plot([xcrd(jp1,1) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,2)],'-k')
   		plot([xcrd(jp2,3) xcrd(jp2,3)],[ycrd(jp,2) ycrd(jp,1)],'-k')
   		plot([xcrd(jp2,3) xcrd(jp1,1)],[ycrd(jp,1) ycrd(jp,1)],'-k')
 		end
	end;
% ----------------------------------------
	
	
% =============================
% Draw "Arrow heads" at the end 
% of the colorbar
% ============================
% ----------------------------------------
	if (poscb=='right')|(poscb=='RIGHT'); 

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
	
	  set(gca,'ylim',[yhead(1,3) yhead(2,3)]);
	
	end;
% ----------------------------------------

	if (poscb=='under')|(poscb=='UNDER'); 
		
		xhead(1,1)=xcrd(1,1);
		yhead(1,1)=ycrd(1,2);
		xhead(1,2)=xcrd(1,1);
		yhead(1,2)=ycrd(1,1);
		xhead(1,3)=xcrd(1,1)-abs(ddy*nsint);
		yhead(1,3)=ycrd(1,1)+0.5*ycrd(1,2);
		xhead(1,4)=xcrd(1,1);
		yhead(1,4)=ycrd(1,2);

		it=(nn-1)*nsint;
		xhead(2,1)=xcrd(it,3);
		yhead(2,1)=ycrd(it,2);
		xhead(2,2)=xcrd(it,3);
		yhead(2,2)=ycrd(it,1);
		xhead(2,3)=xcrd(it,3)+abs(ddy*nsint);
		yhead(2,3)=ycrd(it,1)+0.5*ycrd(it,2);
		xhead(2,4)=xcrd(it,3);
		yhead(2,4)=ycrd(it,2);
	
	  set(gca,'xlim',[xhead(1,3) xhead(2,3)]);

	end;
		
% ----------------------------------------


	pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
	pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
	
%keyboard
% ===================================
% Labels on the colorbar
% ===================================
% ----------------------------------------
	if (poscb=='right')|(poscb=='RIGHT'); 
	
 		xtxt=.6;
 		for ll=1:mint:nn,
   		jw=(ll-1)*nsint+1;
	 		if (ll==nn); jw=(ll-1)*nsint; end;
   		ytxt=ycrd(jw,1);
			if (ll==nn); ytxt=ycrd(jw,2); end;
   		text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz,...
			'HorizontalAlignment','Left');
 		end;

	end;
% ----------------------------------------

	if (poscb=='under')|(poscb=='UNDER'); 

 		ytxt=-0.5;
 		for ll=1:mint:nn,
   		jw=(ll-1)*nsint+1;
	 		if (ll==nn); jw=(ll-1)*nsint; end;
   		xtxt=xcrd(jw,1);
			if (ll==nn); xtxt=xcrd(jw,3); end;
   		text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz,...
			'HorizontalAlignment','Center');
 		end;
		
	end;
% ----------------------------------------
 	set(gca,'visible','off','box','on');
 	az= gca; % handle of the colorbar axes;
		
% ============
% Go back to original axes
% ============
      set(gcf,'currentaxes',hplt);


