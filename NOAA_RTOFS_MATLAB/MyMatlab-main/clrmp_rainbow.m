	function jjs = clrmp_rainbow (cct,cmp,nsint2)   
% 
% =================================================================
% jjs = clrmp_rainbow (cct,cmp,nsint2)
% Program creates its own colormap for specified major intervals (CCT)
% with NSINT2 subintervals within each interval
%
% Recommended number of major intervals (CCT) is 7 = 7 rainbow colors
% If number of intervals is greater than 7, the colors will repeat
%
% Each major interval is assigned its own major color in default order 
% unless it is specified by user (CMP).
% CCT - specified major intervals
% CMP - array of characters specifying color orders ('G','B',...);
% Maximum Number of NSINT = max # of major colors which is 8:
%   B  - Blue
%   M  - Magenta
%   C  - Cyen
%   O - orange
%   Y	 - Yellow
%   G  - Green
%   R  - Red
%   W  - white
%   K  - Black
%
% Inputs: 
%   CCT - specified (major) intervals for contouring
%   CMP - color order specifed by user (e.g., CMP=['R','W','C',...];
%   nsint2 - number of subintervals in each color
%            NOTE: nsint2 must be even number !!!
%	Output:
% 	jjs - colormap array: NSINT2+2 x 3, there are 2 more colors
%         added to jjs for out-of the specified range values
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
%	cct2(1)=min(min(min(ZZ)), cct(2)-dltD); 
  cct2(1)=cct(2)-dltD; 
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

