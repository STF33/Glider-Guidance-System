	function jjs = clrmp_zebra (cct,cmp,nsint2)   
% 
% =================================================================
% jjs = clrmp_zebra (cct,cmp,nsint2)
% Program creates its own colormap for specified major intervals (CCT)
% with NSINT2 subintervals within each interval
% Total # of colors (length of jjs) = (cct-1)*nsint2+1
% For the last interval of cct (maximum value) no subintervals 
% are created - just one color
%
% Note: colors are chosen from specified array of colors in CCT
%    if number of major intervals is greater than specified
%    Colors, the colors are reused in cyclic order
%
% Each major interval is assigned its own major color in default order 
% unless it is specified by user (CMP).
% CCT - specified major intervals
% CMP - array of characters specifying color orders ('G','B',...);
%       if number of colors is greater than the number of intervals
%       the last colors are ignored (length(CMP)>length(CCT))
% Default colors:
%   M  - Magenta
%   B  - Blue
%   C  - Cyen
%   G  - Green
%   Y	 - Yellow
%   O - orange
%   R  - Red
% If "white" is specified in CMP:
%   W - white out the whole interval, no shadings
%
% Inputs: 
%   CCT - specified (major) intervals for contouring
%   CMP - color order specifed by user (e.g., CMP=['R','W','C',...];
%   nsint2 - number of subintervals in each color
%            NOTE: nsint2 must be even number !!!
%	Output:
% 	jjs - colormap array: (cct-1)*nsint2+1 x 3
% ===============================================================

nn=length(cct)-1;          % # of intervals
	
% ----------------------
% Assign default values	
% ----------------------
if (isempty(nsint2)); nsint2=6; end;
if (nsint2>1),
  if (floor(nsint2/2)*2~=nsint2),
    error('Number of subintervals (nsint2) must be even');
  end
end;
%	  nsint2=floor(nsint2/2)*2;  % make sure that nsint2 is even
%	end;
	
	

% -------------------------------
% Create a colormap array CLRS
% if # of intervals > 7, colors
% will repeat 
% ------------------------------
if (isempty(cmp));
  cmp=[''];
  col=['M';'B';'C';'G';'Y';'O';'R'];
  kcheck=0;  % counter, if # of intervals > 7
  for jj=1:nn;  % 
    kcheck=kcheck+1;
    if (kcheck>7); kcheck=1; end;
    cmp(jj,1)=col(kcheck);
  end;
end;
% 
% ---------------------------------------------
% Check that each interval has assigned color
% If # of colors less than # of intervals - use in cyclic order
% ---------------------------------------------
%keyboard
  lcmp=length(cmp);
  if (lcmp<nn)
    nclr=0;
    for j=lcmp+1:nn
      nclr=nclr+1;
      if (nclr>lcmp), nclr=1; end;
      cmp(j)=cmp(nclr);
    end;
  end;
  lcmp=length(cmp);
% If # of colors is overspecified:
  if (lcmp>nn)
    cmp=cmp(1:nn);
  end;
  lcmp=length(cmp);

  if (lcmp~=nn), 
    error('Check CMP: # of colors must be equal # of intervals');
  end;
%	clrs=zeros(nn,3);
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
	
	i1=find(cmp=='W');
        iW=i1;
	if (~isempty(iW)); 
    for ll=1:length(i1),
	    clrs(iW(ll),:)=[1,1,1]; 
    end;
	end;

 if (length(clrs)~=nn),
   error('# of colors CLRS must be = # of intervals'); 
 end;

% ---------------------------------------------------------
% Creat new arrays of data intervals (cct2) and colormap (jjs)
% with specified subintervals
% ---------------------------------------------------------

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
for kk=1:nn; % loop for all major intervals  
%%%%%		
  mm=0;
  for jj=1:nsint;
	nn2=nn2+1;
	mm=mm+1;
	cct2(nn2+1)=cct(1)+dltD*(nn2-1); % data array with subintervals
	j1=clrs(kk,1)+0.8-dltC*(mm-1); % color changes from light to dark
	j2=clrs(kk,2)+0.8-dltC*(mm-1);
	j3=clrs(kk,3)+0.8-dltC*(mm-1);

	if (~isempty(iW)) & (kk==iW); 
	  j1=1; j2=1; j3=1; % White region does not change
	end; 
		 
	if (j1>1); j1=1; end;
	if (j2>1); j2=1; end;
	if (j3>1); j3=1; end;
	jjs(nn2,:)=[j1,j2,j3];
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
		 
    if (~isempty(iW)) & (kk==iW); 
      j1=1; j2=1; j3=1; % White region does not change
    end; 
	 
    if (j1<0); j1=0; end;
    if (j2<0); j2=0; end;
    if (j3<0); j3=0; end;
    jjs(nn2,:)=[j1,j2,j3];
  end;

end;  % for kk 

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
		      jjs((ior(ll)-1)*nsint2+mm,:)=[j1,j2,j3]; % +1 is offset for outcolr
		    end;
        for mm=1:nsint;
		      j1=1-dltC*(mm-1);
			    j3=0;
			    j2=(j1+j3)/2.;	   
		      jjs((ior(ll)-1)*nsint2+mm+nsint,:)=[j1,j2,j3];
		    end;


      end; % for ll=... several orange intervals
    end; % if ~isempty

%keyboard
% =========================
  nsint=nsint2;  % for colorbar 

% ----------------------------------------------
% flip up-down colors within each major interval,
% so that, the shades go : dark - to - light
% ----------------------------------------------
%keyboard
for kkk=1:1:nn,
  i1=(kkk-1)*nsint2+1;
  i2=kkk*nsint2;
  AA=jjs(i1:i2,:);
  AA=flipud(AA);
  jjs(i1:i2,:)=AA;
end;
		
% H=fill([0 0 1 1],[0 1 1 0],'r');
% set(H,'FaceColor',[.9 .9 .9]);

%keyboard

