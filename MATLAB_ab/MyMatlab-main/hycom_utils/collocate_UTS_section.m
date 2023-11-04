function UTS=collocate_UTS_section(I,J,U,V,T,S,dH,HH,LON,LAT,nin);
% This is for NOT COLLOCATED U,V components !!!!
% Check if the components have been collocated
% in the archv files
% 
%
% Collocate U,V,T,S variables
% along section or contour (I,J)
% in HYCOM grid
% dH - layer thicknesses
% HH - bathymetry 2D field
% nin - specify norm directed inside the contour (=1) or outside (=-1)
%     - if nin>1 - then this is an index that 
%                  indicates positive direction of the norm
%                  along the section
%                  
%  nin=1 or -1 is better to use for the contours/closed polygons
%  nin=index >1 - better to use for sections where inside/outside
%                cannot be properly defined
%         be cautious using index for polygons
%         as some segments on the contours may have positive direction
%         towards this point, when pointing outside the contour
%
%  
%      =[] - skip norm - not recommended for sections/contours
%            with "stepping" segments as along the steps
%            positive U or V component should in fact give
%            a negative flux !!!
%
% Return, velocity component normal to the segment
% T,S collocated with normal component
% and norm vectors at each segment (if nin specified)
%
% The contour always goes through u(i,j) or v(i,j) points
% So that the direction you follow the contour doesn't matter
% The problem is "corner" points: in some cases (depending on the
% angle the contour changes) both u(i,j) and v(i,j) have to be
% added (to close off the contour), in others none are needed
% (interior corner grid cell)
% This is determined based on the orientation of P-vector (local
% direction) and u,v points wrt p-vector
%
%
%  HYCOM grid:
%
%          ------------   | V(i,j+1) ------------------
%          |                             |
%          |                             |
%        --- U(i,j)      *P(i,j)        ---  U(i+1,j)
%
%          |                             |
%          |                             |
%          ------------   | V(i,j) ------------------
%
%  Section follows U-V points
%  not centroids (P-points)
%   
%   
btx='collocate_UTS_section.m';

fprintf('collocate_UTS_section\n');


UTS = struct;
UTS.Norm=[];

HH(isnan(HH))=1e-10;
HH(HH>0)=1e-10;
dH(isnan(dH))=0;

ux=[1,0];
uy=[0,1];

nsc=length(I);
cc=0;
for is=1:nsc
  if mod(is,200)==0, fprintf(' %5.2f done ...\n',is/nsc*100); end

  i0=I(is);
  j0=J(is);
  if is<nsc
    ip1=I(is+1);
    jp1=J(is+1);
  else
    di=I(is)-I(is-1);
    dj=J(is)-J(is-1);
    ip1=I(is)+di;
    jp1=J(is)+dj;
  end
  if is>1
    im1=I(is-1);
    jm1=J(is-1);
  else
    di=I(is+1)-I(is);
    dj=J(is+1)-J(is);
    im1=I(is)-di;
    jm1=J(is)-dj;
  end
   
  
% Get local direction pn vector and 
% un orientation wrt to pn
  pnx=ip1-i0;
  pny=jp1-j0;
  pn=[pnx;pny;0]./(sqrt(pnx^2+pny^2));

% un = u or v component:
% en is vector pointing towards U or V component
% from the p-point wrt to p-vector direction (contour)
  if pnx==0   % vertical segment
    en=[-1;0;0];
  else
    en=[0;-1;0];
  end
  
  en_pn=(en(1)*pn(2)-en(2)*pn(1));
  
% Same for the previous step:
  if is==1
    en_pnL=en_pn;
    pnL=pn;
  else
    pnxL=i0-im1;
    pnyL=j0-jm1;
    pnL=[pnxL;pnyL;0]./(sqrt(pnxL^2+pnyL^2));
    if pnxL==0
      enL=[-1;0;0];
    else
      enL=[0;-1;0];
    end
  
    en_pnL=(enL(1)*pnL(2)-enL(2)*pnL(1));
  end
%
% Has the local orientation of the contour changed?
  pn_pnL=(pn(1)*pnL(2)-pn(2)*pnL(1));
  
%  keyboard
  
  if HH(j0,i0)>=1e20,  % turn this off 
    cc=cc+1;
    UTS.Hb(is)=HH(j0,i0);
    UTS.Unrm(:,cc)=V(:,j0,i0);
    UTS.Tnrm(:,cc)=T(:,j0,i0);
    UTS.Snrm(:,cc)=S(:,j0,i0);
    UTS.dH(:,cc)=dH(:,j0,i0);

%  ==============    
% Checking:
    ic=i0;
    jc=j0;
    if ic==i0
      ic1=ic-0.5;
      ic2=ic+0.5;
      jc1=jc;
      jc2=jc;
    else
      ic1=ic;
      ic2=ic;
      jc1=jc-0.5;
      jc2=jc+0.5;
    end
    
    UTS.sgmX(cc,:)=[ic1,ic2];
    UTS.sgmY(cc,:)=[jc1,jc2];
    UTS.Ictr(cc)=is;
%  ==============    
    continue; 
  end;
%
  
% Decide which components is needed:
% on a straight segments - only 1 components needed
%   section goes either through u or v-point
%   corner points - open corner - both components needed
%   closed corner - none
% u,v, both or none - depending
% on the direction of the contour & corner points
% when the contour turns, the uv-contour may cross
% the p-contour (going through p-points)
% this results in switching of the en-vector direction
% that point local direction from u/v point to the contour
% going through the p-points
  if en_pn==en_pnL & pn_pnL==0  % nothing changed, keep going
    if pny==0  % horizontal segment
      j1=j0-1;
      i1=i0;
      CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
    else
      j1=j0;
      i1=i0-1;
      CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
    end
    
    cc=cc+1;

    In=find(CLC.dHn<1e-3);
    CLC.Tn(In)=nan;
    CLC.Sn(In)=nan;
    CLC.Un(In)=nan;
  
    UTS.Hb(cc)=-CLC.Hn;
    UTS.Unrm(:,cc)=CLC.Un;
    UTS.Tnrm(:,cc)=CLC.Tn;
    UTS.Snrm(:,cc)=CLC.Sn;
    UTS.dH(:,cc)=CLC.dHn;
%  ==============    
% Checking:
    ic=0.5*(i0+i1);
    jc=0.5*(j0+j1);
    if ic==i0
      ic1=ic-0.5;
      ic2=ic+0.5;
      jc1=jc;
      jc2=jc;
    else
      ic1=ic;
      ic2=ic;
      jc1=jc-0.5;
      jc2=jc+0.5;
    end
    
    UTS.sgmX(cc,:)=[ic1,ic2];
    UTS.sgmY(cc,:)=[jc1,jc2];
    UTS.Ictr(cc)=is;
%  ==============    
    
% Corner points:
% contour changed direction, U/V still on same side
% corner point - check both components to close the contour
%  Example of an "open" corner point:
%  |________________|
%  |                |
%  |                |
%  -       * p-point| 
%  |                |
%  |                |
%  |                |
%  |_______|________|___________
%  |       ^        |
%  | en    | pn-vect|   pn-vector indicates direction following the contour p-pnts
%  - <--   * p-point|  <--* p-pnt contour
%  |       |        |     |  
%  |       |        |     | en-vector pointing towards V-norm component
%  |       v en     |     v  on the u.v contour
%  |                |        here only 1 component needed to close contour
%  |                |
%  |_______|________|_____|______  u/v pnt contour
%         V comp          V comp
%
% Note direction of pn (clockwise or cntr/clckwise) is unimportant
% for detecting the segments and norms, however
% chosen direction probably should be consistent for 1 contour
%
  elseif en_pn==en_pnL & pn_pnL~=0
    % Check u-comp if needed
    % if u-comp at (i0,j0) points to a contour p-point - then not needed
    % closed corner
    i1=i0-1;
    j1=j0;
    du=min(sqrt((I-i1).^2+(J-j1).^2));
    % if v-comp at (i0,j0) points to a contour p-point - then not needed
    % closed corner
    i1=i0;
    j1=j0-1;
    dv=min(sqrt((I-i1).^2+(J-j1).^2));
% Sanity check: this cannot be, violation of a closed corner point
    if (du==0 & dv>0) | (du>0 & dv==0)
      fprintf('***ERR: collocate_UTS_sections: incosistency corner pnts\n');
      keyboard;
    end
    
    if du==0 & dv==0 % "closed" corner point, no U/V transport 
      continue;
    end
    
% For open corner points
% need to decide which side close first
% to make contour continuously connected
% with the previous segment
    if du>0 & dv>0  % both components are needed, "open" corner point
% U-normal
      j1=j0;
      i1=i0-1;
      CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);

%if HH(j0,i0)>1e10
%   fprintf('ok\n'); 
%   keyboard; 
%end
      
      ic1=0.5*(i1+i0);
      jc1=j0;
      if cc>1
        dsgm=min(sqrt((UTS.sgmX(cc,:)-ic1).^2+(UTS.sgmY(cc,:)-jc1).^2));
      else
	dsgm=0;
      end
      
      if dsgm<1
        cE=0;
      else
	cE=1;
      end
      cc=cc+1;
      

      In=find(CLC.dHn<1e-3);
      CLC.Tn(In)=nan;
      CLC.Sn(In)=nan;
      CLC.Un(In)=nan;
  
      UTS.Hb(cc+cE)=-CLC.Hn;
      UTS.Unrm(:,cc+cE)=CLC.Un;
      UTS.Tnrm(:,cc+cE)=CLC.Tn;
      UTS.Snrm(:,cc+cE)=CLC.Sn;
      UTS.dH(:,cc+cE)=CLC.dHn;

% ===================
% Checking:
      ic=0.5*(i0+i1);
      jc=0.5*(j0+j1);
      if ic==i0
	ic1=ic-0.5;
	ic2=ic+0.5;
	jc1=jc;
	jc2=jc;
      else
	ic1=ic;
	ic2=ic;
	jc1=jc-0.5;
	jc2=jc+0.5;
      end

      UTS.sgmX(cc+cE,:)=[ic1,ic2];
      UTS.sgmY(cc+cE,:)=[jc1,jc2];
      UTS.Ictr(cc+cE)=is;
% ===================
      
% V-normal      
      j1=j0-1;
      i1=i0;
      CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
      
      cE=-cE; % adjust cc slot if needed, either 0 or -1
      
      cc=cc+1;

      In=find(CLC.dHn<1e-3);
      CLC.Tn(In)=nan;
      CLC.Sn(In)=nan;
      CLC.Un(In)=nan;
  
      UTS.Hb(cc+cE)=-CLC.Hn;
      UTS.Unrm(:,cc+cE)=CLC.Un;
      UTS.Tnrm(:,cc+cE)=CLC.Tn;
      UTS.Snrm(:,cc+cE)=CLC.Sn;
      UTS.dH(:,cc+cE)=CLC.dHn;

% ===================
% Checking:
      ic=0.5*(i0+i1);
      jc=0.5*(j0+j1);
      if ic==i0
	ic1=ic-0.5;
	ic2=ic+0.5;
	jc1=jc;
	jc2=jc;
      else
	ic1=ic;
	ic2=ic;
	jc1=jc-0.5;
	jc2=jc+0.5;
      end

      UTS.sgmX(cc+cE,:)=[ic1,ic2];
      UTS.sgmY(cc+cE,:)=[jc1,jc2];
      UTS.Ictr(cc+cE)=is;
      
% ===================
    end % if du>0
    
      
% Complicated case when the contour crosses
% the line where U/V transport is actually estimated
% This is reflected in change of orientation of en-vector wrt pn
% i.e. normal component from left-side switches to right-sight
% wrt to the local contour direction
% In this case choose normal component that 
% points along the contour (either backward or fwd)
  elseif en_pn~=en_pnL & pn_pnL~=0
    i1=i0-1;
    j1=j0;
    du=min(sqrt((I-i1).^2+(J-j1).^2));
    i1=i0;
    j1=j0-1;
    dv=min(sqrt((I-i1).^2+(J-j1).^2));
% Sanity check: this cannot be:    
    if (du==0 & dv==0) | (du>0 & dv>0)
      fprintf('***ERR: collocate_UTS_sections: 2. Incosistency corner pnts\n');
      keyboard;
    end

    if du==0 % need u-normal 
      j1=j0;
      i1=i0-1;
      CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
    else   % v-normal
      j1=j0-1;
      i1=i0;
      CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
    end

    cc=cc+1;

    In=find(CLC.dHn<1e-3);
    CLC.Tn(In)=nan;
    CLC.Sn(In)=nan;
    CLC.Un(In)=nan;

    UTS.Hb(cc)=-CLC.Hn;
    UTS.Unrm(:,cc)=CLC.Un;
    UTS.Tnrm(:,cc)=CLC.Tn;
    UTS.Snrm(:,cc)=CLC.Sn;
    UTS.dH(:,cc)=CLC.dHn;

% ===================
% Checking:
    ic=0.5*(i0+i1);
    jc=0.5*(j0+j1);
    if ic==i0
      ic1=ic-0.5;
      ic2=ic+0.5;
      jc1=jc;
      jc2=jc;
    else
      ic1=ic;
      ic2=ic;
      jc1=jc-0.5;
      jc2=jc+0.5;
    end

    UTS.sgmX(cc,:)=[ic1,ic2];
    UTS.sgmY(cc,:)=[jc1,jc2];
    UTS.Ictr(cc)=is;
% ===================

  end  % if
  
end  % sect

% 
% Calculate segment distances
%keyboard
sgmx=UTS.sgmX;
sgmy=UTS.sgmY;
ncc=length(sgmx);
for ik=1:ncc
  i1=floor(sgmx(ik,1));
  i2=floor(sgmx(ik,2));
  j1=floor(sgmy(ik,1));
  j2=floor(sgmy(ik,2));
  ln1=LON(j1,i1);
  lt1=LAT(j1,i1);
  ln2=LON(j2,i2);
  lt2=LAT(j2,i2);
  dL=distance_spheric_coord(lt1,ln1,lt2,ln2);
  UTS.segmL(ik,1)=dL;
end


%
% Find norms if needed (for simple sections, 
% not needed and +/- X-Y orientatino is taken) 
% Typically, this is the case
% of some polygons/contours/ boxes
% Norms are
% directed inside or outisde the contour
% First, construct a contour where fluxes 
% are calculated
% Second, go around the contour and 
% determine the direction of the norm 
% the end of the normal vector should be inside /outisde
% of the contour wrt to the centroid
% 

lpnt = 0;
if nin>1
  [jP,iP]=ind2sub(size(HH),nin);
  lpnt=1;
end

  
sgmx=UTS.sgmX;
sgmy=UTS.sgmY;

if ~isempty(nin); 

  % Combine contour going through u/v points
  ns=length(sgmx);
  Iuv=[];
  Juv=[];
  x1=sgmx(1,1)-10*(sgmx(1,2)-sgmx(1,1)); % extend segment to guarantee norm is in at the endpoints
  y1=sgmy(1,1)-10*(sgmy(1,2)-sgmy(1,1));
  Iuv(1,1)=x1;
  Juv(1,1)=y1;
  for ik=2:ns-1
    x1=sgmx(ik,1);
    x2=sgmx(ik,2);
    y1=sgmy(ik,1);
    y2=sgmy(ik,2);
    Iuv=[Iuv;0.5*(x1+x2)];
    Juv=[Juv;0.5*(y1+y2)];
  end
  x1=sgmx(end,1)+10*(sgmx(1,end)-sgmx(1,end-1)); % extend segment to guarantee norm is in at the endpoints
  y1=sgmy(end,1)+10*(sgmy(1,end)-sgmy(1,end-1));
  Iuv(ns,1)=x1;
  Juv(ns,1)=y1;
  if nin>1, % for sections, create a polygon, inside - norm is >0
    Iuv(ns+1)=iP;
    Juv(ns+1)=jP;
  end
  
  %keyboard

%CC = Centroid([I,J]);
  for ik=1:ns
    x1=sgmx(ik,1);
    x2=sgmx(ik,2);
    y1=sgmy(ik,1);
    y2=sgmy(ik,2);
    xs0=0.5*(x1+x2); % midpoint
    ys0=0.5*(y1+y2);
  % Rotate by 90 degrees:
  % around mid point
    aa=x2-xs0;
    bb=y2-ys0;
    rx=-bb;
    ry=aa;
    ll=sqrt(rx*rx+ry*ry);
    nx=rx./ll;
    ny=ry./ll;
    rxx=xs0+0.7*nx; % shorten norm to make it stay within the grid cell
    ryy=ys0+0.7*ny;

    [inn,onn]=inpolygon(rxx,ryy,Iuv,Juv);
    if onn, inn=~inn; end % point on the edge, happens at corner points
    if abs(nin==1)
      if (inn & nin==1) | ...
	 (~inn & nin==-1)
	UTS.Norm(ik,1)=nx;
	UTS.Norm(ik,2)=ny;
      else
	UTS.Norm(ik,1)=-nx;
	UTS.Norm(ik,2)=-ny;
      end
    else
      if inn 
	UTS.Norm(ik,1)=nx;
	UTS.Norm(ik,2)=ny;
      else
	UTS.Norm(ik,1)=-nx;
	UTS.Norm(ik,2)=-ny;
      end
    end
    
    
  %  keyboard
  end
  
  


end % nin



f_chck=0;
if f_chck==1
  fprintf('Plotting: Checking contour and u/v norms\n');
  figure(10); clf;
  hold on;
  contour(HH,[0 0],'k');
  axis('equal');
%  set(gca,'xlim',[800 1150],...
%	  'ylim',[800 1100]);
  set(gca,'xlim',[min(sgmx(:,1))-20 max(sgmx(:,1))+20],...
	  'ylim',[min(sgmy(:,1))-20 max(sgmy(:,1))+20]);
  plot(I,J,'.-');
  for ik=1:cc
    x=UTS.sgmX(ik,:);
    y=UTS.sgmY(ik,:);
    plot(x,y,'r-');

%
% Plot norms if defined
    if ~isempty(nin)
      nx=UTS.Norm(ik,1);
      ny=UTS.Norm(ik,2);

      xs0=mean(x); % midpoint
      ys0=mean(y);

      rxx=xs0+nx;
      ryy=ys0+ny;
      plot([xs0 rxx],[ys0 ryy],'Color',[0 1 0.5]);
      
    end
  
  end
  
  bottom_text(btx,'pwd',1);
  
  keyboard
  
  
end

% CHeck segments


  
return
