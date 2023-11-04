  function vectors(uu,vv,cf,beta,X,Y,col,lwd,unt);

%  function vectors(uu,vv,cf,beta,X,Y,col,lwd,unt)
% program plots vectors using X and Y components
% X,Y - meshgrid coordinates of the vectors 
% uu, vv - components
% cf - scale coefficient: if cf=10, then 1 m/s is shown in 10 units of X or Y
%      the larger cf, the smaller velocities can be seen
%      e.g., if X, and Y are in meters and model has a large domain
%      appropriate scale for wind could be cf=0.5e3
% beta - angle (degrees) between vector 
%        line and arrow head beam (arrow head width)
%        the smaller the angle, the narrower the head
% col - color of the vectors:
%       y - yellow
%       b - blue
%       m - magenta
%       c - cinamon
%       g - green
%       r - red
%       w - white
%       k  - black
%       l - light gray
%       d - dark gray
%       a - gray
% lwd - line width, default = 1. pnt
% if u or v is NaN vector is not plotted
% unt >0 - plot unit vector which magnitude = unt in the right upper corner
%          if unt = 0 or unt=[], no unit vector

%  hold on

% ======================
% Set colors:
% =========================
if (~isempty(col));
  if (col=='Y')|(col=='y'); v_col=[1,1,0]; end;
  if (col=='B')|(col=='b'); v_col=[0,0,1]; end;
  if (col=='M')|(col=='m'); v_col=[1,0,1]; end;
  if (col=='C')|(col=='c'); v_col=[0,1,1]; end;
  if (col=='G')|(col=='g'); v_col=[0,1,0]; end;
  if (col=='R')|(col=='r'); v_col=[1,0,0]; end;
  if (col=='W')|(col=='w'); v_col=[1,1,1]; end;
  if (col=='K')|(col=='k'); v_col=[0,0,0]; end;
  if (col=='A')|(col=='a'); v_col=[0.5,0.5,0.5]; end;
  if (col=='L')|(col=='l'); v_col=[0.8,0.8,0.8]; end;
  if (col=='D')|(col=='d'); v_col=[0.3,0.3,0.3]; end;
end;
if (isempty(col)); v_col=[0,0,0]; end; % default color is black
if (isempty(lwd)); lwd=1.; end; % default line width
				
  [mn,ks]=size(uu);	
  [mn2,ks2]=size(uu);	
%  if (mn~=mn2) | (ks~=ks2), error('U and V matrices must be same size'), end;
  if (size(uu)~=size(vv)), error('U and V matrices must be same size'), end;

for ii=1:mn,
  for kk=1:ks,
    if (~isnan(uu(ii,kk)))&(~isnan(vv(ii,kk))); 
      sp=(uu(ii,kk)^2+vv(ii,kk)^2)^0.5;	% speed
      alfa=asin(uu(ii,kk)/sp)*180/pi;		% vector angle from Y
      if (vv(ii,kk)<0)&(uu(ii,kk)<0), 
        alfa=-alfa+180; 
      end;
      if (vv(ii,kk)>0)&(uu(ii,kk)<0), 
        alfa=alfa+360; 
      end;
      if (vv(ii,kk)<0)&(uu(ii,kk)>0), 
        alfa=180-alfa; 
      end;
      u=cf*uu(ii,kk); v=cf*vv(ii,kk);
      var=cf*.2*sp;				% scaling of the arrow head
      dX2=var*sin((alfa-beta)*pi/180);		% arrow head coordinates
      dX3=var*sin((alfa+beta)*pi/180);		% 
      dY2=var*cos((alfa-beta)*pi/180);
      dY3=var*cos((alfa+beta)*pi/180);
      x2=X(ii,kk)+u-dX2;
      x3=X(ii,kk)+u-dX3;
      y2=Y(ii,kk)+v-dY2;
      y3=Y(ii,kk)+v-dY3;
%     pause
	
%p1=plot([X(ii,kk) X(ii,kk)+u],[Y(ii,kk) Y(ii,kk)+v],'k-'); % vector
%p22=plot([X(ii,kk)+u x2],[Y(ii,kk)+v y2],'k-');		% arrow head
%p3=plot([X(ii,kk)+u x3],[Y(ii,kk)+v y3],'k-');
%p2=plot([X(ii,kk)+u x2],[Y(ii,kk)+v y2],'k-');		% arrow head
%set([p1,p2,p3],'color',v_col);
%clear p1, p2, p3;
%keyboard		
      p1=plot([X(ii,kk) X(ii,kk)+u],[Y(ii,kk) Y(ii,kk)+v],'Color',v_col); %vector
      p2=plot([X(ii,kk)+u x2],[Y(ii,kk)+v y2],'Color',v_col);		% arrow head
      p3=plot([X(ii,kk)+u x3],[Y(ii,kk)+v y3],'Color',v_col);
 
      set([p1,p2,p3],'linewidth',lwd);
			
    end; % if not isnan
  end;
end;

if ~isempty(unt) & unt>0
  gca1=gca;
  pos=get(gca1,'position');
  hplt=get(gcf,'currentaxes');% handles for the plot axes 
  x1=pos(1)+pos(3)*0.8;
  y1=pos(2)+pos(4)*0.8;
  hhx=0.2*pos(3);
  hhy=0.2*pos(4);
  hx = axes('position',[x1 y1 hhx hhy]);
  xx1=(x1-pos(1))*max(max(X));
  xx2=(x1+hhx-pos(1))*max(max(X));
  yy1=(y1-pos(2))*max(max(Y));
  yy2=(y1+hhy-pos(2))*max(max(Y));
  set(gca,'xlim',[xx1 xx2],'ylim',[yy1 yy2]);
  set(gca,'xticklabel',' ','yticklabel',' ','xtick',[],'ytick',[]);
  lx1=xx1+0.05*(xx2-xx1);
  lx2=lx1+cf*unt;
  ly1=0.5*(yy1+yy2);
  ly2=ly1;
  carh=0.2;   % scaling of the arrowhead
  draw_arrow(lx1,lx2,ly1,ly2,carh,beta,v_col,lwd);
%  hold on;
  text(0.5*(lx1+lx2),ly1+0.05*ly1,[num2str(unt)]);
  axis('equal');
% ============
% Go back to original axes
% ============
      set(gcf,'currentaxes',hplt);
end;
%  hold off
