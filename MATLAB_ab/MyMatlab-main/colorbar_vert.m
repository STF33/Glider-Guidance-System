	function [az,axc]  = colorbar_vert (jjs,cct,wdth,hght,mint,fsz,bxc,posc,mbx,aend)

% ==============================================================	
% [az,axc]  = colorbar_vert (jjs,cct,hght,lngth,mint,fsz,bxc,posc,mbx,aend)
% The program draws a simple horizontal colorbar under the plot
% for given array of colors jjs
% number of intervals on the colorbar = length of jjs
% with no subintervals within each color
% Input:
%  JJS   - array of colors for colormap
%  CCT   - array with intervals
%  WDTH - specified colorbar widht (default = 0.05 Normal Units)
%  HGHT - specified colorbar height (default = 0.95) as fraction
%          of the hight of the main figure
%  MINT   - interval for labeling the colorbar:
%               mint=1: every interval will be labelled
%  FSZ    - fontsize of the colorbar labels
%  BXC    - box color of the colorbar
%  POSC  - position of the colorbar: [xcolb ycolb hcolb wcolb], 
%              default - placed just to the right of the
%               current gca
%  MBX    - intervals for contouring color boxes on the colorbar
%  aend=1 - end of colorbar is drawn as an arrow, first and last intervals are not 
%           labeled
%           the option is available ONLY if mbx = mint
%           if not 1 or by default, the first interval is labeled and drawn as "box"
%
% Output: 
%      az     - gca of the colorbar
%      axc    - axes of the colorbar
% ====================================

% -----------------------------------
% If parameters are not specified
% use the default
% -----------------------------------
if (isempty(wdth));  wdth=.05;  end;
if (isempty(hght));  hght=.95;  end;
if (isempty(mint));  mint=1;    end;
if (isempty(fsz));   fsz=8;     end;
if (isempty(bxc));   bxc='k';   end;
if (isempty(mbx));   mbx=1;     end;
	
% get gca and position of the main plot
gca1=gca;
pos=get(gca1,'position');
hplt=get(gcf,'currentaxes');% handles for the plot axes 

% ----------------------------------
% Set off the main plot to the left 
% to place the colorbar 
% ----------------------------------
if (isempty(posc));
  set(gca1,'position',[pos(1),pos(2),pos(3)-1.15*wdth,pos(4)]);
  wcolb=wdth;
  hcolb=hght*pos(4);
  xcolb=pos(1)+pos(3)-.85*wdth;
  ycolb=pos(2)+(pos(4)-hcolb)/2;
else
  xcolb=posc(1);
  ycolb=posc(2);
  hcolb=posc(3);
  wcolb=posc(4);
end
hx = axes('position',[xcolb ycolb wcolb hcolb]);
set(gca,'xlim',[0 1],'ylim',[0 1])
		
		
% Plot the colorbar
mm1=length(jjs);
ddy=1/mm1; % delta Y

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
		
%keyboard				
% =============================================
% filling boxes with corresponding patches
% =============================================
hold on
%keyboard
if aend~=1
   for jp=1:mm1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp,:));
     set(ph(jp),'edgecolor','none');
   end
	 
  for ll=1:mbx:mm1
    jp1=ll;
    jp2=ll+mbx-1;
    jp2=min([mm1,jp2]);
    p1=plot([xcrd(jp1,1) xcrd(jp2,2)],[ycrd(jp1,1) ycrd(jp2,2)],'-k');
    p2=plot([xcrd(jp2,2) xcrd(jp2,3)],[ycrd(jp2,2) ycrd(jp2,3)],'-k');
    p3=plot([xcrd(jp2,3) xcrd(jp1,4)],[ycrd(jp2,3) ycrd(jp1,4)],'-k');
    p4=plot([xcrd(jp1,4) xcrd(jp1,5)],[ycrd(jp1,4) ycrd(jp1,5)],'-k');
    set([p1,p2,p3,p4],'Color',bxc);
  end

else   % draw arrowheads at the ends of colorbar
  for jp=mbx:mm1-mbx
    ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp,:));
    set(ph(jp),'edgecolor','none');
  end
	 	
  for ll=mbx:mbx:mm1-2*mbx
    jp1=ll;
    jp2=ll+mbx-1;
    jp2=min(jp2,mm1-1);
    p1=plot([xcrd(jp1,1) xcrd(jp2,2)],[ycrd(jp1,1) ycrd(jp2,2)],'-k');
    p2=plot([xcrd(jp2,2) xcrd(jp2,3)],[ycrd(jp2,2) ycrd(jp2,3)],'-k');
    p3=plot([xcrd(jp2,3) xcrd(jp1,4)],[ycrd(jp2,3) ycrd(jp1,4)],'-k');
    p4=plot([xcrd(jp1,4) xcrd(jp1,5)],[ycrd(jp1,4) ycrd(jp1,5)],'-k');
    set([p1,p2,p3,p4],'Color',bxc);
  end

% Arrowhead
    yhead(1,1)=min(ycrd(1,:));
    xhead(1,1)=0.5*(min(xcrd(1,:))+max(xcrd(1,:)));
    yhead(1,2)=max(ycrd(mbx-1,:));
    xhead(1,2)=min(xcrd(1,:));
    yhead(1,3)=max(ycrd(mbx-1,:));
    xhead(1,3)=max(xcrd(1,:));
    xhead(1,4)=xhead(1,1);
    yhead(1,4)=yhead(1,1);

    yhead(2,1)=min(ycrd(mm1-mbx,:));
    xhead(2,1)=min(xcrd(mm1,:));
    yhead(2,2)=max(ycrd(mm1,:));
    xhead(2,2)=0.5*(min(xcrd(1,:))+max(xcrd(1,:)));
    yhead(2,3)=min(ycrd(mm1-mbx,:));
    xhead(2,3)=max(xcrd(mm1,:));
    xhead(2,4)=xhead(2,1);
    yhead(2,4)=yhead(2,1);
    
    pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
    pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
     set(pp1,'edgecolor',bxc);
     set(pp2,'edgecolor',bxc);

end;   % if aend
      
% ===================================
% Labels on the colorbar
% ===================================
xtxt=0.5;  % x position of the text
lstrt=1;

if aend==1, lstrt=1+mbx;  end;   % start labeling from second interval
%keyboard
for ll=lstrt:mint:mm1,
  jw=ll;
  ytxt=ycrd(jw,1);
  text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz,...
	'HorizontalAlignment','Left');
end;
set(gca,'visible','off','box','on');
	
az  = gca; % handle of the colorbar axes;
axc = get(gcf,'currentaxes');  % axes of the colorbar
% ============
% Go back to original axes
% ============
set(gcf,'currentaxes',hplt);





