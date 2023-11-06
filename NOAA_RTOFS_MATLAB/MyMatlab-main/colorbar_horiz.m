	function [az,axc]  = colorbar_horiz (jjs,cct,hght,lngth,mint,fsz,bxc,posc,mbx,aend)

% ==============================================================	
%	function [az,axc]  = colorbar_horiz (jjs,cct,hght,lngth,mint,fsz,bxc,posc,mbx,aend)
% The program draws a simple horizontal colorbar under the plot
% for given array of colors jjs
% number of intervals on the colorbar = length of jjs
% with no subintervals within each color
% Input:
%       jjs   - array of colors for colormap
%       cct   - array with intervals
%       hght  - height (thickness) of the colorbar
%			 lngth  - length of the colorbar
%      mint   - interval for labeling the colorbar:
%               mint=1: every interval will be labelled
%      fsz    - fontsize of the colorbar labels
%      bxc    - box color of the colorbar
%      posc  - position of the colorbar: [xcolb ycolb hcolb wcolb], 
%              default - placed just below the
%               current gca
%     mbx - intervals for contouring color boxes on the colorbar
%    aend=1 - end of colorbar is drawn as an arrow, first and last intervals are not 
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
if (isempty(hght));  hght=0.05; end;
if (isempty(lngth)); lngth=0.8; end;
if (isempty(mint));  mint=1;    end;
if (isempty(fsz));   fsz=8;     end;
if (isempty(bxc));   bxc='k';   end;
if (isempty(mbx));   mbx=1;     end;
if (isempty(aend));  aend=0; end;
%if isempty(aend) | mint~=mbx;    aend=0;    end;
	
% get gca and position of the main plot
  gca1=gca;
  pos=get(gca1,'position');
  hplt=get(gcf,'currentaxes');% handles for the plot axes 

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
		
		
% Plot the colorbar
mm1=length(jjs);
jjs=reshape(jjs,mm1,3);  
ddy=1/mm1; % delta Y
for jw=1:mm1
  ycrd(jw,1:5)=[0.,0.4,0.4,0.,0.];
end
xcrd(1,1:5)=[0.,0.,ddy,ddy,0.];
for jw=2:mm1+1;
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
if aend~=1
   for jp=1:mm1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp,:));
     set(ph(jp),'edgecolor','none');
   end
	 	
  for ll=1:mbx:mm1
    jp1=ll;
    jp2=ll+mbx-1;
    jp2=min(jp2,mm1);
    p1=plot([xcrd(jp1,1) xcrd(jp1,2)],[ycrd(jp1,1) ycrd(jp1,2)],'-k');
    p2=plot([xcrd(jp1,2) xcrd(jp2,3)],[ycrd(jp1,2) ycrd(jp2,3)],'-k');
    p3=plot([xcrd(jp2,3) xcrd(jp2,4)],[ycrd(jp2,3) ycrd(jp2,4)],'-k');
    p4=plot([xcrd(jp2,4) xcrd(jp1,5)],[ycrd(jp2,4) ycrd(jp1,5)],'-k');
    set([p1,p2,p3,p4],'Color',bxc);
%    keyboard
  end

else   % draw arrowheads at the ends of colorbar
   for jp=2:mm1-1
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jjs(jp,:));
     set(ph(jp),'edgecolor','none');
   end
	 	
  for ll=2:mbx:mm1-1
    jp1=ll;
    jp2=ll+mbx-1;
    jp2=min(jp2,mm1-1);
    p1=plot([xcrd(jp1,1) xcrd(jp1,2)],[ycrd(jp1,1) ycrd(jp1,2)],'-k');
    p2=plot([xcrd(jp1,2) xcrd(jp2,3)],[ycrd(jp1,2) ycrd(jp2,3)],'-k');
    p3=plot([xcrd(jp2,3) xcrd(jp2,4)],[ycrd(jp2,3) ycrd(jp2,4)],'-k');
    p4=plot([xcrd(jp2,4) xcrd(jp1,5)],[ycrd(jp2,4) ycrd(jp1,5)],'-k');
    set([p1,p2,p3,p4],'Color',bxc);
  end

% Arrowhead
    xhead(1,1)=min(xcrd(1,:));
    yhead(1,1)=0.5*(min(ycrd(1,:))+max(ycrd(1,:)));
    xhead(1,2)=max(xcrd(1,:));
    yhead(1,2)=min(ycrd(1,:));
    xhead(1,3)=max(xcrd(1,:));
    yhead(1,3)=max(ycrd(1,:));
    xhead(1,4)=min(xcrd(1,:));
    yhead(1,4)=yhead(1,1);

    xhead(2,1)=min(xcrd(mm1,:));
    yhead(2,1)=min(ycrd(mm1,:));
    xhead(2,2)=max(xcrd(mm1,:));
    yhead(2,2)=0.5*(min(ycrd(1,:))+max(ycrd(1,:)));
    xhead(2,3)=min(xcrd(mm1,:));
    yhead(2,3)=max(ycrd(mm1,:));
    xhead(2,4)=min(xcrd(mm1,:));
    yhead(2,4)=min(ycrd(mm1,:));
    
    pp1=patch(xhead(1,:),yhead(1,:),jjs(1,:));
    pp2=patch(xhead(2,:),yhead(2,:),jjs(length(jjs),:));
     set(pp1,'edgecolor',bxc);
     set(pp2,'edgecolor',bxc);
	
%    set(gca,'ylim',[yhead(1,3) yhead(2,3)]);

end;  % if aend

% ===================================
% Labels on the colorbar
% ===================================
ytxt=-0.5;  % x position of the text
lstrt=1;

if aend==1, lstrt=2;  end;   % start labeling from 1st interval
  for ll=lstrt:mint:mm1+1,
    jw=ll;
    xtxt=xcrd(jw,1);
    if (ll<=length(cct))
    text(xtxt,ytxt,num2str(cct(ll)),'fontsize',fsz,...
	'HorizontalAlignment','Center');
    end
  end;
  set(gca,'visible','off','box','on');
	
  az  = gca; % handle of the colorbar axes;
  axc = get(gcf,'currentaxes');  % axes of the colorbar
% ============
% Go back to original axes
% ============
  set(gcf,'currentaxes',hplt);





