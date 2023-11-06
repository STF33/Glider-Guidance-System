function [az,ph]=draw_colorbar_vert(cct,jj,mint,fsz)      
% function [az,ph]=draw_colorbar_vert(cct,jj,mint,fsz)      
% draw vertical colorbar
% colorbar is placed in (0,0) - left lower corner, 
% width = 0.4, height = 1.
% input: cct - data interval of contouring (cct(1) - min value)
%        jj - colormap handle with color codes, 
%             colormap used for filling contour intervals
%        jj=colormap; gives jj
%        mint - step in labeling on colorbar, 
%        fsz - font size for colorbar labels

mm1=length(cct); % # of intervals on colorbar
%======================
% Determine coordinates of a 'box'
%  x1,y2 ---- x2,y2
%    |       |
%    |       |
%    |       |
%  x1,y1 ---- x2,y1
% 
%======================
 ddy=1/mm1; % delta Y
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
     ph(jp) = patch(xcrd(jp,:),ycrd(jp,:),jj(jp,:));
   end
 set(ph,'edgecolor',[0 0 0]);
% creating boxes
% for jp=1:mm1
%   plot([xcrd(jp,1) xcrd(jp,2)],[ycrd(jp,1) ycrd(jp,2)],'-k')
%   plot([xcrd(jp,2) xcrd(jp,3)],[ycrd(jp,2) ycrd(jp,3)],'-k')
%   plot([xcrd(jp,3) xcrd(jp,4)],[ycrd(jp,3) ycrd(jp,4)],'-k')
%   plot([xcrd(jp,4) xcrd(jp,5)],[ycrd(jp,4) ycrd(jp,5)],'-k')
% end


% ===================================
% Labels on the colorbar
% ===================================
 xtxt=0.45;
 ytxt=(ycrd(1,2)-ycrd(1,1))/3;
 text(xtxt,ytxt,['<',num2str(cct(2))],'fontsize',fsz);

 ytxt=ycrd(mm1,1)+(ycrd(mm1,2)-ycrd(mm1,1))*2/3;
 text(xtxt,ytxt,['> ',num2str(cct(mm1))],'fontsize',fsz);
 xtxt=0.5
 for jw=2:mint:mm1,
   ytxt=ycrd(jw,1);
   text(xtxt,ytxt,num2str(cct(jw)),'fontsize',fsz);
 end;
 set(gca,'visible','off','box','on');
 az= gca; % handle of the colorbar axes;


