function sub_plot_dTz_nemo(pos,ssh,LONH,LATH,LON,LAT,...
                  c1,c2,stt,IN,xl1,xl2,yl1,yl2,clb,Lmsk,varargin);
% Plot dlt T at Z depth
%c1=-0.5;
%c2=0.5;

cmp = [];
nV = length(varargin);
if nV>0
  for k=1:nV
    vfld = varargin{k};
    if strmatch(vfld,'cmp') % 
      cmp = varargin{k+1}; %
    end
  end
end

if isempty(cmp);
		nc=100;
		clB=flipud(colormap_blue(nc));
		clR=colormap_red(nc);
		clB(nc-3,:)=[1 1 1];
		clB(nc-2,:)=[1 1 1];
		clB(nc-1,:)=[1 1 1];
		clB(nc,:)  =[1 1 1];
		clR(4,:)   =[1 1 1];
		clR(3,:)   =[1 1 1];
		clR(2,:)   =[1 1 1];
		clR(1,:)   =[1 1 1];
		cmp=[clB;clR];
		cmp=smooth_colormap(cmp,5,3);
end
nint=length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % 


cmpL=[0 0 0; 1 1 1];

axes('Position',pos);
pcolor(LONH,LATH,Lmsk); shading flat;
colormap(cmpL);
hold on;
freezeColors;

pcolor(LON,LAT,ssh); shading flat;
caxis([c1 c2]);
colormap(cmp);
title(stt,'Interpreter','none');

% ==================
% Plot dT
% =================
smm=ssh;
smm(IN==0)=nan;

if clb==1
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
	         'Position',[0.25 0.11 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
elseif clb==2
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
          'Position',[0.25 0.08 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
elseif clb==3
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
          'Position',[0.25 0.06 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
end
axis('equal');

%set(gca,'xlim',[-98 -56.08],'ylim',[7.05 31.92]);
%set(gca,'xlim',[-97.71 -78],'ylim',[15 31.1]);
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
set(gca,'tickdir','out',...
	'xtick',[-96:2:-50],...
	'ytick',[4:2:40],...
	'Fontsize',12);

return
