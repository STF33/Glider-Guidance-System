% Example how to plot RMSE for perturbation experiments 1a 
% HYCOM vs NEMO
% There are 7 experiments within 1a (started with a 7-day shift)
% and there are 2 Time periods when the experiments are initialized: May/Jun 2011
% and Jan/Feb 2012

% Experiments 1a

% Load grid/ bathymetry, indices of grid points inside the
% GoM and deep ocean <200m
flmat = 'HYCOM_grid_Iocn0200m.mat';
load(flmat);

LAT = GRID.LAT;
LON = GRID.LON;
HH  = GRID.HH;
Iocn= GRID.Iocn;

for iFcst=10:16
  for itime=1:2
    flnm = sprintf('RMSE0200m_fcst%2.2i-%2.2i01.mat',iFcst,itime);
    fprintf('Loading %s\n',flnm);
    load(flnm);
    
% Note that saved is (RMSE)^2 - squarred error
    RMSE     = sqrt(RMSERR.ERR_squared);
    RMSEprst = sqrt(RMSERR.ERRprst_squared); % persistance

% 2D RMSE:
    A=HH*nan;
    A(Iocn)=RMSE;

% To plot RMSE:
				figure(1); clf;
				axes('Position',[0.09 0.22 0.86 0.7]);
				hold on;
				%pcolor(LON,LAT,Lmsk); shading flat;
				%colormap(lcmp);
				%caxis([0 1]);
				%freezeColors;

				axis('equal');
				set(gca,'tickdir','out',...
								'xlim',[-98 -80],...
								'xtick',[-100:2:-70],...
								'ylim',[18 31],...
								'ytick',[18:34],...
								'Fontsize',12);

				pcolor(LON,LAT,A); shading flat;
				%caxis([c1 c2]);
				%colormap(cmp);

				contour(LON,LAT,HH,[-200 -200],'k-','Color',[0.4 0.4 0.4],...
												'Linewidth',1.6);

				clb=colorbar('SouthOutside');
				set(clb,'Position',[0.18 0.1 0.66 0.025],...
												'Fontsize',12,...
												'Ticklength',0.025);

  end
end



