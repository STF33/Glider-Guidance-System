function plot_map(HH);
%  Quickly plot map of the region
%
figure(20); clf;
hold on
contour(HH,[0 0],'k');
contour(HH,[-500 -500],'b');



return