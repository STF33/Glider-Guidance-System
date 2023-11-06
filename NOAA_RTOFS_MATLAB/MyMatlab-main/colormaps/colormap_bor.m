  function cmp = colormap_bor;
% Blue - orange - red colormap
% 15 color shades
% for colormaps that have different
% # of intervals below and above 0, 
% e.g. skewed t anomalies

cmp=[];
p=255;
cmp(1,:)     = [25,25,112];
cmp(end+1,:) = [0,0,255];
cmp(end+1,:) = [65,105,225];
cmp(end+1,:) = [100,149,255];
cmp(end+1,:) = [202,225,255]; % lightsteelblue
cmp(end+1,:) = [255,255,224];
cmp(end+1,:) = [255,250,205];
cmp(end+1,:) = [255,246,143];
cmp(end+1,:) = [255,255,0];
cmp(end+1,:) = [200,200,0];

cmp(end+1,:) = [255,218,185];
cmp(end+1,:) = [255,125,96];
cmp(end+1,:) = [255,127,0];
cmp(end+1,:) = [255,50,0];
cmp(end+1,:) = [205,20,0];
cmp=cmp./p;
