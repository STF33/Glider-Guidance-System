  function cmp = cmp1;
%
% discrete colormap 
% from grey - blue - green - red - brown
% 15 colors

  clear cmp
  cmp(1,:)     =[0.9, 0.9, 0.9];
  cmp(end+1,:) =[0.7, 0.7, 0.7];
  cmp(end+1,:) =[0.5, 0.5, 0.5];
  cmp(end+1,:) =[0.35, 0, 0.5];
  cmp(end+1,:) =[0, 0, 0.7];
  cmp(end+1,:) =[0, 0.7, 1];      % sky blue
  cmp(end+1,:) =[0, 1, 1];        % Light blue
  cmp(end+1,:) =[0, 0.8, 0];      % Green
  cmp(end+1,:) =[0.5, 1, 0.25];  % Yellow-green
  cmp(end+1,:) =[1, 1, 0];        % Yellow
  cmp(end+1,:) =[1, 0.7, 0];        % Orange
  cmp(end+1,:) =[1, 0.45, 0];        % Orange/Red
  cmp(end+1,:) =[1, 0, 0];        % Red
  cmp(end+1,:) =[0.6, 0., 0.];  % Brown
  cmp(end+1,:) =[0.4, 0.2, 0.1];  % Brown
