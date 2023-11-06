   function cmp = colormap_PHC;
% Create colormap similar to that used in PHC website
% There are 12 intervals
clear cmp
cmp(1,:) =[0, 0, 0];
cmp(2,:) =[0.35, 0, 0.5];
cmp(3,:) =[0, 0, 0.7];
cmp(4,:) =[0, 0.7, 1];      % sky blue
cmp(5,:) =[0, 1, 1];        % Light blue
cmp(6,:) =[0, 0.8, 0];      % Green
cmp(7,:) =[0.5, 1, 0.25];  % Yellow-green
cmp(8,:) =[1, 1, 0];        % Yellow
cmp(9,:) =[1, 0.7, 0];        % Orange
cmp(10,:)=[1, 0.45, 0];        % Orange/Red
cmp(11,:)=[1, 0, 0];        % Red
cmp(12,:)=[0.5, 0.3, 0.1];  % Brown
