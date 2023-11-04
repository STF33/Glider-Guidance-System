  function y0 = lin_interp(x1, x2, y1, y2, x0)
% ===================================
% Linear interpolation between two points y1 and y2:
% Find value y0 correponding to x0, such that x1<=x0<=x2
% y1=f(x1), y2=f(x2)
% ===================================
  y0 = y1+(y2-y1)/(x2-x1)*(x0-x1);
