    function A = trapezoid_area(a,b,c,d);
% In general, area A=0.5*(a+b)*h, a, b - lengths of parallel sides
% a < b
%  h - height
% for a trapezoid with 2 parallel sides of different lengths, 
% the area can be determined from the length of all sides:

s1=(-a+b+c+d);
s2=(a-b+c+d);
s3=(a-b+c-d);
s4=(a-b-c+d);
A = 0.25*(a+b)/(b-a)*sqrt(s1*s2*s3*s4);
