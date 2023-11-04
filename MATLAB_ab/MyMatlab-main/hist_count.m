  function [N,xb] = hist_count(X,bin);
%
%  Histogram counts for vector X, counts the numer of values in X 
% that fall between the elements in the bin vector
% similar to matlab function histc  
% COunts are done for the segments bin(k)<= & <bin(k)
% Note that # of counts is 1 less than specified # of bin intervals
%
% Output: N - hist. counts, xb - middle of intervals (for plotting, finding mode, etc);
%
Xs = sort(X);
nl = length(bin);
dx = bin(2)-bin(1);
xb = bin(1:nl-1)+dx*0.5;

for ik=1:nl-1,
  x1=bin(ik);
  x2=bin(ik+1);
  I=find(Xs>=x1 & Xs<x2);
  N(ik)=length(I);
end;



