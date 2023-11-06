function C = Centroid(P)
% Centroid(P), P=[X,Y]    
% This function computes the centroid, C, of the polygon P.  The polygon P
% must be entered as a n x 2 matrix, where row k denotes the (x,y) 
% coordinates of the k-th vertex of P.


% Redefine P
if P(1,:) ~= P(end,:)
    P=[P; P(1,:)];
end

% Define variables
n = length(P);
A = zeros(n,1);
Sx = zeros(n,1);
Sy = zeros(n,1);

% Compute the area of P.
for k = 1:n-1,
    A(k,1) = P(k,1)*P(k+1,2)-P(k+1,1)*P(k,2);
end

poly_area = .5*sum(A);


% Compute the centroid.
for k = 1:n-1,
    Sx(k,1) = (P(k,1)+P(k+1,1))*(P(k,1)*P(k+1,2)-P(k+1,1)*P(k,2));
end

for k = 1:n-1,
    Sy(k,1) = (P(k,2)+P(k+1,2))*(P(k,1)*P(k+1,2)-P(k+1,1)*P(k,2));
end

x = 1/(6*(poly_area))*sum(Sx);
y = 1/(6*(poly_area))*sum(Sy);
C = [x,y];

end

