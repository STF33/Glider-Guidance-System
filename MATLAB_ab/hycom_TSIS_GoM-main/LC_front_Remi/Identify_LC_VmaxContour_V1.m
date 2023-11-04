function [LC_struct]  = Identify_LC_VmaxContour_V1(Longitude_Vec,Latitude_Vec,SSH,u,v,SSH_resolution,YC,FS)
% [LC_struct]  = Identify_LC_VmaxContour_V1(Longitude_Vec,Latitude_Vec,SSH,u,v,SSH_resolution,YC,FS)
%
%  Identify the Loop Current (LC) as the SSH contour associated with the
%  maximum averaged curvilinear velocity
%
%  Input Parameters :
%  X and Y column vectors
%
%  SSH,u,v in m and m/s. size(SSH,1)=size(X,2) and size(SSH,2)=size(Y,2)
%
%  SSH_resolution (can be left empty) = SSH resolution parameter to define
%  the LC. The smaller it is, the higher the accuracy will be, but the longer 
%  the run time will be.
%
%  YC and FS  (can be left empty) = sections of the Yucatan and Florida straits.
%  2x2 matrices containing the longitude (1st column) and latitude (2nd
%  column) of the two ends of the section
%
%  Return the structure LC_struct that contains longitude, latitude vectors
%  and SSH value of the contour. It is empty if no LC contour has been identified
%  (i.e. no SSH contour cross both YC ans FS)
%
%
%  FSU COAPS, Rémi Laxenaire and Dmitry Dukhovskoy

%%%%%%% Parameters
if isempty(YC)
    YC=[-87, 21.5; -84, 22]; % X and Y cordinates of the Yucatan Channel
end

if isempty(FS)
    FS=[-81, 25.5; -81,23.];% X and Y cordinates of the Straits of Florida
end

if isempty(SSH_resolution)
    SSH_resolution = 1e-3; % SSH resolution to define a contour [m]
end

%%%%%%% Define the check sections in the straits:

%%% Yucatan:
D=distance_spheric_coord(YC(1,2),YC(1,1),Latitude_Vec,Longitude_Vec);
[~,id_min]=min(D);
[j1,i1] = ind2sub(size(D),id_min(1));
IJY(1,1)=i1(1);
IJY(1,2)=j1(1);

D=distance_spheric_coord(YC(2,2),YC(2,1),Latitude_Vec,Longitude_Vec);
[~,id_min]=min(D);
[j1,i1] = ind2sub(size(D),id_min(1));
IJY(2,1)=i1(1);
IJY(2,2)=j1(1);

% Get Yucatan section:
jy0=IJY(1,2);
iy0=IJY(1,1);
iy1=IJY(2,1);
if iy0>iy1; dmm=iy1; iy1=iy0; iy0=dmm; end


%%% Str. Florida
D=distance_spheric_coord(FS(1,2),FS(1,1),Latitude_Vec,Longitude_Vec);
[~,id_min]=min(D);
[j1,i1] = ind2sub(size(D),id_min(1));
IJF(1,1)=i1(1);
IJF(1,2)=j1(1);
D=distance_spheric_coord(FS(2,2),FS(2,1),Latitude_Vec,Longitude_Vec);
[~,id_min]=min(D);
[j1,i1] = ind2sub(size(D),id_min(1));
IJF(2,1)=i1(1);
IJF(2,2)=j1(1);

% Get Str. Florida section:
jf0=IJF(1,2);
if0=IJF(1,1);
jf1=IJF(2,2);
if jf0>jf1; dmm=jf1; jf1=jf0; jf0=dmm; end;

%%%%%%% Identify innermost contour (i.e. lowest SSH contour connecting
% the two sections)

% Define the limits to scan
Outmost= double(nanmax([nanmin(SSH(jy0,iy0:iy1)),nanmin(SSH(jf0:jf1,if0))]))-SSH_resolution;
Inmost =  double(nanmin([nanmax(SSH(jy0,iy0:iy1)),nanmax(SSH(jf0:jf1,if0))]))+SSH_resolution;
Inmost_tmp = Inmost;

% Dichotomic research of the innermost contour
Bisol=NaN;
while Inmost-Outmost>SSH_resolution
    
    kk=1;
    
    cc1 = contourc(Longitude_Vec(1,:),Latitude_Vec(:,1), SSH, [Inmost Inmost]);
    np=size(cc1,2);
    
    LCfound=0;
    while kk<np & LCfound==0
        nrd=cc1(2,kk);
        iiS=kk+1;
        iiE=kk+nrd;
        kk=iiE+1;
        
        xx=cc1(1,iiS:iiE)';
        yy=cc1(2,iiS:iiE)';
        
        if length(xx)>1
            [X0,Y0,i0,~] = my_intersections(yy,xx,YC(:,2),YC(:,1));
            
            if i0<length(xx)/2
                i0=1;
            else
                i0=length(xx);
            end
            if length(X0)==1 & yy(i0)<Y0
                [X0,~] = my_intersections(yy,xx,FS(:,2),FS(:,1));
                if length(X0)==1
                    
                    LCfound=1;
                    Bisol = Inmost;
                    
                end
            end
        end
        
    end; % while kk<np
    
    % Condition on the detection of contour corresponding to criterions
    if LCfound == 1 % One contour was detected
        
        % Divide and shift the area of research toward exterior
        % (toward fin). ie Divide distance between contour and
        % exterior by 2. Reset IsGood to 0.
        Outmost = Inmost;
        
    elseif LCfound == 0 % No contour was detected
        
        % Divide and shift the area of research toward interior
        % (toward center). ie Divide distance between contour and
        % center by 2.
        
        Inmost_tmp = Inmost;
    end
    
    Inmost = Outmost + (Inmost_tmp-Outmost)/2;
end
Inmost_All=Bisol;


%%%%%%% Identify outermost contour (i.e. highest SSH contour connecting
% the two sections)

% Define the limits to scan
Outmost= double(nanmax([nanmin(SSH(jy0,iy0:iy1)),nanmin(SSH(jf0:jf1,if0))]))-SSH_resolution;
Inmost =  Inmost_All;
Inmost_tmp = Inmost;

% Dichotomic research of the outermost contour
Bisol=NaN;
while Inmost-Outmost>SSH_resolution
    
    kk=1;
    
    cc1 = contourc(Longitude_Vec(1,:),Latitude_Vec(:,1), SSH, [Inmost Inmost]);
    np=size(cc1,2);
    
    LCfound=0;
    while kk<np & LCfound==0
        nrd=cc1(2,kk);
        iiS=kk+1;
        iiE=kk+nrd;
        kk=iiE+1;
        
        xx=cc1(1,iiS:iiE)';
        yy=cc1(2,iiS:iiE)';
        
        if length(xx)>1
            [X0,Y0,i0,~] = my_intersections(yy,xx,YC(:,2),YC(:,1));
            
            if i0<length(xx)/2
                i0=1;
            else
                i0=length(xx);
            end
            if length(X0)==1 & yy(i0)<Y0
                [X0,~] = my_intersections(yy,xx,FS(:,2),FS(:,1));
                if length(X0)==1
                    
                    LCfound=1;
                    Bisol = Inmost;
                    
                end
            end
        end
        
    end; % while kk<np
    
    % Condition on the detection of contour corresponding to criterions
    if LCfound == 1 % One contour was detected
        
        % Divide and shift the area of research toward interior
        % (toward center). ie Divide distance between contour and
        % center by 2.
        Inmost_tmp = Inmost;
        
    elseif LCfound == 0 % No contour was detected
        
        
        % Divide and shift the area of research toward exterior
        % (toward fin). ie Divide distance between contour and
        % exterior by 2. Reset IsGood to 0.
        Outmost = Inmost;
    end
    
    Inmost = Outmost + (Inmost_tmp-Outmost)/2;
end
Outmost_All=Bisol;

% If none contour was found, exit the function
if isnan(Outmost_All) |  isnan(Inmost_All)
    Bisol_Flag=[NaN,0];
    LC_struct=struct;
    Id_IN_LC = NaN(size(SSH));
    return
end

%%%%%%% First scan of up to 20 calculated speeds on the contours between 
%%%%%%% the identified terminals
Line_tmp=linspace(Outmost_All,Inmost_All,20);
dp=max([Line_tmp(2)-Line_tmp(1),SSH_resolution*2]);
pos_all=Outmost_All:dp:Inmost_All;
vit_all=zeros(size(pos_all));
for loop_tmp=1:length(vit_all)
    
    kk=1;
    cc1 = contourc(Longitude_Vec(1,:),Latitude_Vec(:,1), SSH, [pos_all(loop_tmp) pos_all(loop_tmp)]);
    np=size(cc1,2);
    
    LCfound=0;
    while kk<np & LCfound==0
        nrd=cc1(2,kk);
        iiS=kk+1;
        iiE=kk+nrd;
        kk=iiE+1;
        
        xx=cc1(1,iiS:iiE)';
        yy=cc1(2,iiS:iiE)';
        
        if length(xx)>1
            [X0,~] = my_intersections(yy,xx,YC(:,2),YC(:,1));
            if length(X0)==1
                [X0,~] = my_intersections(yy,xx,FS(:,2),FS(:,1));
                if length(X0)==1
                    
                    LCfound=1;
                    
                    v_tmp = interp2(Longitude_Vec,Latitude_Vec,v,xx,yy,'linear');
                    u_tmp = interp2(Longitude_Vec,Latitude_Vec,u,xx,yy,'linear');
                    
                    dy=diff(yy).*6370e3*pi/180;
                    dx=diff(xx).*6370e3*pi/180.*cosd((yy(2:end)+yy(1:end-1))/2);
                    a12=atand(dy./dx);
                    a12(dx>0 & dy>=0)=90-a12(dx>0 & dy>=0);
                    a12(dx<0 & dy>=0)=270-a12(dx<0 & dy>=0);
                    a12(dx<0 & dy<=0)=270-a12(dx<0 & dy<=0);
                    a12(dx>0 & dy<0)=90-a12(dx>0 & dy<0);
                    vit_all(loop_tmp)=nanmean(abs(cosd(a12).*v_tmp(1:end-1)+sind(a12).*u_tmp(1:end-1)));
                    
                end
            end
        end
        
    end; 
    
end
[~,id_max] = nanmax(movmean(vit_all,2));

%%%%%%% Second scan of up to 50 calculated speeds on the contours close to
%%%%%%% the maximum identified during the first scan
Line_tmp=linspace(pos_all(max([id_max-2,1])),pos_all(min([id_max+2,length(pos_all)])),50);
dp=max([Line_tmp(2)-Line_tmp(1),SSH_resolution]);
pos_all=Line_tmp(1):dp:Line_tmp(end);
vit_all=nan(size(pos_all));
for loop_tmp=1:length(vit_all)
    
    kk=1;
    cc1 = contourc(Longitude_Vec(1,:),Latitude_Vec(:,1), SSH, [pos_all(loop_tmp) pos_all(loop_tmp)]);
    np=size(cc1,2);
    
    LCfound=0;
    while kk<np & LCfound==0
        nrd=cc1(2,kk);
        iiS=kk+1;
        iiE=kk+nrd;
        kk=iiE+1;
        
        xx=cc1(1,iiS:iiE)';
        yy=cc1(2,iiS:iiE)';
        
        if length(xx)>1
            [X0,~] = my_intersections(yy,xx,YC(:,2),YC(:,1));
            if length(X0)==1
                [X0,~] = my_intersections(yy,xx,FS(:,2),FS(:,1));
                if length(X0)==1
                    
                    LCfound=1;
                    
                    v_tmp = interp2(Longitude_Vec,Latitude_Vec,v,xx,yy,'linear');
                    u_tmp = interp2(Longitude_Vec,Latitude_Vec,u,xx,yy,'linear');
                    % If curvilinear
                    dy=diff(yy).*6370e3*pi/180;
                    dx=diff(xx).*6370e3*pi/180.*cosd((yy(2:end)+yy(1:end-1))/2);
                    a12=atand(dy./dx);
                    a12(dx>0 & dy>=0)=90-a12(dx>0 & dy>=0);
                    a12(dx<0 & dy>=0)=270-a12(dx<0 & dy>=0);
                    a12(dx<0 & dy<=0)=270-a12(dx<0 & dy<=0);
                    a12(dx>0 & dy<0)=90-a12(dx>0 & dy<0);
                    vit_all(loop_tmp)=nanmean(abs(cosd(a12).*v_tmp(1:end-1)+sind(a12).*u_tmp(1:end-1)));
                    
                    
                end
            end
        end
        
    end; % while kk<np
    
end
[~,id_max] = nanmax(vit_all);

%%%%%%% Store the LC as the contour associated with the maximum averaged speed
Bisol=pos_all(id_max);

kk=1;
cc1 = contourc(Longitude_Vec(1,:),Latitude_Vec(:,1), SSH, [Bisol Bisol]);
np=size(cc1,2);

LCfound=0;
while kk<np & LCfound==0
    nrd=cc1(2,kk);
    iiS=kk+1;
    iiE=kk+nrd;
    kk=iiE+1;
    
    xx=cc1(1,iiS:iiE)';
    yy=cc1(2,iiS:iiE)';
    
    if length(xx)>1
        [X0,~] = my_intersections(yy,xx,YC(:,2),YC(:,1));
        if length(X0)==1
            [X0,~] = my_intersections(yy,xx,FS(:,2),FS(:,1));
            if length(X0)==1
                LCfound=1;
            end
        end
    end 
end % while kk<np


if LCfound==1
    % Store Loop Current Informations
    LC_struct(1).xx=xx';
    LC_struct(1).yy=yy';
    LC_struct(1).SSH=Bisol;
else
    LC_struct=struct;
end

end

function dist = distance_spheric_coord(xla1,xlo1,xla2,xlo2)
%
%   function dist = distance_spheric_coord(LAT1,LON1,LAT2,LON2)
%   units = m
%  this procedure calculates the great-circle distance between two
%  geographical locations on a spheriod given it
%  lat-lon coordinates with its appropiate trigonometric
%  signs.
%  INPUT: xla1, xlo1 - first point coordinates (latitude, longitude)
%         xla2, xlo2 - second point
% all input coordinates are in DEGREES: latitude from 90 (N) to -90,
% longitudes: from -180 to 180 or 0 to 360,
% LAT2, LON2 can be either coordinates of 1 point or N points (array)
% in the latter case, distances from Pnt 1 (LAT1,LON1) to all pnts (LAT2,LON2)
% are calculated
%  OUTPUT - distance (in m)
%  R of the earth is taken 6371.0 km
%
%
%
% FSU COAPS, Dmitry Dukhovskoy
% Vincenty formula is used

if (abs(xla1)>90); error('Latitude 1 is > 90'); end;
if (abs(xla2)>90); error('Latitude 2 is > 90'); end;

R=6371.0e3;
cf=pi/180;
phi1=xla1*cf;
phi2=xla2*cf;
lmb1=xlo1*cf;
lmb2=xlo2*cf;
dphi=phi2-phi1;
dlmb=lmb2-lmb1;
%keyboard
% Central angle between 2 pnts:
dmm1=(cos(phi1).*sin(dlmb)).^2;
dmm2=(cos(phi2).*sin(phi1)-sin(phi2).*cos(phi1).*cos(dlmb)).^2;
dmm3=abs(sin(phi2).*sin(phi1)+cos(phi2).*cos(phi1).*cos(dlmb));

dsgm = atan(sqrt((dmm1+dmm2)./dmm3));

% The great-circle distance:
dist1 = R*dsgm;

% Another formula:
dmm = sin(phi1).*sin(phi2)+cos(phi1).*cos(phi2).*cos(lmb2-lmb1);
dmm(dmm>1)=1;
dist2 = acos(dmm).*R;

dist=0.5*(dist1+dist2);

%keyboard

end

function [x0,y0,iout,jout] = my_intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Input checks.
%error(narginchk(2,5,nargin))

% Adjustments when fewer than five arguments are supplied.
switch nargin
    case 2
        robust = true;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 3
        robust = x2;
        x2 = x1;
        y2 = y1;
        self_intersect = true;
    case 4
        robust = true;
        self_intersect = false;
    case 5
        self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
        length(x1) ~= length(y1)
    error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
        length(x2) ~= length(y2)
    error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
    repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
    repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
    repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
    repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
    repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
    repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
    repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Force i and j to be column vectors, even when their length is zero, i.e.,
% we want them to be 0-by-1 instead of 0-by-0.
i = reshape(i,[],1);
j = reshape(j,[],1);

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
    remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

if robust
    overlap = false(n,1);
    warning_state = warning('off','MATLAB:singularMatrix');
    % Use try-catch to guarantee original warning state is restored.
    try
        lastwarn('')
        for k = 1:n
            T(:,k) = AA(:,:,k)\B(:,k);
            [unused,last_warn] = lastwarn;
            lastwarn('')
            if strcmp(last_warn,'MATLAB:singularMatrix')
                % Force in_range(k) to be false.
                T(1,k) = NaN;
                % Determine if these segments overlap or are just parallel.
                overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
            end
        end
        warning(warning_state)
    catch err
        warning(warning_state)
        rethrow(err)
    end
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
    % For overlapping segment pairs the algorithm will return an
    % intersection point that is at the center of the overlapping region.
    if any(overlap)
        ia = i(overlap);
        ja = j(overlap);
        % set x0 and y0 to middle of overlapping region.
        T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
            min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
        T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
            min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
        selected = in_range | overlap;
    else
        selected = in_range;
    end
    xy0 = T(3:4,selected).';
    
    % Remove duplicate intersection points.
    [xy0,index] = unique(xy0,'rows');
    x0 = xy0(:,1);
    y0 = xy0(:,2);
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        sel_index = find(selected);
        sel = sel_index(index);
        iout = i(sel) + T(1,sel).';
        jout = j(sel) + T(2,sel).';
    end
else % non-robust option
    for k = 1:n
        [L,U] = lu(AA(:,:,k));
        T(:,k) = U\(L\B(:,k));
    end
    
    % Find where t1 and t2 are between 0 and 1 and return the corresponding
    % x0 and y0 values.
    in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
    x0 = T(3,in_range).';
    y0 = T(4,in_range).';
    
    % Compute how far along each line segment the intersections are.
    if nargout > 2
        iout = i(in_range) + T(1,in_range).';
        jout = j(in_range) + T(2,in_range).';
    end
end
end