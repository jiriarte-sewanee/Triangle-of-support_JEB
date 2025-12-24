function [hits,intersectPts] = rayTriangleIntersection_mlc(origin,rayDirs,...
    p0,p1,p2)

% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).
% Copyright 2014 The MathWorks, Inc.
%
% TRIANGLERAYINTERSECTION_MLC Ray/triangle intersection.
%    TriangleRayIntersection(ORIG, DIR, VERT1, VERT2, VERT3)
%      calculates ray/triangle intersections using the algorithm proposed
%      BY Möller and Trumbore (1997). The ray starts at ORIGIN and points
%      toward RAYDIRS. The triangle is defined by vertix points: p0, p1,
%      p2. All input arrays are in Nx3 or 1x3 format, where N is number of
%      triangles or rays.
% 
%   [hit, intersectPts] = TriangleRayIntersection(...) 
%     Returns:
%     * hits         - boolean array of length N informing which line and
%                      triangle pair intersect
%     * intersectPts - carthesian coordinates of the intersection point


% initialize the output variables
hits = false( size(p0,1), size(rayDirs,1) ) ;
intersectPts = nan( size(p0,1), 3, size(rayDirs,1) ) ;

for i = 1:size(rayDirs,1)
    for j = 1:size(p0,1)

        [flag,pos] = rayTriangleIntersection_inner(...
            origin, rayDirs(i,:), p0(j,:), p1(j,:), p2(j,:)) ;

        hits(j,i) = flag ;
        intersectPts(j,:,i) = pos ;

    end
end
% numhits = sum(hits,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag,pos] = rayTriangleIntersection_inner (o, d, p0, p1, p2)

% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).

epsilon = -1e-5 ; % Inclusive boder. More intersections are found.
% epsilon = 1e-5 ; % Exclusive border. Fewer intersections are found.

e1 = p1-p0;
e2 = p2-p0;
q  = cross(d,e2);
a  = dot(e1,q); % determinant of the matrix M

if (a>-epsilon && a<epsilon)
    % the vector is parallel to the plane (the intersection is at infinity)
    flag = false;
    pos = nan(1,3) ;
    return
end

f = 1/a;
s = o-p0;
u = f*dot(s,q);

if (u<0.0)
    % the intersection is outside of the triangle
    flag = false;
    pos = nan(1,3) ;
    return
end

r = cross(s,e1);
v = f*dot(d,r);

if (v<0.0 || u+v>1.0)
    % the intersection is outside of the triangle
    flag = false ;
    pos = nan(1,3) ;
    return
end
flag = true ;

% calculate intersection coordinates
pos = p0 + e1.*repmat(u,1,3) + e2.*repmat(v,1,3);

return
end
