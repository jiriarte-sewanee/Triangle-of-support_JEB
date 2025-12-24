function [dist2vertex, verticesPts] = closestPoint2Mesh(TR, points)
%DISTANCEPOINTMESH Shortest distance between a (3D) point and a triangle mesh.
%
%   DIST = closestPoint2Mesh(TR,POINTS)
%   Returns the shortest distance between the query point(s) POINTS and the
%   triangular mesh defined by the triangulation TR. POINTS is a NP-by-3
%   array, TR is created by the 'triangulation' command.
%
%   [DIST, PROJ] = closestPoint2Mesh(...)
%   Also returns the NP-by-3 projection of the query point(s) on the 
%   triangular mesh.
%

verticesIdx = nearestNeighbor(TR,points) ;
verticesPts = TR.Points(verticesIdx,:) ;

% calculate the distance of the query points to the closest vertices
dist2vertex = sqrt(sum((verticesPts - points).^2,2)) ;
