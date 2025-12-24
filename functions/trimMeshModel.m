function newTR = trimMeshModel(TR, verts2Remove)

% Modified from:
% REMOVEMESHVERTICES Remove vertices and associated faces from a mesh.
%
%   [V2, F2] = removeMeshVertices(VERTS, FACES, VERTINDS)
%   Removes the vertices specified by the vertex indices VERTINDS, and
%   remove the faces containing one of the removed vertices.
%
%   Example
%     % remove some vertices from a soccerball polyhedron
%     [v, f] = createSoccerBall; 
%     plane = createPlane([.6 0 0], [1 0 0]);
%     indAbove = find(~isBelowPlane(v, plane));
%     [v2, f2] = removeMeshVertices(v, f, indAbove);
%     drawMesh(v2, f2);
%     axis equal; hold on;
%
%   See also 
%     meshes3d, trimMesh

% ------
% Author: David Legland
% E-mail: david.legland@nantes.inra.fr
% Created: 2016-02-03, using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016-2022 INRA - Cepia Software Platform

% parse input
if any(size(verts2Remove)==1)
    if size(verts2Remove,1)==1 
        idx = verts2Remove ;
    else
        idx = verts2Remove' ;
    end
elseif size(verts2Remove,1)==3
    idx = ismember(TR.Points,verts2Remove,'rows') ;
else
    error('verts2remove must be either a Nx1 list of inidices or a Nx3 list of coordinates')
end

vertices = TR.Points ;
faces    = TR.ConnectivityList ;

% create array of indices to keep
nVertices = size(vertices, 1);
newInds = (1:nVertices)';
newInds(idx) = [];

% create new vertex array
newVertices = vertices(newInds, :);

% compute map from old indices to new indices
oldNewMap = zeros(nVertices, 1);
for iIndex = 1:size(newInds, 1)
   oldNewMap(newInds(iIndex)) = iIndex; 
end

% change labels of vertices referenced by faces
newFaces = oldNewMap(faces);
if size(newFaces,2)==1
    newFaces=newFaces'; 
end

% keep only faces with valid vertices
newFaces = newFaces(sum(newFaces == 0, 2) == 0, :);


% format output arguments
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
newTR = triangulation(newFaces,newVertices) ;
warning('on','MATLAB:triangulation:PtsNotInTriWarnId')