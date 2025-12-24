function newTR = cleanMesh(TR)
% CLEANMESH Reduces memory footprint of a polygonal mesh.
%
%   newTR = cleanMesh(TR)
%   Unreferenced vertices are removed.
%   Following functions are implemented only for numeric faces:
%       Duplicate vertices are removed.
%       Duplicate faces are removed.
%
%   Modified from trimmesh (from MatGeom)
% ------


% parse input data
if isobject(TR) % it is a triangulation
    isTriangulation = true ;
    vertices = TR. Points ;
    faces    = TR.ConnectivityList ;
else
    isTriangulation = false ;
    vertices = TR.vertices ;
    faces    = TR.faces ;
end

% Delete duplicate vertices
[tempVertices, ~, tempFaceVertexIdx] = unique(vertices, 'rows') ;
tempFaces = tempFaceVertexIdx(faces) ;

% Delete unindexed/unreferenced vertices
usedVertexIdx = ismember(1:length(tempVertices),unique(tempFaces(:))) ;
newVertexIdx = cumsum(usedVertexIdx) ;
faceVertexIdx = 1:length(tempVertices) ;
faceVertexIdx(usedVertexIdx) = newVertexIdx(usedVertexIdx) ;
faceVertexIdx(~usedVertexIdx) = nan ;
tempFaces2 = faceVertexIdx(tempFaces) ;
tempVertices2 = tempVertices(usedVertexIdx,:) ;

% Delete duplicate faces
[~, uniqueFaceIdx, ~] = unique(tempFaces2, 'rows') ;
duplicateFaceIdx=~ismember(1:size(tempFaces2,1),uniqueFaceIdx) ;
[vertices2, faces2] = removeMeshFaces(tempVertices2, tempFaces2, duplicateFaceIdx) ;

% format output arguments
if isTriangulation
    newTR = triangulation(faces2,vertices2) ;
else
    newTR.faces    = faces2 ;
    newTR.vertices = vertices2 ;
end
