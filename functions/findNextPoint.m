function newPoint = findNextPoint (TR, originPt, targetPt, ...
    recursionCriteria, critValue, directionCriteria, onlyEdges)

% FINDNEXTPOINT finds the next point on the surface of the mesh TR, from
% the ORIGINPT towards the TARGETPT. The code has two parts. First it finds
% a number of attached vertices to the origin point and then it finds the
% closest of those points that points towards the target point accordind to
% a direction criteria.
%
% RECURSIONCRITERIA indicates the criterion to be used to determine how many
%   attached vertices to find based on the CRITVALUE. The options are:
%       'recursions':     number of recursions (for-loops) to do before
%                         stopping
%       'distance':       stops when an attached points reaches this 3D 
%                         distance from the reference points
%       'referencePoint': stops when it reaches any  vertices defined by the
%                         perpendicular projection of the referencePoint
%
% DIRECTIONCRITERIA indicates the criterion used to select the closest
%   point towards the target point. The options are:
%       'distance2d': uses the point that is closest to the target point
%                     based on the 2D distance in the XZ plane.
%       'distance3d': uses the point that is closest to the target point
%                     based on the 3D distanc.
%       'angle2d':    uses the point that has the smallest angle between
%                     the current point and the origin and target point in
%                     the XZ plane.
%       'angle3d':    uses the point that has the smallest 3D angle between
%                     the current point and the origin and target point.
%
% ONLYEDGES is a flag that indicates if we only want to keep the outer
% edges of each recursion (TRUE) or we want to keep all the attached points
% (FALSE)
% 
% OUTPUT:
%   NEWPOINT: 3D cartesian coordinates of the following point.

refLine = createLine3d(originPt,targetPt) ;

switch recursionCriteria
    case 'recursions'
        critRecursion = critValue ;
        critIdx = 1 ;
    case 'distance'
        critDistance  = critValue ;
        critIdx = 2 ;
    case 'referencePoint'
        refPt = critValue ;
        critIdx = 3 ;

        % find all the indices that match the criteria
        [~,refIdx] = isProjPointOnMeshes(TR,refPt,0.6,3) ;
        refIdx = refIdx{1} ;

        % keep only the vertices with the correct normals
        refNormal = [0 0 -1] ;
        normalFaces = vertexNormal(TR,refIdx) ;
        verticesNormalAngle = ang2vec(normalFaces,refNormal,'deg') ;
        isBad = verticesNormalAngle>100 ;
        refIdx(isBad) = [] ;

        if isempty(refIdx)
            newPoint = [] ;
            return
        end
        selRefIdx = refIdx ;
end

%-------------------------------------------------------------------------%
% find the neighboring vertices
%-------------------------------------------------------------------------%

% find the vertex of the last muscle point
selVertexIdx   = nearestNeighbor(TR,originPt) ;
oldVerticesIdx = selVertexIdx ;

% define the current point
currentPt = TR.Points(selVertexIdx,:) ;

% find the normal of the origin and target vertices
refNormal = [0 0 -1] ;

isDone = false ;
nRecursions = 0 ;

while ~isDone
    nRecursions = nRecursions + 1 ;

    % find the neighboring triangles of the vertex/vertices
    trianglesIdx = vertexAttachments(TR,selVertexIdx) ;
    trianglesIdx = unique(cell2mat(trianglesIdx')) ;

    % find the neighboring vertices of the selected triangles
    neighborsIdx = TR.ConnectivityList(trianglesIdx(:),:) ;

    % compile the selected vertices
    selVertexIdx = unique([selVertexIdx;neighborsIdx(:)]) ;

    % if OnlyEdges is selected, make sure you don't use previously selected
    % vertices
    if onlyEdges
        selVertexIdx = setdiff(selVertexIdx,oldVerticesIdx) ;
        oldVerticesIdx = [oldVerticesIdx; selVertexIdx] ;
    end

    % check if the selected points reach the criteria to terminate the loop
    if critIdx==2 % DISTANCE as criterion
        distPts = vecmag(TR.Points(selVertexIdx,:) - originPt) ;
        if distPts>=critDistance
            isDone = true ;
        end
    elseif critIdx==1 % RECURSIONS as criterion
        if nRecursions>=critRecursion
            % make sure that we only consider points that have normals in
            % the right direction
            normalFaces = vertexNormal(TR,selVertexIdx) ;
            verticesNormalAngle = ang2vec(normalFaces,refNormal,'deg') ;
            isBad = verticesNormalAngle>100 ;

            % make sure that the points are more posterior than the
            % current point
            currentPos = line3dPosition(currentPt,refLine) ;
            linePos    = line3dPosition(TR.Points(selVertexIdx,:),refLine) ;
            isMoreAnterior = linePos<currentPos ;

            if ~isempty(selVertexIdx(~isBad & ~isMoreAnterior))
                selVertexIdx(isBad & isMoreAnterior) = [] ;
                isDone = true ;
            else
                newPoint = [] ;
                return
            end
            
        end
    elseif critIdx==3 % REFERENCE POINT as criterion
        % check if at least one of the reference points has been reached
        selIdx = find(ismember(selVertexIdx,selRefIdx)) ;
        if ~isempty(selIdx)
            % find the closest point to the seed point
            selPts = TR.Points(selVertexIdx(selIdx),:) ;
            dist2seed = vecmag(selPts - originPt) ;
            [~,newPointIdx] = min(dist2seed) ;
            newPoint = selPts(newPointIdx,:) ;
            return
        end
    end
end

% find the 3D position of the selected vertices
selPoints = TR.Points(selVertexIdx,:) ;

if size(selPoints,1)==1
    newPoint = selPoints ;
    return
end
%-------------------------------------------------------------------------%
% Find the one point that matches the criteria
%-------------------------------------------------------------------------%

% % make sure that none of the new points is more posterior than the last
% % selected muscle point
% if ~isempty(refLine)
%     currentPos = line3dPosition(originPt,refLine) ;
%     linePos    = line3dPosition(selPoints,refLine) ;
%     isBadPt = linePos<currentPos ;
%     selPoints(isBadPt,:) = [] ;
% end

switch directionCriteria
    case 'distance3d'
        % find the closest point to the target
        distPts = vecmag(selPoints - targetPt) ;
        [~,minIdx] = min(distPts) ;
        newPoint = selPoints(minIdx,:) ;

    case 'distance2d' % this uses the distance in the sagittal (XZ) plane
        idx = [1 3] ;
        distPts = vecmag(selPoints(:,idx) - targetPt(idx)) ;
        [~,minIdx] = min(distPts) ;
        newPoint = selPoints(minIdx,:) ;

    case 'angle2d' % this uses the angle in the sagittal (XZ) plane
        idx = [1 2] ;

        % calculate the vectors from the seedPoint to the identified points
        vectors = createVector(selPoints,targetPt) ;

        % calculate the 3D angle between the vectors and the reference vector
        refVector = createVector(originPt,targetPt) ;
        angles = ang2vec(vectors(:,idx),refVector(idx),'deg') ;

        % find the point with the minimal angle between. This is the new
        % selected point
        [~,minIdx] = min(abs(angles)) ;
        newPoint = selPoints(minIdx,:) ;

    case 'angle3d'
        % calculate the vectors from the seedPoint to the identified points
        vectors = createVector(selPoints,targetPt) ;

        % calculate the 3D angle between the vectors and the reference vector
        refVector = createVector(originPt,targetPt) ;
        angles = rad2deg(vectorAngle3d(vectors,refVector)) ;

        % find the point with the minimal angle between. This is the new
        % selected point
        [~,minIdx] = min(abs(angles)) ;
        newPoint = selPoints(minIdx,:) ;
end
