function [selIntersectIdx,selIntersectPts] = directionalRayTriangleIntersect...
    (proxPoint,distPoint,TR,includeFisrtPoint)

% Wrapper for the rayTriangleIntersection_mlc function that calculates
%      ray/triangle intersections using the algorithm proposed BY Möller
%      and Trumbore (1997). The ray begins at PROXPOINT and points towards
%      DISTPOINT. The triangles are defined by the trinagulation TR.
%      INCLUDEFIRSTPOINT is a flag whether or not to include the PROXPOINT
%      in the results
% 
%   [selIntersectIdx,selIntersectPts] = directionalRayTriangleIntersect (...) 
%     Returns:
%     * selIntersectIdx - indices of the triangles intersected by the ray
%     * selIntersectPts - cartesian coordinates of the intersection points


% get the ray position and direction
rayPos = proxPoint ;
rayDir = normalizeVector3d(distPoint-rayPos) ;

% get the vertices of the triangulation
vertices = TR.Points ;
faces    = TR.ConnectivityList ;
t1 = vertices(faces(:,1),:) ;
t2 = vertices(faces(:,2),:) ;
t3 = vertices(faces(:,3),:) ;

% find the intersecting faces using Möller and Trumbore (1997) algorithm
[isIntersect,pts] = rayTriangleIntersection_mlc_mex (rayPos, rayDir,...
    t1, t2, t3) ;
intersectIdx = find(isIntersect) ;

if isempty(intersectIdx)
    selIntersectIdx = [] ;
    selIntersectPts = [] ;
    return
else
    % get the intersection points
    intersectPts = pts(intersectIdx,:) ;

    % remove any intersecting points that are NOT between the queary
    % point and the distal point
    posLine = line3dPosition(intersectPts,[rayPos rayDir]) ;
    distProx2Dist = line3dPosition(distPoint,[rayPos rayDir]) ;
    isOutofLine = posLine < -1e-10 | posLine > distProx2Dist ;

    intersectIdx(isOutofLine)   = [] ;
    intersectPts(isOutofLine,:) = [] ;
    posLine(isOutofLine) = [] ;

    if isempty(intersectIdx)
        selIntersectIdx = [] ;
        selIntersectPts = [] ;
        return
    end

    % make sure to order the points in increasing order
    [~,sortIdx] = sort(posLine) ;
    intersectPts = intersectPts(sortIdx,:) ;
    intersectIdx = intersectIdx(sortIdx) ;

    % Because we are using the inclusive method, there can be points that
    % are associated with multiple faces. Let's find the closest face, so
    % that there is only one per point
    uniquePts = removeRepeatedPoints(intersectPts) ;

    if size(uniquePts,1)~=size(intersectPts,1)
        new_intersectIdx = nan(size(uniquePts,1),1) ;

        for i = 1:size(uniquePts,1)
            refPoint = uniquePts(i,:) ;

            % calculate the distance of all points to the reference point
            dist2ref = vecmag(intersectPts - refPoint) ;

            % find the points that are very close to the reference point
            isSelPoint = find(dist2ref<1e-6) ;

            dist2pt = nan(numel(isSelPoint),1) ;
            for j = 1:numel(isSelPoint)
                selFace = intersectIdx(isSelPoint(j)) ;
                selPtsIdx = TR.ConnectivityList(selFace,:) ;
                selPts = TR.Points(selPtsIdx,:) ;
                tmpDist = vecmag(selPts - refPoint) ;
                dist2pt(j) = mean(tmpDist) ;
            end
            [~,selIdx] = min(dist2pt) ;
            new_intersectIdx(i) = intersectIdx(isSelPoint(selIdx)) ;
        end
    else
        new_intersectIdx = intersectIdx ;
    end

    % let's make sure that the points are ordered from proximal to
    % distal point
    distFromRef = vecmag(uniquePts - rayPos) ;
    [~,sortIdx] = sort(distFromRef) ;

    selIntersectPts = uniquePts(sortIdx,:) ;
    selIntersectIdx = new_intersectIdx(sortIdx,:) ;


    % check if we need to add the first point or remove it
    if includeFisrtPoint
        % make sure that the proximal point is included
        distPt = vecmag(selIntersectPts(1,:)-proxPoint) ;

        if distPt > 1e-10 
            selIntersectPts = [proxPoint; selIntersectPts] ;
            selIntersectIdx = [nan;selIntersectIdx] ;
        end
    else
        % make sure that the proximal point is NOT included
        distPt = vecmag(intersectPts(1,:)-proxPoint) ;

        if distPt < 1e-10 
            selIntersectPts(1,:) = [] ;
            selIntersectIdx(1)   = [] ;
        end
    end
end