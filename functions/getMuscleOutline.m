function selMusclePts = getMuscleOutline(TR, refPts,...
    muscle, isCranium, fossaPts)

% GETMUSCLEOUTLINE creates the outline of a muscle segment on the mesh TR
% between the REFPTS that wraps around the bone, if needed. REFPTS is a 2*3
% matrix that contains the cranium and mandible attachment of the muscle
% segment. MUSCLE is a character array that indicates the muscle being
% processed. ISCRANIUM is a logical flag that indicates whether the model
% corresponds to a cranium or not. FOSSAPTS is the cartesian coordinates of
% the posterior edge of the temporal fossa (only needed for the temporalis
% muscle).

%-------------------------------------------------------------------------%
% if we are working with the cranium attachment of the medial pterygoid,
% we'll just use cranium attachment point.
if ismember(muscle,{'MedPterygoid','Masseter'}) && isCranium 
    selMusclePts = refPts(1,:) ;
    return
end

isTemporalis = ismember({muscle},'Temporalis') ;
if isTemporalis && ~isCranium
    selMusclePts = refPts(1,:) ;
    return
end

%-------------------------------------------------------------------------%
% Find if the muscle axis crosses the bone

% define the reference points
originPt = refPts(1,:) ;
attachmentPt = refPts(2,:) ;

% if the point is on the right side, we need to correct the direction
if attachmentPt(2)<0
    cf = -1 ;
else
    cf = 1 ;
end

% Find if the axisRef crosses any bone. If not the muscle is a straight
% line
[~,refPts] = directionalRayTriangleIntersect...
    (originPt, attachmentPt, TR, 0) ; 

if isempty(refPts)
    selMusclePts = originPt ;
    return
end
exitPoint = refPts(end,:) ;

altOrigin = [] ;
if isTemporalis
    % Find if the muscle axis passes through inside the temporal fossa ---%

    % create the line of the muscle axis
    lineAxis = createLine3d(originPt,attachmentPt) ;

    % calculate the perpendicular distance of the fossa points to the line
    % and find the closest fossa point
    dist1 = distancePointLine3d(fossaPts,lineAxis) ;
    [~,minIdx] = min(dist1) ;
    selPt = fossaPts(minIdx,:) ;

    % Determine if the selected point is to the left of the line projected
    % on the X-Z plane. If the point is to the left, then the lines passes
    % below the fossa and passes outside of it.
    isBelow = isLeftOriented(selPt([1 3]),lineAxis([1 3 4 6])) ;
    if isBelow
        isOutsideFossa = true ;
    else
        isOutsideFossa = false ;
    end

    % calculate the most external point around midway between the origin
    % and exit point. This will be used a the target point from the edge of
    % the fossa
    line1 = createLine3d(originPt,exitPoint) ;
    pos1 = line3dPoint(line1,0.5) ;
    tmpPts = isProjPointOnMeshes(TR,pos1,0.5,2) ;
    [~,idx1] = max(tmpPts{1}(:,2).*cf) ;
    refPos2 = tmpPts{1}(idx1,:) ;

    % If the muscle axis doesn't pass through the temporal fossa, we need
    % move the exit point to the closes point of the muscle axis to the
    % fossa points
    [~,face_size] = faceSize(TR.ConnectivityList,TR.Points) ;
    if isOutsideFossa
        altOrigin = refPos2 ;
        v1 = attachmentPt - fossaPts ;
        angle1 = atan2d(v1(:,2),v1(:,1)) ;
        v2 = attachmentPt - refPos2 ;
        angle2 = atan2d(v2(2),v2(1)) ;
        dist = angle1 - angle2 ;
        [~,minIdx] = min(abs(dist)) ;
        exitPoint = fossaPts(minIdx,:) ;
    else
        % even if it passes through the fossa, if the exit point is close
        % to the fossa points, then still use an intermediate point
        dist2fossa = vecmag(fossaPts - exitPoint) ;
        
        if any(dist2fossa<face_size*2)
            altOrigin = refPos2 ;
        end
    end
end

%-------------------------------------------------------------------------%
% MUSCLE TRACING
%-------------------------------------------------------------------------%
muscleAxis = [originPt; exitPoint] ;
if ismember(muscle,'Temporalis')
    isTemporalis = true ;
else
    isTemporalis = false;
end
musclePts = traceMuscle2point(TR, muscleAxis, attachmentPt,isTemporalis,altOrigin) ;

%-------------------------------------------------------------------------%
% IDENTIFY CONCAVE SECTIONS
% (where the muscle is not attached to the bone)
%-------------------------------------------------------------------------%

% get the muscle position to test
testPts = [musclePts; attachmentPt] ;

removeLastPoint = true ;
if isTemporalis
    % find if any of the resultant points is close to the fossa. If so,
    % that is going to be the terminal point to calculate the convex hull.
    % Otherwise, the use the attachment point.
    [~,dist2fossa2] = findClosestPoints(musclePts,fossaPts,0) ;
    [~,minDistIdx] = min(dist2fossa2) ;
    if dist2fossa2(minDistIdx)<face_size/2
        testPts = musclePts(1:minDistIdx,:) ;
        removeLastPoint = false ;
    end
end

if size(testPts,1)<3
    selMusclePts = testPts ; 
    return
end

% rotate the points so that parallel to the XY plane
[~,pp] = lsplane(testPts) ;
rotPts = (pp'*testPts')' ;

% define a straight line between the origin and attachment points 
refLine = createLine(rotPts(1,1:2),rotPts(end,1:2)) ;

% determine which points are to either side of the line
whichSide = isLeftOriented(rotPts(:,1:2),refLine) ;

% determine which side of the corresponds the origin (0,0). This will be
% the reference value to determine concavity
refSign = isLeftOriented([0 0],refLine) ;
% refSign = whichSide(2) ;

% determine which points are on the "good size". If something is on the
% same side as the origin, that means that the point is in a concavity and
% we need to get rid of it. Also, we need to make sure that the origin and
% attachment points are not included.
isBad = whichSide==refSign ;
isBad(1) = 0 ; isBad(end) = 0 ;

rotPts(isBad,:) = [] ;
testPts(isBad,:) = [] ;
% redPoints = [rotPts(1,:); rotPts(isGood,:); rotPts(end,:)] ;
% testPts = [testPts(1,:); testPts(isGood,:); testPts(end,:)] ;

% Finally, use the convex hull to remove further concave points
if size(rotPts,1)<=4
    selMusclePts = testPts(1:end-1,:) ; 
    return
end
k = convhull(rotPts(:,1:2),"Simplify",true) ;
selIdx = unique(k) ;
if removeLastPoint
    selIdx(end,:) = [] ; 
end
selMusclePts = testPts(selIdx,:) ;
