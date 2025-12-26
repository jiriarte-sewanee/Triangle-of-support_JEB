function musclePts = traceMuscle2point(TR, muscleAxis, attachmentPt, ...
    isTemporalis, altRefPos)

% get the points of the muscle axis
originPt = muscleAxis(1,:) ;
exitPoint = muscleAxis(2,:) ; 

% determine if the muscle is on the left or the right side
if attachmentPt(2)<0
    cf = -1 ; 
else
    cf = 1 ;
end

%-------------------------------------------------------------------------%
% DETERMINATION OF OVERALL MUSCLE DIRECTION
%-------------------------------------------------------------------------%
muscleAxisDir = normalizeVector3d(diff(muscleAxis)) ;

% add a lateral normal component
if ~isempty(altRefPos)
    rotAngle = 5 ;
else
    rotAngle = 0 ;
end
rX = rotx(-rotAngle,'deg') ;
altRefPos = rotateT(altRefPos,rX) ;
ptsNormal = repmat(rotateT([0 1 0].*cf,rX),2,1) ;
% ptsNormal = repmat([0 1 0].*cf,2,1) ;

% find the normal perpendicular to the muscle axis
tmpAxis = cross(repmat(muscleAxisDir,2,1),ptsNormal) ;
normal2axis = normalizeVector3d(cross(tmpAxis,repmat(muscleAxisDir,2,1))) ;

% calculate the average between the two normals
avgNormal = mean(normal2axis) ;

% create a rotation matrix that is oriented based on the muscle axis
% direction and the average normal
axis1 = muscleAxisDir ;
axis2 = normalizeVector3d(cross(axis1,avgNormal)) ;
axis3 = normalizeVector3d(cross(axis1,axis2)) ;
R = [axis1;axis2;axis3] ;

% rotate the model and the associated points
newTR         = rotatePatch(TR      ,R) ;
rotOrigin     = rotateT(originPt    ,R) ;
rotExit       = rotateT(exitPoint   ,R) ;
rotAttachment = rotateT(attachmentPt,R) ;
rotTarget2    = rotateT(altRefPos   ,R) ;

refY = rotOrigin(2) ;
refZ = rotOrigin(3) ;

%-------------------------------------------------------------------------%
% MODEL TRIMMIING
%-------------------------------------------------------------------------%
% calculate the maximum face size and use to create a minimum addition to
% limits
[face_size] = faceSize(newTR.ConnectivityList,newTR.Points) ;
step = face_size*4 ;

% determine the vertices of the model that are within the desired values
pts = newTR.Points ;
isBetweenY = pts(:,2)>refY-step & pts(:,2)<refY+step ;
isBelowZ   = pts(:,3)<refZ+face_size*3 ;
isGood = isBelowZ & isBetweenY ;

% create a trimmed model
rotTR = trimMeshModel(newTR,~isGood) ;

% find the intersecting points and associated meshes
[selTRs,intersectPts] = removeNonIntersectingMesh(rotTR,rotOrigin,rotExit) ;
if isempty(intersectPts)
    musclePts = originPt ;
    return
end

%-------------------------------------------------------------------------%
%   TRACE THE MUSCLES PORTIONS BETWEEN INTERSECTION POINTS 
%-------------------------------------------------------------------------%
nParts = numel(selTRs) ;
useAltTarget = false ;
musclePtsParts = cell(nParts,1) ;
for q = 1:nParts

    % get the selected intersection pts and associated mesh
    selPts = removeRepeatedPoints(intersectPts{q}) ;
    
    % if there is only a single intersection point, then there is nothing
    % to trace
    if size(selPts,1)==1
        musclePtsParts{q} = selPts ;
        continue;
    end

    % get the selected mesh
    selMesh = selTRs{q} ;

    if q==nParts && isTemporalis
        selPts = flipud(selPts) ;
        if ~isempty(altRefPos)
            useAltTarget = true ;
        end
    end

    % calculate the distance between origin and target
    totalDist = vecmag(diff(selPts)) ;

    % determine the number of reference points as a function of the face size
    nPts = round(totalDist/(face_size*3)) ;

    % if points are too close, then just use the intersection point
    if nPts<3; musclePtsParts{q} = selPts ; continue; end

    % create a reference line between the origin and target point. This line
    % will be used as reference for the location of new points
    seed   = selPts(1,:) ;
    target = selPts(2,:) ;
    if useAltTarget
        altTarget = rotTarget2 ;
        altRefLine = createLine3d(seed,altTarget) ;
    end
    refLine = createLine3d(seed,target) ;
    refPos = linspace(0,1,nPts) ;
    refPts = line3dPoint(refLine,refPos) ;
    
    % initiate the output matrix
    tmpPts        = nan(nPts,3) ;
    tmpPts(1,:)   = seed ;
    tmpPts(end,:) = target ;
    currentPt = seed ;
    startIteration = 2 ;

    onlyEdges = true ;
    growingPts = seed ;
    hasReachedAltPt = false ;
    while useAltTarget & currentPt(1)>0 & ~hasReachedAltPt
        oldDist = vecmag(currentPt - altTarget) ;
        newPoint = findNextPoint(selMesh, currentPt, altTarget, ...
            'recursions',3,'distance3d',onlyEdges) ;
        % newPoint = findAttachedVertices2(selMesh, currentPt, altTarget, ...
        %     altRefLine,'distance3d',3,'recursions',onlyEdges) ;
        if isempty(newPoint); hasReachedAltPt = true ; continue; end
        
        newDist = vecmag(newPoint - altTarget) ;
        if newDist>oldDist
            hasReachedAltPt = true ;
        else
            growingPts(end+1,:) = newPoint ;
            currentPt = newPoint ;
        end
    end

    if size(growingPts,1)>1
        [~,currentIteration] = min(abs(refPts(:,1) - currentPt(1))) ;
        tmpPts(1:currentIteration,:) = curvspace(growingPts,currentIteration) ;
        startIteration = currentIteration +1 ;
    end

    for i = startIteration:nPts-1        
        %-----------------------------------------------------------------%
        % Find the points of the vertices associated to the last selected
        % muscle point
        newPoint = findNextPoint(selMesh, currentPt, target, ...
            'referencePoint',refPts(i,:),'angle2d',onlyEdges) ;
        % newPoint = findAttachedVertices2(selMesh, currentPt, target, ...
        %     refLine,'angle2d',refPos(i),'referencePoint',onlyEdges) ;

        if isempty(newPoint); continue; end
        tmpPts(i,:) = newPoint ;
        currentPt = newPoint ;
    end
    tmpPts(isnan(tmpPts(:,1)),:) = [] ;
    
    % if the data is reversed, flip the calculated muscle points
    if q==nParts && isTemporalis
        tmpPts = flipud(tmpPts) ;
    end
    
    musclePtsParts{q} = tmpPts ;
end

% combine all different parts into a single muscle variable pts
rotMusclePts = cell2mat(musclePtsParts) ;


%% Remove or add points at the end ---------------------------------------%
% check if the new point intersects the bone. This will determine what do
% next. If it is empty, then we need to remove points until we find the
% end. If it is not empty, then we need to add more points, until the
% intersection is empty.

% first, we need to create another trimmed model, but that includes the
% attachment point as well
target = rotAttachment ;
refPoints = [rotMusclePts; target] ;
limits = [min(refPoints);max(refPoints)] ;

pts = newTR.Points ;
isWithin = pts>(limits(1,:)-step*2) &  pts<(limits(2,:)+step*2) ;
toRemove = any(~isWithin,2) ;

testTR = trimMeshModel(newTR,toRemove) ;

idx = directionalRayTriangleIntersect...
    (rotMusclePts(end,:), target, testTR, 0) ;

endFound = false ;
if isempty(idx) % We might need to REMOVE points
    while ~endFound && size(rotMusclePts,1)>1
        % check if the new point intersects the bone, otherwise, this is
        % the end
        idx = directionalRayTriangleIntersect...
            (rotMusclePts(end-1,:), target, testTR, 0) ;

        if isempty(idx)
            rotMusclePts(end,:) = [] ;
        else
            endFound = true ;
        end
    end
else % we need to ADD points
    while ~endFound
        % Find the points of the vertices associated to the last selected
        % muscle point
        seedPt = rotMusclePts(end,:) ;
        onlyEdges = true ;
        % refLine = createLine3d(seedPt,target) ;
        newPoint = findNextPoint(rotTR, seedPt, target,...
            'recursions', 3, 'angle3d', onlyEdges) ;
        % newPoint = findAttachedVertices2(rotTR,seedPt, target, refLine,...
        %     'angle3d',3,'recursions',onlyEdges) ;

        % check if the new point is closer than the previous point
        oldDist = vecmag(seedPt - target) ;
        newDist = vecmag(newPoint - target) ;

        if newDist>=oldDist
            endFound = true ;
        else
            rotMusclePts = [rotMusclePts; newPoint] ;
            idx = directionalRayTriangleIntersect...
                (newPoint, target, rotTR, 0, 'exclusive') ;
            if isempty(idx)
                endFound = true ;
            end
        end
    end
end

%% Bring the points back to their original coordinate system
% rotate the musclePts to the original coordinate system
tmpMusclePts = rotateT(rotMusclePts,R') ;

% smooth the muscles
smooth_musclePts = smoothPts(tmpMusclePts,1-(3e-4)) ;
if isnan(smooth_musclePts)
    smooth_musclePts = tmpMusclePts ;
end

% make sure that the points are on the surface
idx = nearestNeighbor(TR,smooth_musclePts) ;
musclePts = TR.Points(idx,:) ;