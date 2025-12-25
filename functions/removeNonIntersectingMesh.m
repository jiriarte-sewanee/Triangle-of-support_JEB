function [newTR,intersectPts] = removeNonIntersectingMesh(TR,originPt,targetPt)

% Only keeps continous meshes (TR) that intersect the ray defined by the
% ORIGINPT and the TARGETPT.


% find the different contiguous meshes
TR = cleanMesh(TR) ;
meshes = splitMesh(TR.Points,TR.ConnectivityList) ;
nMesh = size(meshes,1) ;

% remove any empty meshes
toRemove = false(1,nMesh) ;
for i = 1:nMesh
    if isempty(meshes(i).faces)
        toRemove(i) = true ;
    end
end
meshes(toRemove) = [] ;
nMesh = size(meshes,1) ;

% get info size of the triangulation
modelSize = mean(max(TR.Points)-min(TR.Points)) ;
if modelSize > 1 % then the model is in mm
    modelTol = 1 ;
else % then the model is in m
    modelTol = 1*1e-3 ;
end

% find the intersecting points
aDir = normalizeVector3d(targetPt-originPt) ;
[faceIdx,refPts] = directionalRayTriangleIntersect...
    (originPt, targetPt+(aDir*1e3), TR, 0) ;     

% make sure that no point beyond the target point is included
critXpos = targetPt(1) ;
isBeyond = refPts(:,1) > critXpos+modelTol;
refPts(isBeyond,:) = [] ;
faceIdx(isBeyond)  = [] ;

% also make sure to include the origin point
if min(vecmag(refPts - originPt))>1e-13 | isempty(faceIdx)
    refPts = [originPt; refPts] ;
    faceIdx = [nan; faceIdx] ;
end
% if min(vecmag(refPts - targetPt))>1e-13 | isempty(faceIdx)
%     refPts = [refPts; targetPt] ;
%     faceIdx = [faceIdx; nan] ;
% end


% find the identified faces or, alternatively, the closest three vertices
% to a reference point
nPts = size(refPts,1) ;
makeContact = false(nPts,nMesh) ;

for i = 1:nPts
    testPt = refPts(i,:) ;
    if i==1 % find the closest vertex
        D = nan(nMesh,1) ;
        for j = 1:nMesh
            dist = vecmag(meshes(j).vertices - testPt) ;
            D(j) = min(dist) ;
        end
        [minD,minDidx] = min(D) ;
        if minD<modelTol
             makeContact(i,minDidx) = true ;
        end
    elseif ~isnan(faceIdx(i))
        vertIdx = TR.ConnectivityList(faceIdx(i),:) ;
        vertPts = TR.Points(vertIdx,:) ;
        D = nan(nMesh,1) ;
        for j = 1:nMesh
            [~,dist] = findClosestPoints(vertPts,meshes(j).vertices,0) ;
            D(j) = mean(dist) ;
        end
        [minD,minDidx] = min(D) ;
        if minD<modelTol
            makeContact(i,minDidx) = true ;
        end
    end
end

% remove any mesh that doesn't intersect with a point
isBad = ~any(makeContact,1) ;
makeContact(:,isBad) = [] ;
meshes(isBad) = [] ;
nMeshes = size(meshes,1) ;

% find if there is any point that is not associated to a mesh. If that is
% the case, then the reference point is the selected point
singlePointIdx = find(all(makeContact==0,2)) ;

% find the indices of the first an last intersections of each mesh
intersIdx = []; %nan(2,nMeshes) ;
selPts = {} ; %cell(1,nMeshes) ; 
selMeshes = [] ;
ii=1 ;
for q = 1:nMeshes
    cIdx = find(makeContact(:,q)) ;

    % if only we want first entry and last exit
    if isempty(cIdx); continue; end
    intersIdx(:,q) = cIdx([1 end])' ;
    selPts{q} = refPts(cIdx([1 end]),:) ;
    selMeshes(q) = q ;
end

% remove the internal intersections 
selMeshIdx = [] ;
selPtsIdx = [] ;
idx = 1:size(intersIdx,2) ;
while ~isempty(intersIdx)
    [~,maxIdx] = max(intersIdx(2,:)) ;
    % selMeshIdx = [idx(maxIdx) selMeshIdx] ;
    selMeshIdx = [selMeshes(maxIdx) selMeshIdx] ;
    selPtsIdx = [intersIdx(:,maxIdx) selPtsIdx] ;

    % remove any intersection between the selected indices
    i1 = intersIdx(1,maxIdx); i2 = intersIdx(2,maxIdx) ;
    intersIdx(:,maxIdx) = [] ;
    idx(maxIdx) = [] ;
    selMeshes(maxIdx) = [] ;
    if isempty(intersIdx); break; end
    isBad = all(intersIdx>i1 & intersIdx<i2,1) ;
    intersIdx(:,isBad) = [] ;
    idx(isBad) = [] ;
    selMeshes(isBad) = [] ;
end

% make sure to add any single point without associated mesh to the list
for t = numel(singlePointIdx):-1:1
    idx1 = singlePointIdx(t) ;
    selCol = find(selPtsIdx(1,:)==idx1+1) ;
    if isempty(selCol) % then the selPoint is at the end
        selPtsIdx = [selPtsIdx repmat(idx1,2,1)] ;
        selMeshIdx = [selMeshIdx nan] ;
    elseif selCol==1 % then the selPoint is at the beginning
        selPtsIdx = [repmat(idx1,2,1) selPtsIdx] ;
        selMeshIdx = [nan selMeshIdx] ;
    else 
        selPtsIdx = [selPtsIdx(:,1:selCol-1) repmat(idx1,2,1) selPtsIdx(:,selCol:end)] ;
        selMeshIdx = [selMeshIdx(1:selCol-1) nan selMeshIdx(selCol:end)] ;
    end
end

% finalize the output variable
nOutputs = numel(selMeshIdx) ;
newTR         = cell(1,nOutputs) ;
intersectPts  = cell(1,nOutputs) ;
for w = 1:nOutputs
    if ~isnan(selMeshIdx(w))
        tmpTR = trimMesh(meshes(selMeshIdx(w))) ;
        newTR{w} = triangulation(tmpTR.faces,tmpTR.vertices) ;
        
    end
    
    outIdx = unique(selPtsIdx(:,w)) ;
    intersectPts{w} = refPts(outIdx,:) ;
end

