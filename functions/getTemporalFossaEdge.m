function fossaPts = getTemporalFossaEdge(TR_mandible,TR_cranium,isRight) 

% This function finds a series of points on the posterior edge of the
% temporal fossa, which are needed to determine the proper position of the
% posterior fibers of the temporalis muscle.
%
% Inputs are:
%   TR_mandible & TR_cranium: triangulation representation of the mandible
%                             and cranium models, respectively.
%   isRight: flag that indicates whether the fossa is on the right side of
%            the model or not.

% Determine the side of the cranium 
if isRight
    cf = -1 ;
else
    cf = 1 ; 
end

%-------------------------------------------------------------------------%
% Get the position of the fossa
%-------------------------------------------------------------------------%
pts = TR_cranium.Points ;
sizeCranium = max(TR_mandible.Points) - min(TR_mandible.Points) ;
isInFrontX = pts(:,1)>0 & pts(:,1)<sizeCranium(1)/2 ;
isBelowZ   = pts(:,3)<max(pts(:,3)/4);
isLateral = pts(:,2).*cf > 0 ;

isSelected = isInFrontX & isBelowZ & isLateral ;    
crTR = trimMeshModel(TR_cranium,~isSelected) ;
pts = crTR.Points ;

% find the best cross section to calculate the lateral location of the
% temporal fossa
sizeLims = [min(pts);max(pts)] ;
nCuts = 20 ;
xTest = linspace(sizeLims(1,1),sizeLims(2,1),nCuts) ;

innerLimits = nan(nCuts,2) ;
for i = 2:nCuts-2
    x1 = xTest(i) ; x2 = xTest(i+1) ;
    is1 = pts(:,1)>=x1 & pts(:,1)<=x2 ;
    a = alphaShape(pts(is1,1),pts(is1,2)) ;
    a.Alpha = a.Alpha*2 ;
    
    if a.numRegions==2
        [~,aPts1] = boundaryFacets(a,1) ;
        [~,aPts2] = boundaryFacets(a,2) ;
        if abs(max(aPts1(:,2))) > abs(max(aPts2(:,2)))
            innerLimits(i,:) = [max(abs(aPts2(:,2))) min(abs(aPts1(:,2)))].*cf ;
        else
            innerLimits(i,:) = [max(abs(aPts1(:,2))) min(abs(aPts2(:,2)))].*cf ;
        end
    end
end

% find the most outside outer edge
[~,selIdxX] = max(abs(innerLimits(:,2))) ;
yLimIn  = innerLimits(selIdxX,1) ;
yLimOut = innerLimits(selIdxX,2) ;
maxX = xTest(selIdxX) ;

% trim the cranium model
pts3 = TR_cranium.Points ;

isBetweenX = pts3(:,1)>=0 & pts3(:,1)<maxX ;
isBelowZ = pts3(:,3) < sizeCranium(3)/5 ;
if isRight
    isBetweenY = pts3(:,2)<yLimIn*1.1 & pts3(:,2)>yLimOut*1.1 ;
else
    isBetweenY = pts3(:,2)>yLimIn*1.1 & pts3(:,2)<yLimOut*1.1 ;
end
isGood = isBetweenX & isBetweenY & isBelowZ ;
         
redTR = trimMeshModel(TR_cranium,~isGood) ;

%-------------------------------------------------------------------------%
% Get the position of the condyle
%-------------------------------------------------------------------------%

% get the mandible points
pts1 = TR_mandible.Points ;

% determine the mandible size
mdSize = max(pts1) - min(pts1) ;

% find the points close to where the condyle should be
isGood = pts1(:,3)>-mdSize(3)/10 & pts1(:,1)<mdSize(1)/20 & pts1(:,2).*cf>0 ;
condylePos = pts1(isGood,:) ;

% get the reference points
yRef = max(abs(condylePos(:,2)))*cf ;
xRef = mean(condylePos(:,1)) ;
zRef = -mdSize(2)/2 ;
refCondyle = [xRef yRef zRef] ;

% find the closest point on the cranium
idx = nearestNeighbor(redTR,refCondyle) ;
seedPos = redTR.Points(idx,:) ;

% split the different non-connected meshes and find the one closest to the
% seed position on the cranium. Eliminate all other meshes.
m = splitMesh(redTR.Points,redTR.ConnectivityList);

D = nan(size(m,1),1) ;
for i=1:size(m,1)
    m(i) = trimMesh(m(i)) ;
    if isempty(m(i).faces)
        D(i) = nan ;
        continue
    end
    pts2test = m(i).vertices ;
    dist = vecmag(pts2test - seedPos) ;
    D(i) = min(dist) ;
end
[~,ii] = min(D) ;

% get the selected mesh
selMesh = trimMesh(m(ii))  ;

testY = linspace(yLimIn,yLimOut,22) ;

tmpPts = nan(20,3) ;
for i = 2:21
    pp = intersectPlaneSurf(selMesh,[0 testY(i) 0],[0 1 0]) ;
    if isempty(pp)
        continue
    elseif numel(pp)>1
        N = nan(numel(pp),1) ;
        for j = 1:numel(pp)
            if any(pp{j}(1,:)>maxX)
                N(j) = nan ;
            else
                N(j) = size(pp{j},2) ;
            end
        end
        [~,jj] = max(N) ;
        pp = pp(jj) ;
    end
    points = cell2mat(pp)' ;
    
    % find the most posterior point and make it the first point
    [~,minIdx] = min(points(:,1)) ;
    points = points([minIdx:end 1:minIdx-1],:) ;

    if points(1,3)>points(end,3)
        points = flipud(points) ;
    end
    dirPts = diff(points) ;
    angle1 = atan2d(dirPts(:,1),abs(dirPts(:,3))) ;
    idx = find(diff(sign(angle1))) ;
    [~,idx1] = max(points(idx,1)) ;
    selIdx = idx(idx1) ;
    if isempty(selIdx) || selIdx==1 || any(points(:,1)>maxX); continue; end
    tmpPts(i-1,:) = points(selIdx+1,:) ;
end
tmpPts(isnan(tmpPts(:,1)),:) = [] ;
fossaPts = tmpPts ;

% find outliers and remove them, if needed
z = zscore(fossaPts(:,3)) ;
idx2 = find(abs(z)>1) ;
fossaPts(idx2,:) = [] ;