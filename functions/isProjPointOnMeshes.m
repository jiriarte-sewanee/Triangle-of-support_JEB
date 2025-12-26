function [ptsOnMesh,ptsIdx] = isProjPointOnMeshes ...
(TR, refPts, faceFraction, fixedAxis)

% ISPROJPOINTONMESHES finds the closest vertices of the mesh TR to the
% points (REFPTS) projected along the axis FIXAXIS. FACEFRACTION indicates
% the fraction of the mean face size to be used as the critical distance to
% define if a point close enough to a vertex.


% get the size of the search width as a fraction of the maximum face size
[~,sizeTR] = faceSize(TR.ConnectivityList,TR.Points) ;
w = sizeTR*faceFraction ;

% initialize output
nRefPts = size(refPts,1) ;
ptsOnMesh = cell(nRefPts,1) ;
ptsIdx    = cell(nRefPts,1) ;

ii = 1:3 ;
ii(fixedAxis) = [] ;

pts = TR.Points ;
for i = 1:nRefPts
    testPt = refPts(i,:) ;

    x = testPt(ii(1)) ; y = testPt(ii(2)) ;

    isGood = pts(:,ii(1)) > x-w & pts(:,ii(1)) < x+w & ...
        pts(:,ii(2)) > y-w & pts(:,ii(2)) < y+w ;
    selIdx = find(isGood) ;
    selPts = TR.Points(selIdx,:) ;
    ptsOnMesh{i} = selPts ;
    ptsIdx{i}    = selIdx ;
end
