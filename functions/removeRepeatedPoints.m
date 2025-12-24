function output = removeRepeatedPoints(pts)

% make sure that there no repeated consecutive points. Sometimes "unique"
% command fails to find repeated points

threshold = 1e-10 ;

if iscell(pts)
    numSets = numel(pts) ;

    output = cell(numSets,1) ;

    for i = 1:numSets
        tmpOutput = pts{i} ;

        diffPts = vecmag(diff(pts{i})) ;
        idxSame = find(diffPts<threshold) ;
        tmpOutput(idxSame,:) = [] ;
        output{i} = tmpOutput ;
    end
else
    output = pts ;

    if size(pts,1)>1
        diffPts = vecmag(diff(pts)) ;
        idxSame = find(diffPts<threshold) ;
        output(idxSame,:) = [] ;
    end
end