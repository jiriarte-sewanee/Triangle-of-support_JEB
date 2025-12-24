function [index,minDist] = findClosestPoints(queryPoints,allPoints,doWaitbar)

if nargin<3
    doWaitbar = false ;
end

% number of point in first input to process
np = size(queryPoints, 1);
nVars = size(queryPoints,2) ;

% allocate memory for result
index = zeros(np, 1);
minDist = zeros(np, 1);

if doWaitbar
    for i = 1:np
        textwaitbar(i, np,'Finding the closest points');

        % compute squared distance between current point and all point in array
        tmp = (allPoints - queryPoints(i,:)).^2 ;

        if nVars==2
            dist = sqrt(tmp(:,1) + tmp(:,2)) ;
        else
            dist = sqrt(tmp(:,1) + tmp(:,2) + tmp(:,3)) ;
        end
        
        % keep index of closest point
        [minDist(i), index(i)] = min(dist);   
    end
else
    for i = 1:np
        % compute squared distance between current point and all point in array
        tmp = (allPoints - queryPoints(i,:)).^2 ;

        if nVars==2
            dist = sqrt(tmp(:,1) + tmp(:,2)) ;
        else
            dist = sqrt(tmp(:,1) + tmp(:,2) + tmp(:,3)) ;
        end

        % keep index of closest point
        [minDist(i), index(i)] = min(dist);        
    end
end