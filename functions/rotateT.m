function rotData = rotateT(data,T)

% Function that rotates and translates the input matrix 'data' using the
% homogeneous transformation matrix 'T'.
%  Input:   data:   XYZ position of the N markers to be rotated. Data must
%                   be entered as a nPoints X (3*nMarkers) matrix
%           T:      4*4*t transformation matrix obtained from the soder
%                   function. It represents the matrix [R t ; 0 0 0 1],
%                   where R is a 3x3 rotation matrix and P is 3x1
%                   translation matrix
%  Output:  XYZrot: Nx3*t matrix with the rotated and translated position of
%                   the markers

% get the number of transformation matrices
if size(T,3)>1
    nTranforms = size(T,3) ;
    fT = T ;
else
    if size(T)==4
        fT = T ;
        nTranforms = 1 ;
    elseif size(T,2)==16
        nTranforms = size(T,1) ;
        fT = reshape(T',4,4,nTranforms) ;
    elseif size(T)==3
        fT = T ;
        fT(4,4) = 1 ;
        nTranforms = 1 ;
    else
        error('T matrix is wrongly formatted')
    end
end

% get the number of frames and markers
nPoints  = size(data,1)   ;
nMarkers = size(data,2)/3 ;

if nTranforms==1
    newData = [reshape(data',3,nPoints*nMarkers) ; ones(1,nPoints*nMarkers)] ;
    rotNewData = (fT * newData) ;
    rotData = reshape(rotNewData(1:3,:),3*nMarkers,nPoints)' ;
else
    if nTranforms==nPoints
        % this assumes that every point will be rotated by their unique
        % rotation matrix
        rotData = nan(nPoints,3*nMarkers) ;
        for i = 1:nTranforms
            newData = [reshape(data(i,:),3,nMarkers) ; ones(1,nMarkers)] ;
            rotNewData = (fT(:,:,i) * newData) ;
            rotData(i,:) = reshape(rotNewData(1:3,:),1,3*nMarkers) ;
        end
    else
        % This assumes that the points are aligned in a single row
        if nPoints==1 && nMarkers>1
            newData = [reshape(data',3,nMarkers) ; ...
                ones(1,nMarkers)] ;
            rotData = nan(nTranforms,3*nMarkers) ;
            for i = 1:nTranforms
                rotNewData = (fT(:,:,i) * newData) ;
                rotData(i,:) = reshape(rotNewData(1:3,:),1,3*nMarkers) ;
            end
        % this assumes that data are in a single n*3 matrix
        elseif nMarkers==1 && nPoints>1
            newData = [data ones(nPoints,1)]' ;
            rotData = nan(nTranforms,3*nPoints) ;
            for i = 1:nTranforms
                rotNewData = (fT(:,:,i) * newData) ;
                rotData(i,:) = reshape(rotNewData(1:3,:),1,3*nPoints) ;
            end
        else
            rotData = nan(nPoints,3*nMarkers,nTranforms) ;
            for i = 1:nTranforms
                newData = [reshape(data',3,nPoints*nMarkers) ; ...
                    ones(1,nPoints*nMarkers)] ;
                rotNewData = (fT(:,:,i) * newData) ;
                rotData(:,:,i) = reshape(rotNewData(1:3,:),nPoints,3*nMarkers) ;
            end
        end
    end
end