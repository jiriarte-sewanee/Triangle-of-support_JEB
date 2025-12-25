function [mMoments, mLength] = muscleLenghtMoment(musclePoints,refPoint)

% MUSCLELENGTHMOMENT calculates the length of the muscle-tendon unit (MTU)
% and the muscle moment arms with respect to the REFPOINT. MUSCLEPOINTS is
% a nPts*(3*nSegments) matrix of the cartesian coordinates of the muscle
% segments, where nPts indicates the number of points per segments and
% nSegments indicates the number of segments.

nSegments = size(musclePoints,2)/3 ;

if mod(nSegments,1)~=0
    error('Muscle position data should be multiple of 3 columns')
end

%% calculate muscle moments

% determine the direction of the each individual muscle segment
directionMuscle = -diff(musclePoints,1,1) ;

% for each muscle fiber, find the attachment point that has no change in
% muscle direction, and use this as the point as the point where muscle
% force will be applied to
musclePosition     = nan(nSegments,3) ;
muscleDirection = nan(nSegments,3) ;
for i = 1:nSegments
    % calculate the differences in direction between segments
    diffDir = ang2vec(directionMuscle(1:end-1,3*i-2:3*i),...
        directionMuscle(2:end,3*i-2:3*i),'deg') ;
    
    % determine the segments that do not change direction
    isSameAngle = diffDir<0.1 ;
    
    % find the last point that does not change direction. This will be the
    % point where moment will be calculated
    idxPoint = find(isSameAngle,1,'last') ;
    
    % determine the position of the muscle point
    selPoint = musclePoints(idxPoint+2,3*i-2:3*i) ;

    if isempty(selPoint)
        musclePosition(i,:) = musclePoints(end,3*i-2:3*i) ;
        muscleDirection(i,:) = normalizeVector3d(directionMuscle(end,3*i-2:3*i)) ;
    else

        musclePosition(i,:) = selPoint ;

        % determine the direction of the muscle fiber at this point
        selDir = directionMuscle(idxPoint+1,3*i-2:3*i) ;
        selDir = selDir./norm(selDir) ;
        muscleDirection(i,:) = selDir ;
    end
end

% determine the vector between the muscle position and the reference point
rVector = musclePosition - refPoint ;

% calculate the muscle moments. Because the muscle force used here is a
% unit force, we are calculating the moment arm instead of the moment per se.
% The units of the moment is just "m" (or whatever units are
% used in the model) instead of "Nm"
mMoments = cross(rVector,muscleDirection) ;

%% calculate the muscle lenghts
distBetweenSegments = vecmag(directionMuscle,3) ;
mLength = sum(distBetweenSegments) ;