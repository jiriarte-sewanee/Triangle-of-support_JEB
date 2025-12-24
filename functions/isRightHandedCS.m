function [isRight,newCS] = isRightHandedCS(CS)

% function that checks whether a coordinate system is right-handed or not.
% If not, it inverts the last axis to make it right-handed. Each row
% represents an axis.

newCS = CS ;

axis3 = cross(CS(1,:),CS(2,:)) ;

if sign(axis3)~=sign(CS(3,:))
    isRight = false ;
    newCS(3,:) = -CS(3,:) ;
else
    isRight = true ;
end