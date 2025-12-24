function [F_bite,F_ws,F_bs,M] = biteForce(F_muscle,p_muscle,p_bite,p_TMJ)

% This functions calculates the bite and joint reaction forces assuming
% that the bite force is directed perpendicular to the vector between the
% bite point and the TMJ axis.

% calculate the angle between the horizontal and the vector connecting the
% bite point and the TMJ axis
alpha = -atan2d(p_bite(3),p_bite(1)) ;

% Determine the balancing and working side
if sign(p_TMJ(1,2))==sign(p_bite(2))
    p_bs = p_TMJ(2,:) ;
    p_ws = p_TMJ(1,:) ;
else
    p_bs = p_TMJ(1,:) ;
    p_ws = p_TMJ(2,:) ;
end

% For all calculations, we'll assume that the origin is at the balancing
% side TMJ
r_ws = p_ws - p_bs ;
r_muscle = p_muscle - p_bs ;
r_bite = p_bite - p_bs ;

% calculate the muscle moments around the balancing side TMJ
M = cross(r_muscle,F_muscle,2) ;
nMuscles = size(M,1) ;

%-------------------------------------------------------------------------%
% Calculation of forces if the bite force is perpendicular to the axis of
% the bite point
Fb = -M(:,2)./( r_bite(3).*sind(alpha)-r_bite(1).*cosd(alpha) ) ;
Fws_z = (-M(:,1) - r_bite(2).*Fb.*cosd(alpha)) ./ r_ws(2) ;
Fws_x = (M(:,3) - r_bite(2).*Fb.*sind(alpha)) ./ r_ws(2) ;
Fbs_x = -F_muscle(:,1) - Fws_x - Fb.*sind(alpha) ;
Fbs_z = -F_muscle(:,3) - Fws_z - Fb.*cosd(alpha) ;

F_bite = [Fb.*sind(alpha) zeros(nMuscles,1) Fb.*cosd(alpha)] ;
F_ws = [Fws_x zeros(nMuscles,1) Fws_z] ;
F_bs = [Fbs_x zeros(nMuscles,1) Fbs_z] ;