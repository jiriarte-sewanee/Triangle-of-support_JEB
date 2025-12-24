function [fiberStrain,pennAngle] = calculateFiberStrainFromMuscleLength...
    (diffMuscleLength, FiberLength, pennationAngle)

% Script to calculate the muscle fiber strain from changes in muscle
% length.
%
% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% diffLenghtMuscle --> changes in muscle length (in any length units, but
%                      it must match the units for "FiberLength"
% FiberLenght      --> Length of the muscle fiber at rest
% pennationAngle   --> Pennation angle of the muscle at rest (in degrees)
%
% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fiberStrain: structure with the muscle fiber strains for all conditions
% pennAngle  : structure with the resultant pennation angles for all the
%              conditions indicated below
%
% Conditions:
% constantAngle --> muscle fiber strain assuming that the pennation angle
%                     does not change (but the muscle width does change)
% constantWidth --> muscle fiber strain assuming that the fiber width
%                     does not change (but pennation angle changes)
% variableWidth --> muscle fiber strain assuming that both pennation
%                     angle and muscle width change (but assumes that the
%                     sectional area of the muscle remains constant)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nPoints,nVars] = size(diffMuscleLength) ;

%-------------------------------------------------------------------------%
% CONSTANT PENNATION ANGLE 
e_constantAngle = diffMuscleLength./(FiberLength.*cosd(pennationAngle))+1 ;
pa_constantAngle = ones(nPoints,nVars).*pennationAngle ;

%-------------------------------------------------------------------------%
% CONSTANT MUSCLE WIDTH
e_constantWidth = sqrt( ((FiberLength.*cosd(pennationAngle)+diffMuscleLength).^2 ...
    + (FiberLength.*sind(pennationAngle)).^2)./(FiberLength.^2) ) ;
pa_constantWidth = acosd((FiberLength.*cosd(pennationAngle)...
    + diffMuscleLength)./(e_constantWidth.*FiberLength) ) ;

%-------------------------------------------------------------------------%
% VARIABLE MUSCLE WIDTH 
b_squared = (cosd(pennationAngle) + (diffMuscleLength./FiberLength)).^2 ;
e_variableWidth = sqrt(b_squared + (sind(pennationAngle.*2)).^2./(4.*b_squared)) ;
pa_variableWidth = acosd((FiberLength.*cosd(pennationAngle)...
    + diffMuscleLength)./(e_constantArea.*FiberLength) ) ;

%-------------------------------------------------------------------------%
% Create the output files
fiberStrain.constantAngle = e_constantAngle ;
fiberStrain.constantWidth = e_constantWidth ;
fiberStrain.variableWidth = e_variableWidth  ;

pennAngle.constantAngle = pa_constantAngle ;
pennAngle.constantWidth = pa_constantWidth ;
pennAngle.variableWidth = pa_variableWidth  ;