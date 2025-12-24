function [angleValue]=ang2vec(v1,v2,varargin)

%ANG2VEC Angle between two vectors.
%
%   THETA = vectorAngle3d(V1, V2)
%   Computes the angle between the two vectors V1 and V2. The result ANGLEVALUE
%   is given in either radians (between 0 and PI) or in degrees (between 0 and 180).

uv1 = v1 ./ repmat (vecmag(v1) , 1, size(v1,2));
uv2 = v2 ./ repmat (vecmag(v2) , 1, size(v2,2));

if isempty(varargin)
    angleType = 'rad' ; % default
else
    angleType = varargin{1} ;
end

switch angleType
    case 'rad'
        angleValue = 2 .* atan2 (vecmag(uv1-uv2) , vecmag(uv1+uv2));
    case 'deg'
        angleValue = 2 .* atan2d (vecmag(uv1-uv2) , vecmag(uv1+uv2));
end