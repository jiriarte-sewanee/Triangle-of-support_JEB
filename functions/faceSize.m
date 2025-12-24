function [meanLength, maxLength, minLength] = faceSize(faces,vertices) 

% this function calculates the average distance of the vertices to the
% circumcenter of each face. This will gives an idea of the average size
% the all faces

p1 = vertices(faces(:,1),:) ;
p2 = vertices(faces(:,2),:) ;
p3 = vertices(faces(:,3),:) ;

d1 = vecmag(p1-p2) ;
d2 = vecmag(p2-p3) ;
d3 = vecmag(p1-p3) ;

meanLength = mean([d1;d2;d3]) ;
maxLength = max([d1;d2;d3]) ;
minLength = min([d1;d2;d3]) ;
