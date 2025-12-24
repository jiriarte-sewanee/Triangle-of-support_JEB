function magnitude = vecmag(x,nDimensions)

% function to calculate the magnitude of a vector
if nargin<2
    nDimensions = size(x,2) ;
end

if rem(size(x,2),nDimensions)
    disp('error. Wrong number of columns!')
    return
end

nVariables = size(x,2)/nDimensions ;
magnitude = nan(size(x,1),nVariables) ;

for i = 1:nVariables
    if nDimensions==3
        magnitude(:,i) = sqrt( sum (x(:,3*i-2:3*i).^2,2 ) )  ;
    elseif nDimensions==2
        magnitude(:,i) = sqrt( sum( x(:,2*i-1:2*i).^2,2 ) ) ;
    end
end