function intPts = smoothPts(pts,varargin)

% Function to smooth 3D points.
% a smoothing parameter of 1- (3e-7) is a good value. For a smoother
% parameter, increase the value (bring farther from 1)

if nargin==1
    smoothParameter = varargin{1} ;
end

[nRows,nCols] = size(pts) ;

x = linspace(0,1,nRows)' ;

intPts = nan(size(pts)) ;

for i = 1:nCols
    isGood = ~isnan(pts(:,i)) ;
    if sum(isGood)<6; continue; end
    x_in = x(isGood) ; y_in = pts(isGood,i) ;
    if nargin==1
        fittedCurve = fit(x_in,y_in,'smoothingspline',...
            'smoothingparam',smoothParameter) ;
    else
        fittedCurve = fit(x_in,y_in,'smoothingspline') ;
    end
    y_out = feval(fittedCurve,x_in) ;
    intPts(isGood,i) = y_out ;
end
