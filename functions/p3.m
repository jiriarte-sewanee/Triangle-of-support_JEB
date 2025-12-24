function varargout = p3(data,cmap,varargin)

% wrapper for the plot3 function.
% DATA is nPts x (3*nSets), where nPts in the total number of points and
%       nSets is the numbers of sets of points
% CMAP is the colormap used to plot multiple sets. If empty,it would use
%       the default colormap (parula)


if mod(size(data,2),3)>0
    error('matrix must have 3 columns')
end
nSets = size(data,2)/3 ;

if ~isempty(cmap)
    cm = eval( ['colormap( ' cmap '(nSets))']) ;
end

for i = 1:nSets
    if any( strcmp('color',varargin) )
        plotHandle(i) = plot3(data(:,3*i-2),data(:,3*i-1),data(:,3*i),varargin{:}) ;
    elseif exist('cm','var')
        plotHandle(i) = plot3(data(:,3*i-2),data(:,3*i-1),data(:,3*i),varargin{:},'markerfacecolor',cm(i,:)) ;
    else
        plotHandle(i) = plot3(data(:,3*i-2),data(:,3*i-1),data(:,3*i),varargin{:}) ;
    end
    hold on
end

if nargout==1
    varargout = {plotHandle} ;
end