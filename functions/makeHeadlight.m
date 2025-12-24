function makeHeadlight(varargin)

% adds a front light to all images in the current figure

h = get(gcf,'children') ;
hAxis = findobj(h,'type','Axes') ;

az = 30 ; 
el = 0 ;
for i = 1:numel(hAxis)
    set(gcf,'CurrentAxes',hAxis(i))
    
    lightHandle = findobj(hAxis(i),'type','light') ;
    
    if isempty(lightHandle)
        camlight(az,el)
        clear('lightHandle')
    else
        for ii = 1:numel(lightHandle)
            camlight(lightHandle(ii),az,el)
        end
        clear('lightHandle','ii')
    end
    material dull
    lighting gouraud
end