function newModel = rotatePatch(model,RotMatrix)

% Function that applies a rotation to a 3D triangulated structure

if isobject(model)
    if all(size(RotMatrix)==4) || size(RotMatrix,2)==16
        newPoints = rotateT(model.Points,RotMatrix) ;
    else
        newPoints = (RotMatrix*model.Points')' ;
    end
    newModel = triangulation(model.ConnectivityList,newPoints) ;
else    
    newModel = model ;
    
    if all(size(RotMatrix)==4) || size(RotMatrix,2)==16
        newModel.vertices = rotateT(model.vertices,RotMatrix) ;
    else
        newModel.vertices = (RotMatrix*model.vertices')' ;
    end
end