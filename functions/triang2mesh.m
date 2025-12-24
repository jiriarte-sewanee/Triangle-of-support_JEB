function FV = triang2mesh(TR) 

% function transforms a triangulated mesh into a structure mesh

FV.faces    = TR.ConnectivityList ;
FV.vertices = TR.Points ;