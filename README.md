# Triangle-of-support_JEB
Code use for: Triangle of Support and Joint-reaction forces in the primate feeding system

NOTE 1: The code needs several functions from the geometry processing library MatGeom
(https://github.com/mattools/matGeom) 
NOTE 2: To run the 'calculateMuscleMoment' script, we need the 'intersectPlaneSurf' function. You can download the function here:
https://www.mathworks.com/matlabcentral/fileexchange/32256-intersectplanesurf-ii
Make sure to compile the attached cpp file as: 
	mex('IntersectPlaneTriangle.cpp','-v')
NOTE 3: The function 'directionalRayTriangleIntersect' calls the mex function 'rayTriangleIntersection_mlc_mex'. The compiled mex version is substantially faster than the original code. We have included compiled mex files for Windows (.mexw64), for macOS with Intel processor (.mexmaci64), and for macOS with Apple silicon processor (.mexmaca64). We also provide the source code (rayTriangleIntersection_mlc), so that you can create your own compiled mex file, if needed.
