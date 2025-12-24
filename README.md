# Triangle-of-support_JEB
MATLAB code used for: Triangle of Support and Joint-reaction forces in the primate feeding
system

This repository contains the most important files in the main folder and several accessory
functions in the 'functions' folder.

The main folder contains the code used to calculate the joint reaction forces (JRFs) and 
the position of the resultant muscle force (RMF) with respect to the triangle of support 
(ToS). There are two main scripts that are used to calculate the JRFs and the position of 
the RMF and one script that needs to be run before to calculate the moment arms of the jaw 
adductor muscles at different gapes.

- **calculateMuscleMoments:** This script calculates the moment arms and relative lengths 
of the muscle segments of the three jaw addutor muscles at different gape angles. It exports
a MATLAB file that is needed to run the following two scripts.
- **TOS_muscleForceEffect:** This script calculates the the effect of differences in relative 
force of the three jaw adductor muscles on JRFs and the position of the RMF with respect to the ToS.
- **ToS_gapeEffect:** This script calculates the effect of variation in muscle force 
capabilities due to gape (as a consequence of the length-tension curve) on JRFs and the 
position of the RMF with respect to the ToS.

Each of these scripts require some data that needs to be calculated beforehand.
- **BONE MODEL data:** MATLAB file that contains the 3D triangulated cranium and mandible models. 
This file is exported from the 'modelPreparation' MATLAB app available from 
https://jiriarte-sewanee.github.io/#software
- **MODEL ROTATION data:** MATLAB file that contains the homogeneous tranformation matrices 
needed to rotate the mandible for each gape angle. 
This file is exported from the 'modelPreparation' MATLAB app available from 
https://jiriarte-sewanee.github.io/#software
- **MUSCLE data:** MATLAB file that contains the position, length, and moment arms data 
for the three jaw adductor muscles at different gapes. This data is obtained by running 
the "calculateMuscleMoments.m" script.
- **LANDMARK data** pp file that contains the cartesian coordinates of the left and right 
TMJs and the teeth of the 3D model. This file is exported directly from Meshlab (https://www.meshlab.net/).

---
**NOTE 1:** The code needs several functions from the geometry processing library MatGeom
(https://github.com/mattools/matGeom)

**NOTE 2:** To run the 'calculateMuscleMoment' script, we need the 'intersectPlaneSurf'
function. You can download the function here:
https://www.mathworks.com/matlabcentral/fileexchange/32256-intersectplanesurf-ii
Make sure to compile the attached cpp file as: 
`mex('IntersectPlaneTriangle.cpp','-v')`

**NOTE 3:** The function 'directionalRayTriangleIntersect' calls the mex function
'rayTriangleIntersection_mlc_mex'. The compiled mex version is substantially faster than
the original code. We have included compiled mex files for Windows (.mexw64), for macOS
with Intel processor (.mexmaci64), and for macOS with Apple silicon processor
(.mexmaca64). We also provide the source code (rayTriangleIntersection_mlc), so that you
can create your own compiled mex file, if needed.
