# Triangle-of-support_JEB
MATLAB code used for: Triangle of Support and Joint-reaction forces in the primate feeding
system

This repository contains the most important files in the main folder and several accessory
functions in the 'functions' folder. We also provide an 'Examples' folder with data from 
two sample species to test the code.

The main folder contains the code used to calculate the joint reaction forces (JRFs) and 
the position of the resultant muscle force (RMF) with respect to the triangle of support 
(ToS). Two main scripts are used to calculate the JRFs and the position of the RMF:

- **TOS_muscleForceEffect:** This script calculates the effect of differences in relative 
force of the three jaw adductor muscles on JRFs and the position of the RMF with respect to the ToS.
- **ToS_gapeEffect:** This script calculates the effect of variation in muscle force 
capabilities due to gape (as a consequence of the length-tension curve) on JRFs and the 
position of the RMF with respect to the ToS.

%--------------

**IMPORTANT:** Each of these scripts requires some data that needs to be calculated beforehand. 
Also, the code assumes that all required files are in the current MATLAB folder.

- **BONE MODEL data:** MATLAB file that contains the 3D triangulated cranium and mandible models. 
This file is exported from the 'modelPreparation' MATLAB app available from 
https://jiriarte-sewanee.github.io/#software
- **MODEL ROTATION data:** MATLAB file that contains the homogeneous transformation matrices 
needed to rotate the mandible for each gape angle. 
This file is exported from the 'modelPreparation' MATLAB app available from 
https://jiriarte-sewanee.github.io/#software
- **MUSCLE data:** MATLAB file that contains the position, length, and moment arms data 
for the three jaw adductor muscles at different gapes. This data is obtained by running 
the "calculateMuscleMoments.m" script, which is available in the main folder.
- **LANDMARK data** pp file that contains the Cartesian coordinates of the left and right 
TMJs and the teeth of the 3D model. This file is exported directly from Meshlab (https://www.meshlab.net/).

---

**NOTE 1:** The 'ToS_gapeEffect' script requires to load the 'normalized_LTcurve.mat' file
provided in the main folder. You'll have to change the location of the folder on line 94 
of the 'ToS_gapeEffect' script to the specific folder on your computer.


**NOTE 2:** To run the 'calculateMuscleMoment' script, we need the 'intersectPlaneSurf'
function. You can download the function here:
https://www.mathworks.com/matlabcentral/fileexchange/32256-intersectplanesurf-ii
Make sure to compile the attached cpp file as: 
`mex('IntersectPlaneTriangle.cpp','-v')`

**NOTE 3:** The function 'directionalRayTriangleIntersect' calls the mex function
'rayTriangleIntersection_mlc_mex'. The compiled mex version is substantially faster than
the original code. We provide compiled mex files for Windows (.mexw64), for macOS
with Intel processor (.mexmaci64), and for macOS with Apple silicon processor
(.mexmaca64). We also provide the source code (rayTriangleIntersection_mlc), so that you
can create your own compiled mex file, if needed.
