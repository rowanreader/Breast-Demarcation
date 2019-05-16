# Breast-Demarcation

Given a series of points, fits a circle to them (using code from 
https://www.mathworks.com/matlabcentral/fileexchange/55304-best-fit-3d-circle-to-a-set-of-points)

Rotates the circle be parallel to the XY plane, and adds vectors t degrees apart (specified in vectors.m)
Draws dirNum (specified in vectors.m) points along vectors (perpindicular to circle) separated by space (specified in vectors.m)
Unrotates the circle and vectors back to original position

Output:
graph of points used to fit the circle, the circle, and the points for the vectors
nx3 array of points along vectors
