// Gmsh project created on Fri Oct 19 16:42:12 2012
Point(1) = {0.0, 0.0, 0, 100000.0};
Point(2) = {750000.0, 0, 0, 100000.0};
Circle(1) = {2, 1, 2};
Line Loop(2) = {1};
Plane Surface(3) = {2};
