h = 0.025;

Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 2, 0, h};
Point(4) = {0, 2, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1,2,3,4};

Plane Surface(1) = {1};

Physical Line(1) = {1,2,3,4};
Physical Surface(1) = {1};
