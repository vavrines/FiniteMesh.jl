//+
Point(1) = {-1, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {0, -6, 0, 1.0};
//+
Point(5) = {0, 6, 0, 1.0};

//+
Circle(1) = {1, 3, 2};
//+
Circle(2) = {2, 3, 1};
//+
Circle(3) = {5, 3, 4};
//+
Circle(4) = {4, 3, 5};
//+
Curve Loop(1) = {3, 4};
//+
Curve Loop(2) = {2, 1};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("inflow") = {3};
//+
Physical Curve("outflow") = {4};
//+
Physical Curve("cylinder") = {2, 1};
//+
Physical Surface("surface") = {1};
