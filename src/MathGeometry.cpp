#include "../include/MathGeometry.h"

//********************************************************************************
//* Compute the area of the triangle defined by the nodes, 1, 2, 3.
//*
//*              3 (x3,y3)
//*              o 
//*             / \ 
//*            /   \
//* (x1,y1) 1 o-----o 2 (x2,y2)
//*
//* Nodes must be ordered counterclockwise (otherwise it gives negative area)
//*
//********************************************************************************
double tri_area(double x1, double x2, double x3, double y1, double y2, double y3) {
   return 0.5*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
 }

