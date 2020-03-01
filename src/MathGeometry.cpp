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
real tri_area(real x1, real x2, real x3, real y1, real y2, real y3) {
    real result = 0.5*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
    if (result < 1.e-10) {
      cout << "ERROR: triangle with bad area" << endl;
      cout << "triangle area = " << result << endl;
      std::exit(0);
    }
    return result;
 }

// real sumArray2d(Array2D<real> A, int col) {
//   real sum_ = 0.0;
//   for (size_t i = 0; i < A.nrows; i++) {
//     sum_ += A(i,col);
//   }
//   return sum_;
// }


