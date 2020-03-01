
//=================================
// include guard
#ifndef __MATHGEOMETRY_INCLUDED__
#define __MATHGEOMETRY_INCLUDED__

//#define REAL_IS_DOUBLE true
#ifdef REAL_IS_DOUBLE
  typedef double real;
#else
  typedef float real;
#endif

#include <cmath>

#include "../include/array_template.hpp"
#include "../include/EulerUnsteady2D_basic_package.h"

real tri_area(real x1, real x2, real x3, real y1, real y2, real y3);

//real sumArray2d(Array2D<real> A, int col);


#endif