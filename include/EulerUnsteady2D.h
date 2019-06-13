//=================================
// include guard
#ifndef __eulerUnsteady2d_INCLUDED__
#define __eulerUnsteady2d_INCLUDED__


//======================================
// my simple array class template (type)
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"



namespace EulerSolver2D
{


// use an array of structs (may be inefficient//)
struct cell_data{
    float xc;  // Cell-center coordinate
    Array2D<float> u  = Array2D<float>(3,1);  // Conservative variables = [rho, rho*u, rho*E]
    Array2D<float> u0 = Array2D<float>(3,1);  // Conservative variables at the previous time step
    Array2D<float> w  = Array2D<float>(3,1);  // Primitive variables = [rho, u, p]
    Array2D<float> dw = Array2D<float>(3,1);  // Slope (difference) of primitive variables
    Array2D<float> res= Array2D<float>(3,1);  // Residual = f_{j+1/2) - f_{j-1/2)
};



class Solver{

public:

    //constructor
    Solver();
    // destructor
    ~Solver();
};



//=================================
// the driver function
void driverEuler2D();


}//end namespace

#endif 