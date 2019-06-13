//********************************************************************************
//* Educationally-Designed Unstructured 2D (EDU2D) Code
//*
//*
//*                      --- This is EDU2D-Euler-RK2 ---
//*
//*
//* EDU2D-Euler-RK2: An Euler code with
//*
//*    - Node-centered finite-volume discretization
//*    - 2-stage Runge-Kutta explicit time-stepping scheme (RK2)
//*
//*
//*
//*             specially set up for a shock-diffraction problem
//*
//*                                 Wall
//*                         --------------------
//*     Post-shock (inflow) |                  |
//*                         |->Shock           |            o: Corner node
//*                         |  Mach=5.09       |
//*                  .......o                  |Outflow
//*                    Wall |                  |
//*                         |                  |
//*                         |                  |
//*                         --------------------
//*                               Outflow
//*
//* - Node-centered finite-volume method for unstructured grids (quad/tri/mixed)
//* - Roe flux with an entropy fix and Rotated-RHLL flux
//* - Gradient reconstruction by unweighted least-squares method
//* - Van Albada slope limiter to the primitive variable gradients
//* - 2-Stage Runge-Kutta global time-stepping towards the final time
//* - All quantities are nondimensionalized; velocity and pressure are
//*   nondimensionalized based on the free stream speed of sound
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