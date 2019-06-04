//*******************************************************************************
// One-dimensional Euler solver for a Sod's shock tube problem.
//
//        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//        translated to C++ by Luke McCulloch
//
// the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//
// This is Version 1 (2010).
// 
// This F90 program was written and made available for download
// for an educational purpose. Comments are welcome.
//
// This file may be updated in future.
//
// Katate Masatsuka, December 2010. http://www.cfdbooks.com
//*******************************************************************************

//*******************************************************************************
// --- Main program for the 1D Euler shock-tube solver.
//
// This code solves the Sod's shock tube problem which is described
// in Section 7.13.3 of "I do like CFD, VOL.1": Sod's problem 1, Figure 7.12.2.
//
// - t=0                               - t=tf
// Density                             Density
//   ****************|                 *********\
//                   |                           \
//                   |                            \
//                   |                             ****|
//                   |                                 |
//                   |                                 ****|
//                   ***************                       ***********
//
// Methods employed:
//   - Roe's flux
//   - Minmod limiter
//   - Two-stage Runge-Kutta time stepping scheme
//
// Input ---------------------------------------------------
//    ncells = # of cells on a grid.
//        tf = Final time
//       cfl = CFL number (<1)
//      xmin = Left end of the domain
//      xmax = Right end of the domain
// 
// Output --------------------------------------------------
//  "solution.dat" = Data file containing for each cell,
//                   cell-center coordinate, density, velocity, pressure, 
//                   entropy, in the following format:
//
//       write(*,*) ncells
//      do i=1,ncells
//       write(*,*) xc(i), rho(i), u(i), p(i), entropy(i)
//      end do
// 
//     Use the matlab program, oned_euler_v1.m, to plot the solutions.
//
//
//  Note: Explore other shock tube problems by changing the initial condition
//        and the final time (Other problems are described in Section 7.13.3
//        of "I do like CFD, VOL.1").
//
//  Note: Other limiters may be explored (see CFD textbooks).
//
//  Note: Other flux functions may be explored.
//        Various subroutines are available at cfdbooks.com: Osher, Van Leer, etc.
//
//  Note: Boundary condition need to be modified for problems with hard
//        boundaries.
//
//
// Katate Masatsuka, December 2010. http://www.cfdbooks.com
//
//
// 12-29-10: Some compiler warnings fixed.
//
//*******************************************************************************

#include <iostream>     // std::cout, std::fixed
//#include <iomanip>      // std::setprecision - only works for output :(

//======================================
// my simple array class template (type)
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"

//======================================
// my simple vector class template 
#include "../include/vector.h"

//======================================
using namespace std;

//======================================
//fwd declarations
struct cell_data; 
void initialize( cell_data* cell, int ncells, 
                float dx, float xmin, const float gamma );
//float* w2u( float w[3] );
void w2u( float w[3], float u[3]  );


// use an array of structs (very inefficient//)
struct cell_data{
  float xc;  // Cell-center coordinate
  float u[3];   // Conservative variables = [rho, rho*u, rho*E]
  float u0[3];  // Conservative variables at the previous time step
  float w[3];   // Primitive variables = [rho, u, p]
  float dw[3];  // Slope (difference) of primitive variables
  float res[3]; // Residual = f_{j+1/2) - f_{j-1/2)
};

struct constants{
    const float  zero = 0.0;
    const float   one = 1.0;
    const float  half = 0.5;
    const float gamma = 1.4;  //Ratio of specific heats for air
};

void Euler1D(){
    //Numeric parameters: [Note: no Fortran-like way to handle precision?]
    //const int p2 = 10;
    const float  zero = 0.0;
    const float   one = 1.0;
    const float  half = 0.5;
    const float gamma = 1.4;  //Ratio of specific heats for air

    float xmin, xmax; //Left and right ends of the domain
    float dx;         //Cell spacing (uniform grid)
    float t, tf;      //Current time and final time
    float cfl, dt;    //CFL number and global time step
    int   ncells;     //Total number of cells
    int   nsteps;     //Number of time steps
    int   itime;      //Index for time stepping
    int   istage;     //Index for Runge-Kutta stages
    int   i, j;

    //Local variables used for computing numerical fluxes.
    Array2D<float>  dwl(3,1), dwr(3,1); //Slopes between j and j-1, j and j+1
    Array2D<float>  wL(3,1), wR(3,1);   //Extrapolated states at a face
    Array2D<float>  flux(3,1);          //Numerical flux

//--------------------------------------------------------------------------------
// 0. Input parameters and initial condition.

//Parameters
  ncells =  80;   // Number of cells
      tf = 1.7;   // Final time
     cfl = 0.8;   // CFL number
    xmin =-5.0;   // Left boundary coordinate
    xmax = 5.0;   // Right boundary coordinate

    
// Allocate the cell array: 2 ghost cells, 0 on the left and ncells+1 on the right.
//  E.g., data in cell j is accessed by cell(j)%xc, cell(j)%u, cell(j)%w, etc.
    //Array2D<cell_data> cell(ncells+1,1);//experimental
    struct cell_data cell[ncells+1];      //Array of cell-data

// Cell spacing (grid is uniform)
  dx = (xmax-xmin)/float(ncells);

// The initial condition for Sod's shock tube problem (I Do Like CFD, VOL.1, page 199).
// [Note: Change these data (and tf) to solve different problems.]
    initialize(cell, ncells, dx, xmin, gamma);

//--------------------------------------------------------------------------------
// Time stepping loop to reach t = tf 
//--------------------------------------------------------------------------------

    t = zero; //Initialize the current time.
}


void initialize( cell_data* cell, int ncells, 
                float dx, float xmin, const float gamma){
    /* 
        The initial condition for Sod's shock tube problem
    */
    for ( int i = 0; i < ncells+2; ++i ) {
        if (i <= ncells/2) {
            cell[i].w[0] = 1.0;   //Density  on the left
            cell[i].w[1] = 0.0;   //Velocity on the left
            cell[i].w[2] = 1.0;   //Pressure on the left
        } else {
            cell[i].w[0] = 0.125; //Density  on the right
            cell[i].w[1] = 0.0;   //Velocity on the right
            cell[i].w[2] = 0.1;   //Pressure on the right
        }

        //cell[i].u  = w2u( cell[i].w );      //Compute the conservative variables
        w2u( cell[i].w, cell[i].u );
        cell[i].xc = xmin+float(i-1)*dx;    //Cell center coordinate
    }
    return;
}


//********************************************************************************
//* Compute U from W
//*
//* ------------------------------------------------------------------------------
//*  Input:  w = primitive variables (rho, u, p, 0, 0)
//* Output:  u = conservative variables (rho, rho*u, rho*E, 0, 0)
//* ------------------------------------------------------------------------------
//* 
//********************************************************************************
// float* w2u( float w[3] ) {

//     float gamma = 1.4;
//     float half = 0.5;
//     float one  = 1.0;
//     float u[3];

//     u[0] = w[0];
//     u[1] = w[0]*w[1];
//     u[2] = w[2]/(gamma-one)+half*w[0]*w[1]*w[1];
//     return *u;
// }
//--------------------------------------------------------------------------------
//********************************************************************************
void w2u( float w[3], float u[3] ) {

    float gamma = 1.4;
    float half = 0.5;
    float one  = 1.0;

    u[0] = w[0];
    u[1] = w[0]*w[1];
    u[2] = w[2]/(gamma-one)+half*w[0]*w[1]*w[1];
    return;
}
//--------------------------------------------------------------------------------


int main(){
    Euler1D();
    return 0;
}