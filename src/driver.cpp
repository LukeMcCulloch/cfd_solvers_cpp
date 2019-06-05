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
//  [] ***************|                [] ********\
//                   |                           \
//                   |                            \
//                   |                            [] ***|
//                   |                                 |
//                   |                                [] ***|
//                  [] **************                      [] **********
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

//=================================
#include <iostream>     // std::cout, std::fixed
//#include <iomanip>    // std::setprecision - only works for output :(
#include <math.h>       // sqrt 
//=================================
#include <cstring>
#include <string.h>

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
class Solver;
// void initialize( cell_data* cell, int ncells, 
//                 float dx, float xmin, const float gamma );
// void w2u( float w[3], float u[3]  );



// use an array of structs (may be inefficient//)
// struct cell_data{
//     float xc;  // Cell-center coordinate
//     float u[3];   // Conservative variables = [rho, rho*u, rho*E]
//     float u0[3];  // Conservative variables at the previous time step
//     float w[3];   // Primitive variables = [rho, u, p]
//     float dw[3];  // Slope (difference) of primitive variables
//     float res[3]; // Residual = f_{j+1/2) - f_{j-1/2)
// };


// use an array of structs (may be inefficient//)
struct cell_data{
    float xc;  // Cell-center coordinate
    Array2D<float>* u;   // Conservative variables = [rho, rho*u, rho*E]
    Array2D<float>* u0;  // Conservative variables at the previous time step
    Array2D<float>* w;   // Primitive variables = [rho, u, p]
    Array2D<float>* dw;  // Slope (difference) of primitive variables
    Array2D<float>* res; // Residual = f_{j+1/2) - f_{j-1/2)

    cell_data(){
        u = new Array2D<float>(3,1);
        u0 = new Array2D<float>(3,1);
        w = new Array2D<float>(3,1);
        dw = new Array2D<float>(3,1);
        res = new Array2D<float>(3,1);
    }
};



class Solver{

public:

    //constructor
    Solver();
    // destructor
    ~Solver();


    void Euler1D();
    void initialize( int ncells, 
                float dx, float xmin, const float gamma);
    float timestep(float cfl, float dx, float gamma, int ncells);
    //void w2u( float w[3], float u[3] );
    void w2u( Array2D<float>*& w, Array2D<float>*& u );
    float minmod(float a, float b);

    struct constants{
        const float  zero = 0.0;
        const float   one = 1.0;
        const float  half = 0.5;
        const float gamma = 1.4;  //Ratio of specific heats for air
    };

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
    Array2D<float>  dwl = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  dwr = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  wL = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  wR = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  flux = Array2D<float>(3,1); //Numerical flux
    
    struct cell_data* cell;

};

Solver::Solver(){

//--------------------------------------------------------------------------------
// 0. Input parameters and initial condition.

    printf("\n Custom Parameters\n");
//custom Parameters
    ncells =  80;   // Number of cells
      tf = 1.7;   // Final time
     cfl = 0.8;   // CFL number
    xmin =-5.0;   // Left boundary coordinate
    xmax = 5.0;   // Right boundary coordinate

    
    printf("\n Cell Struct array\n");
// Allocate the cell array: 2 ghost cells, 0 on the left and ncells+1 on the right.
//  E.g., data in cell j is accessed by cell(j)%xc, cell(j)%u, cell(j)%w, etc.
    //Array2D<cell_data> cell(ncells+1,1);//experimental
    cell = new cell_data[ncells+1];      //Array of cell-data

// Cell spacing (grid is uniform)
    dx = (xmax-xmin)/float(ncells);


    printf("\n Initialize solver\n");
// The initial condition for Sod's shock tube problem (I Do Like CFD, VOL.1, page 199).
// [Note: Change these data (and tf) to solve different problems.]
    initialize(ncells, dx, xmin, gamma);

    
    printf("\n local arrays\n");
    //Local variables used for computing numerical fluxes.
    // dwl = new Array2D<float>(3,1); 
    // dwr = new Array2D<float>(3,1); 
    // wL = new Array2D<float>(3,1); 
    // wR = new Array2D<float>(3,1); 
    // flux = new Array2D<float>(3,1); 
}



Solver::~Solver(){
    delete[] cell;
    // delete[] dwl;
    // delete[] dwr;
    // delete[] wL;
    // delete[] wR;
    // delete[] flux;
}

void Solver::Euler1D(){
    
//--------------------------------------------------------------------------------
// Time stepping loop to reach t = tf 
//--------------------------------------------------------------------------------
    printf("\nEuler1D\n");
    t = zero;      //Initialize the current time.
    nsteps = 0;    //Initialize the number of time steps.
    //50000 is large enough to reach tf=1.7.
    for ( int itime = 0; itime < 500; ++i ) {
        if (t==tf) break;                   //Finish if the final time is reached.
        //dt = 1.;
        dt = timestep(cfl,dx,gamma,ncells); //Compute the global time step.
        if (t+dt > tf) dt =  tf - t;        //Adjust dt to finish exactly at t=tf.
        t = t + dt;                 //Update the current time.
        nsteps = nsteps + 1;                //Count the number of time steps.

//---------------------------------------------------
// Runge-Kutta Stages
//
// Two-stage Runge-Kutta scheme:
//  1. u^*     = u^n - dt/dx*Res(u^n)
//  2. u^{n+1} = 1/2*u^n + 1/2*[u^*- dt/dx*Res(u^*)]
//---------------------------------------------------
        for ( int istage = 0; istage < 2; ++istage ) {

            //(1) Residual computation: compute cell(:)%res(1:3).

            // Compute the slopes (as difference) at every cell.
            // NB: for uniform meshes, difference (du) can be used in place of gradient (du/dx).
            for ( int j = 1; j < ncells+1; ++j ) {

                //dwl = cell[j].w   - cell[j-1].w // diff( w_{j}   , w_{j-1} )
                //dwr = cell[j+1].w -   cell[j].w // diff( w_{j+1} , w_{j}   )

                //    Apply a slope limiter.
                //    (minmod: zero if opposite sign, otherwise the one of smaller magnitude.)
                for (int i = 0; i < 3; ++i){
                //    cell[j].dw[i] = minimod( dwl[i], dwr[i] );
                }
            }//end reconstruction, j
        } //end rk stages, istage
    } //end time stepping
}


void Solver::initialize( int ncells, float dx, float xmin, const float gamma){
    //
    //The initial condition for Sod's shock tube problem
    // for ( int i = 0; i < ncells+2; ++i ) {
    //     if (i <= ncells/2) {
    //         cell[i].w[0] = 1.0;   //Density  on the left
    //         cell[i].w[1] = 0.0;   //Velocity on the left
    //         cell[i].w[2] = 1.0;   //Pressure on the left
    //     } else {
    //         cell[i].w[0] = 0.125; //Density  on the right
    //         cell[i].w[1] = 0.0;   //Velocity on the right
    //         cell[i].w[2] = 0.1;   //Pressure on the right
    //     }

    //     //w2u( cell[i].w, cell[i].u );        //Compute the conservative variables
    //     //cell[i].xc = xmin+float(i-1)*dx;    //Cell center coordinate
    // }
    return;
}
//******************************************************************************
// Compute the global time step dt restricted by a given CFL number.
//
// ------------------------------------------------------------------------------
//  Input: CFL number
// Output: Global time step, dt.
// ------------------------------------------------------------------------------
//
//******************************************************************************
float Solver::timestep(float cfl, float dx, float gamma, int ncells){
    float dt;             //Output
    //Local variables
    float one = 1.0;
    float u, c, max_speed;
    int   i;

    max_speed = -one;

    for ( int i = 1; i < ncells+1; ++i ) {
        u = cell[i].w->array[1,0];                          //Velocity
        c = sqrt(gamma*cell[i].w->array[2,0]/cell[i].w->array[0,0]); //Speed of sound
        max_speed = max( max_speed, abs(u)+c );

        dt = cfl*dx/max_speed; //CFL condition: dt = CFL*dx/max_wavespeed, CFL <= 1.
    }
    return dt;
}

//***************************************************************************
// Minmod limiter
// --------------------------------------------------------------------------
//  Input: two real values, a and b
// Output: minmod of a and b.
// --------------------------------------------------------------------------
// 
//***************************************************************************
 float Solver::minmod(float a, float b){

    float minmod;

    //Local parameter
    float zero = 0.0;

    if (a*b <= zero) {
        minmod = zero;               // a>0 and b<0; or a<0 and b>0
    }else if (abs(a)<abs(b)) {
        minmod = a;                  // |a| < |b|
    }else if (abs(b)<abs(a)) {
        minmod = b;                  // |a| > |b|
    }else{
        minmod = a;                  // Here, a=b, so just take a or b.
    }
    return minmod;
}
//--------------------------------------------------------------------------------

//********************************************************************************
//* Compute U from W
//*
//* ------------------------------------------------------------------------------
//*  Input:  w = primitive variables (rho, u, p, 0, 0)
//* Output:  u = conservative variables (rho, rho*u, rho*E, 0, 0)
//* ------------------------------------------------------------------------------
//* 
//********************************************************************************
void Solver::w2u( Array2D<float>*& w, Array2D<float>*& u ) {

    float gamma = 1.4;
    float half = 0.5;
    float one  = 1.0;

    u[0] = w[0];
    u[1] = w[0]*w[1];
    u[2] = (w->array[2,0]/(gamma-one))+half*w->array[0]*w->array[1]*w->array[1];
    return;
}
//--------------------------------------------------------------------------------


int main(){
    Solver solver;
    //solver.Euler1D();
    return 0;
}