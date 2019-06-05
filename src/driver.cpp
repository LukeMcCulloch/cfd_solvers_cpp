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
// struct cell_data{
//     float xc;  // Cell-center coordinate
//     Array2D<float>* u;   // Conservative variables = [rho, rho*u, rho*E]
//     Array2D<float>* u0;  // Conservative variables at the previous time step
//     Array2D<float>* w;   // Primitive variables = [rho, u, p]
//     Array2D<float>* dw;  // Slope (difference) of primitive variables
//     Array2D<float>* res; // Residual = f_{j+1/2) - f_{j-1/2)

//     cell_data(){
//         u = new Array2D<float>(3,1);
//         u0 = new Array2D<float>(3,1);
//         w = new Array2D<float>(3,1);
//         dw = new Array2D<float>(3,1);
//         res = new Array2D<float>(3,1);
//     }
// };


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


    void Euler1D();
    void initialize( int ncells, 
                float dx, float xmin, const float gamma);
    float timestep(float cfl, float dx, float gamma, int ncells);
    //void w2u( float w[3], float u[3] );
    void w2u( Array2D<float>& w, Array2D<float>& u );
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
    // init arrays here
    // https://stackoverflow.com/questions/11490988/c-compile-time-error-expected-identifier-before-numeric-constant
    Array2D<float>  dwl = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  dwr = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  wL = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  wR = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  flux = Array2D<float>(3,1); //Numerical flux

    cell_data* cell;

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
    //Array2D<cell_data> cell(ncells+1,1); //experimental
    cell = new cell_data[ncells+2];      //Array of cell-data

// Cell spacing (grid is uniform)
    dx = (xmax-xmin)/float(ncells);


    printf("\n Initialize solver\n");
// The initial condition for Sod's shock tube problem (I Do Like CFD, VOL.1, page 199).
// [Note: Change these data (and tf) to solve different problems.]
    initialize(ncells, dx, xmin, gamma);

    
    printf("\n local arrays\n");
    //Local variables used for computing numerical fluxes.
    // done already via 
    // https://stackoverflow.com/questions/11490988/
    //    c-compile-time-error-expected-identifier-before-numeric-constant
    // 
    // dwl = new Array2D<float>(3,1); 
    // dwr = new Array2D<float>(3,1); 
    // wL = new Array2D<float>(3,1); 
    // wR = new Array2D<float>(3,1); 
    // flux = new Array2D<float>(3,1); 
}



Solver::~Solver(){
    delete[] cell;
    // delete[] dwl; //raw arrays only
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
    for ( int itime = 0; itime < 50000; ++i ) {
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

                dwl = cell[j].w   - cell[j-1].w; // diff( w_{j}   , w_{j-1} )
                dwr = cell[j+1].w -   cell[j].w; // diff( w_{j+1} , w_{j}   )

                //    Apply a slope limiter.
                //    (minmod: zero if opposite sign, otherwise the one of smaller magnitude.)
                for (int i = 0; i < 3; ++i){
                    cell[j].dw(i) = minmod( dwl(i), dwr(i) );
                }
            }//end reconst

            // init the residuals:
            for (int j = 1; j < ncells+1; ++j){
                cell[j].res = zero;
            }

        // Compute the residuals: residual_j = flux_{j+1/2} - flux_{j-1/2}.
        // Here, compute the flux at j+1/2 and add it to the left cell and subtract
        // from the right cell. Only the internal faces are considered; the left
        // and right most faces are considered later.

        //     j+1/2
        //   | wL|   |
        //   |  /|wR |
        //   | / |\  |
        //   |/  | \ |
        //   |   |  \|
        //   |   |   |
        //     j  j+1
        //

        // flux comparison

        for (int j = 1; j < ncells; ++j){

            wL = cell[j  ].w + half*cell[j  ].dw //State extrapolated to j+1/2 from j
            wR = cell[j+1].w - half*cell[j+1].dw //State extrapolated to j+1/2 from j+1
            flux = roe_flux(wL,wR,gamma)           //Numerical flux at j+1/2
            cell[j  ].res = cell[j  ]res + flux   //Add it to the left cell.
            cell[j+1].res = cell[j+1]res - flux   //Subtract from the right cell.

        }

        } //end rk stages, istage
//---------------------------------------------------
// End of Runge-Kutta Stages
//---------------------------------------------------
    } //end time stepping
//--------------------------------------------------------------------------------
// End of time stepping
//--------------------------------------------------------------------------------

}
//********************************************************************************
// End of program
//********************************************************************************



void Solver::initialize( int ncells, float dx, float xmin, const float gamma){
    //
    //The initial condition for Sod's shock tube problem
    for ( int i = 0; i < ncells+2; ++i ) {
        if (i <= ncells/2) {
            cell[i].w(0) = 1.0;   //Density  on the left
            cell[i].w(1) = 0.0;   //Velocity on the left
            cell[i].w(2) = 1.0;   //Pressure on the left
        } else {
            cell[i].w(0) = 0.125; //Density  on the right
            cell[i].w(1) = 0.0;   //Velocity on the right
            cell[i].w(2) = 0.1;   //Pressure on the right
        }

        w2u( cell[i].w, cell[i].u );        //Compute the conservative variables
        cell[i].xc = xmin+float(i-1)*dx;    //Cell center coordinate
    }
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
        u = cell[i].w(1);                          //Velocity
        c = sqrt(gamma*cell[i].w(2)/cell[i].w(0)); //Speed of sound
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
void Solver::w2u( Array2D<float>& w, Array2D<float>& u ) {

    float gamma = 1.4;
    float half = 0.5;
    float one  = 1.0;

    u(0) = w(0);
    u(1) = w(0)*w(1);
    u(2) = ( w(2)/(gamma-one) ) + half*w(0)*w(1)*w(1);
    return;
}
//--------------------------------------------------------------------------------

//*******************************************************************************
// Compute U from W
//
// ------------------------------------------------------------------------------
//  Input:  u = conservative variables (rho, rho*u, rho*E, 0, 0)
// Output:  w = primitive variables (rho, u, p, 0, 0)
// ------------------------------------------------------------------------------
// 
//*******************************************************************************
 function u2w(u,gamma) result(w)
 implicit none
 float u(3), gamma !Input
 float             :: w(3)        !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = (gamma-one)*( u(3) - half*w(1)*w(2)*w(2) )

 end function u2w
//--------------------------------------------------------------------------------



//*******************************************************************************
// -- Roe's Flux Function without entropy fix---
//
// P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
// Schemes, Journal of Computational Physics, 43, pp. 357-372.
//
// ------------------------------------------------------------------------------
//  Input:   wL(1:3) =  left state (rhoL, uL, pL)
//           wR(1:3) = right state (rhoR, uR, pR)
//
// Output:  flux(1:3) = numerical flux for the Euler equations (the Roe flux)
// ------------------------------------------------------------------------------
// 
// Katate Masatsuka, December 2010. http://www.cfdbooks.com
//*******************************************************************************
 float* roe_flux(wL,wR,gamma) result(flux)

 implicit none
 float wL(3), wR(3), gamma //  Input (conservative variables rho*[1, v, E])
 float flux(3)             // Output (numerical flux across L and R states)

//Local parameters
 float    zero = 0.0
 float     one = 1.0
 float    four = 4.0
 float    half = 0.5
 float quarter = 0.25
//Local variables
 float uL(3), uR(3)
 float rhoL, rhoR, vL, vR, pL, pR   // Primitive variables.
 float aL, aR, HL, HR               // Speeds of sound.
 float RT,rho,v,H,a                 // Roe-averages
 float drho,du,dP,dV(3)
 float ws(3),Da, R(3,3)
 integer :: j, k

    uL = w2u(wL,gamma)
    uR = w2u(wR,gamma)

//Primitive and other variables.
//  Left state
    rhoL = wL(1)
      vL = wL(2)
      pL = wL(3)
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
//  Right state
    rhoR = wR(1)
      vR = wR(2)
      pR = wR(3)
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

//First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*v*v) )

//Differences in primitive variables.
   drho = rhoR - rhoL
     du =   vR - vL
     dP =   pR - pL

//Wave strength (Characteristic Variables).
   dV(1) =  half*(dP-rho*a*du)/(a*a)
   dV(2) = -( dP/(a*a) - drho )
   dV(3) =  half*(dP+rho*a*du)/(a*a)

//Absolute values of the wave speeds (Eigenvalues)
   ws(1) = abs(v-a)
   ws(2) = abs(v  )
   ws(3) = abs(v+a)

//Modified wave speeds for nonlinear fields (the so-called entropy fix, which
//is often implemented to remove non-physical expansion shocks).
//There are various ways to implement the entropy fix. This is just one
//example. Try turn this off. The solution may be more accurate.
   Da = max(zero, four*((vR-aR)-(vL-aL)) )
   if (ws(1) < half*Da) { ws(1) = ws(1)*ws(1)/Da + quarter*Da }
   Da = max(zero, four*((vR+aR)-(vL+aL)) )
   if (ws(3) < half*Da) { ws(3) = ws(3)*ws(3)/Da + quarter*Da }

//Right eigenvectors
   R(1,1) = one
   R(2,1) = v - a
   R(3,1) = H - v*a

   R(1,2) = one
   R(2,2) = v
   R(3,2) = half*v*v

   R(1,3) = one
   R(2,3) = v + a
   R(3,3) = H + v*a

//Compute the average flux.
   flux = half*( euler_physical_flux(wL) + euler_physical_flux(wR) )

//Add the matrix dissipation term to complete the Roe flux.
  do j = 1, 3
   do k = 1, 3
    flux(j) = flux(j) - half*ws(k)*dV(k)*R(j,k) 
   end do
  end do

} //end roe_flux
//--------------------------------------------------------------------------------


int main(){
    Solver solver;
    solver.Euler1D();
    return 0;
}