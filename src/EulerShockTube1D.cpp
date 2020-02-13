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
//  ******************|                ***********\
//                    |                            \
//                    |                             \
//                    |                              *****|
//                    |                                   |
//                    |                                   *****|
//                   *****************                         *************
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
//     Use the python program, oned_euler_v1.py, to plot the solutions.
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
//
//*******************************************************************************
//#define CHECKPT {printf("Checkpoint: .s, line .d\n",__FILE__,__LINE__);\
//fflush(stdout);}
#ifdef DEBUG_BUILD
#  define DEBUG(x) fprintf(stderr, x)
#else
#  define DEBUG(x) do {} while (0)
#endif

//=================================
#include <iostream>     // std::cout, std::fixed
#include <fstream>      // write to file
#include <iomanip>    // std::setprecision - only works for output :(
#include <math.h>       // sqrt 
//=================================
#include <cstring> //needed for memset
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
// 1D Euler approximate Riemann sovler
#include "../include/EulerShockTube1D.h"

//======================================
//using namespace std;
//using std::cout;
//using namespace EulerSolver1D;

//======================================
//fwd declarations moved to header



EulerSolver1D::Solver::Solver(){

//--------------------------------------------------------------------------------
// 0. Input parameters and initial condition.

    printf("\n Custom Parameters\n");
    //DEBUG("\n Custom Parameters\n");
//custom Parameters
    ncells =  80;   // Number of cells
      tf = 1.7;   // Final time
     cfl = 0.8;   // CFL number
    xmin =-5.0;   // Left boundary coordinate
    xmax = 5.0;   // Right boundary coordinate

    
    printf("\n Cell Struct array\n");
// Allocate the cell array: 2 ghost cells, 0 on the left and ncells+1 on the right.
//  E.g., data in cell j is accessed by cell[j].xc, cell[j].u, cell[j].w, etc.
    //Array2D<cell_data> cell(ncells+1,1); //experimental
    cell = new cell_data[ncells+2];      //Array of cell-data

// Cell spacing (grid is uniform)
    dx = (xmax-xmin)/float(ncells);


    printf("\n Initialize solver\n");
// The initial condition for Sod's shock tube problem (I Do Like CFD, VOL.1, page 199).
// [Note: Change these data (and tf) to solve different problems.]
    initialize(ncells, dx, xmin, gamma);

}



EulerSolver1D::Solver::~Solver(){
    delete[] cell;
}

void EulerSolver1D::Solver::Euler1D(){
    
//output();
//--------------------------------------------------------------------------------
// Time stepping loop to reach t = tf 
//--------------------------------------------------------------------------------
    printf("\nEuler1D\n");
    t = zero;      //Initialize the current time.
    nsteps = 0;    //Initialize the number of time steps.
    //50000 is large enough to reach tf=1.7.
    for ( int itime = 0; itime < 50000; ++itime ) {
    //for ( int itime = 0; itime < 1; ++itime ) {
        if (t==tf) { 
            //printf("\n breaking \n");
            printf("\n");
            break;
        }                //Finish if the final time is reached.
        dt = timestep(cfl,dx,gamma,ncells); //Compute the global time step.
        printf("\n%f, %f",t,dt);
        if (t+dt > tf){ 
            dt =  tf - t;  
        }    //Adjust dt to finish exactly at t=tf.
        t = t + dt;                         //Update the current time.
        nsteps = nsteps + 1;                //Count the number of time steps.

        //printf("\nRK step \n");
        //---------------------------------------------------
        // Runge-Kutta Stages
        //
        // Two-stage Runge-Kutta scheme:
        //  1. u^*     = u^n - dt/dx*Res(u^n)
        //  2. u^{n+1} = 1/2*u^n + 1/2*[u^*- dt/dx*Res(u^*)]
        //---------------------------------------------------
        for ( int istage = 0; istage < 2; ++istage ) {

            //(1) Residual computation: compute cell(:).res(1:3).

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
                wL = cell[j  ].w + half*cell[j  ].dw; //State extrapolated to j+1/2 from j
                wR = cell[j+1].w - half*cell[j+1].dw; //State extrapolated to j+1/2 from j+1
                flux = roe_flux(wL,wR);           //Numerical flux at j+1/2
                cell[j  ].res = cell[j  ].res + flux;   //Add it to the left cell.
                cell[j+1].res = cell[j+1].res - flux;   //Subtract from the right cell.
            }
            
            // Add boundary fluxes: left end and right end.
            // For the problem considered here, it suffices to simply copy the state
            //  from inside the domain to the ghost cell (no gradient condition).

            //  Left most face: left face of cell i=1.
            wR = cell[1].w - half*cell[1].dw;  //State extrapolated to j-1/2 from j=1
            wL = wR;                           //The same state
            flux = roe_flux(wL,wR);      //Use Roe flux to compute the flux.
            cell[1].res = cell[1].res - flux;  //Subtract the flux: -flux_{j-1/2}.


            //  Right most face: right face of cell i=ncells.
            wL = cell[ncells].w + half*cell[ncells].dw; //State extrapolated to ncells+1/2 from j=ncells
            wR = wL;                                    //The same state
            flux = roe_flux(wL,wR);               //Use Roe flux to compute the flux.
            cell[ncells].res = cell[ncells].res + flux; //Add the flux: +flux_{j+1/2}.

            //(2) Solution update

            if (istage==0){ 
                //  1st Stage of Runge-Kutta: save u^n as u0(:); u^* is stored at u(:).
                //stage01_update : do j = 1, ncells
                for (int j = 1; j < ncells+1; ++j){
                    cell[j].u0 = cell[j].u;            //Save the solution at n for 2nd stage.
                    cell[j].u  = cell[j].u - (dt/dx)*cell[j].res;
                    cell[j].w  = u2w(cell[j].u); //Update primitive variables
                }//end do stage01_update
            }else{
                //  2nd Stage of Runge-Kutta:
                //stage02_update : do j = 1, ncells
                for (int j = 1; j < ncells+1; ++j){
                    cell[j].u = cell[j].u - (dt/dx)*cell[j].res;
                    cell[j].u = half*(cell[j].u0 + cell[j].u ); //sends things to nan
                    cell[j].w = u2w(cell[j].u);  //Update primitive variables
                }//end do stage02_update

            }

            // Copy the solutions to the ghost cells.
            // In this program, the ghost cell values are used only in the reconstruction.
            cell[0].w        = cell[1].w;
            cell[ncells+1].w = cell[ncells].w;
        } 
        //---------------------------------------------------
        // End of Runge-Kutta Stages
        //---------------------------------------------------
    } 
    //--------------------------------------------------------------------------------
    // End of time stepping
    //--------------------------------------------------------------------------------

}
//********************************************************************************
// End of program
//********************************************************************************



void EulerSolver1D::Solver::initialize( int ncells, float dx, float xmin, const float gamma){
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

        //w2u( cell[i].w, cell[i].u );        //Compute the conservative variables
        cell[i].u = w2u( cell[i].w ); 
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
float EulerSolver1D::Solver::timestep(float cfl, float dx, float gamma, int ncells){
    float dt;             //Output
    //Local variables
    float one = 1.0;
    float u, c, max_speed;
    int   i;

    max_speed = -one;

    for ( int i = 1; i < ncells; ++i ) {
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
 float EulerSolver1D::Solver::minmod(float a, float b){

    float minmod;

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
void EulerSolver1D::Solver::w2u_efficient( Array2D<float>& w, Array2D<float>& u ) {

    u(0) = w(0);
    u(1) = w(0)*w(1);
    u(2) = ( w(2)/(gamma-one) ) + half*w(0)*w(1)*w(1);
    return;
}
Array2D<float> EulerSolver1D::Solver::w2u( Array2D<float>& w) {

    Array2D<float> u(3,1);

    u(0) = w(0);
    u(1) = w(0)*w(1);
    u(2) = ( w(2)/(gamma-one) ) + half*w(0)*w(1)*w(1);
    return u;
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
void EulerSolver1D::Solver::u2w_efficient( Array2D<float>& u, Array2D<float>& w ) {
     
    
    w(0) = u(0);
    w(1) = u(1)/u(0);
    w(2) = (gamma-one)*( u(2) - half*w(0)*w(1)*w(1) );
    return;
}
Array2D<float>  EulerSolver1D::Solver::u2w( Array2D<float>& u ) {
     
    Array2D<float> w(3,1);
    
    w(0) = u(0);
    w(1) = u(1)/u(0);
    w(2) = (gamma-one)*( u(2) - half*w(0)*w(1)*w(1) );
    return w;
}
//--------------------------------------------------------------------------------



//*******************************************************************************
// -- Roe's Flux Function without entropy fix---
//
// P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
// Schemes, Journal of Computational Physics, 43, pp. 357-372. (1981)
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
Array2D<float> EulerSolver1D::Solver::roe_flux(Array2D<float>&  wL, Array2D<float>&  wR){

    //  Input:   wL(3), wR(3) =   Input (conservative variables rho*[1, v, E])
    Array2D<float> flux(3,1);  // Output (numerical flux across L and R states)

    //Local parameters
    float    zero = 0.0;
    float     one = 1.0;
    float    four = 4.0;
    float    half = 0.5;
    float quarter = 0.25;
    //Local variables
    Array2D<float> uL(3,1), uR(3,1);
    float rhoL, rhoR, vL, vR, pL, pR;   // Primitive variables.
    float aL, aR, HL, HR;               // Speeds of sound.
    float RT,rho,v,H,a;                 // Roe-averages
    float drho,du,dP;
    float Da;
    Array2D<float> ws(3,1), R(3,3),dV(3,1);
    int j, k;

    DEBUG("\n w2u \n");
    //w2u(wL,uL);
    //w2u(wR,uR);
    uL = w2u(wL);
    uR = w2u(wR);

    DEBUG("\n left state \n");
//Primitive and other variables.
//  Left state
    rhoL = wL(0);
      vL = wL(1);
      pL = wL(2);
      aL = sqrt(gamma*pL/rhoL);
      HL = ( uL(2) + pL ) / rhoL;
    DEBUG("\n right state \n");
//  Right state
    rhoR = wR(0);
      vR = wR(1);
      pR = wR(2);
      aR = sqrt(gamma*pR/rhoR);
      HR = ( uR(2) + pR ) / rhoR;

    DEBUG("\n Roe Averages \n");
//First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL;
     v = (vL+RT*vR)/(one+RT);
     H = (HL+RT*HR)/(one+RT);
     a = sqrt( (gamma-one)*(H-half*v*v) );

    DEBUG("\n diff primitive variables \n");
//Differences in primitive variables.
   drho = rhoR - rhoL;
     du =   vR - vL;
     dP =   pR - pL;


    DEBUG("\n wave strength characteristic variables \n");
//Wave strength (Characteristic Variables).
   dV(0) =  half*(dP-rho*a*du)/(a*a);
   dV(1) = -( dP/(a*a) - drho );
   dV(2) =  half*(dP+rho*a*du)/(a*a);

    DEBUG("\n abs wave speeds (eigenvalues) \n");
//Absolute values of the wave speeds (Eigenvalues)
   ws(0) = abs(v-a);
   ws(1) = abs(v  );
   ws(2) = abs(v+a);



    DEBUG("\n entropy fix \n");
//Modified wave speeds for nonlinear fields (the so-called entropy fix, which
//is often implemented to remove non-physical expansion shocks).
//There are various ways to implement the entropy fix. This is just one
//example. Try turn this off. The solution may be more accurate.
   Da = max(zero, four*((vR-aR)-(vL-aL)) );
   if (ws(0) < half*Da) { ws(0) = ws(0)*ws(0)/Da + quarter*Da; }
   Da = max(zero, four*((vR+aR)-(vL+aL)) );
   if (ws(2) < half*Da) { ws(2) = ws(2)*ws(2)/Da + quarter*Da; }



    DEBUG("\n Right Eigenvectors \n");
//Right eigenvectors
   R(0,0) = one;
   R(1,0) = v - a;
   R(2,0) = H - v*a;

   R(0,1) = one;
   R(1,1) = v;
   R(2,1) = half*v*v;

   R(0,2) = one;
   R(1,2) = v + a;
   R(2,2) = H + v*a;

    DEBUG("\n Average Flux \n");
//Compute the average flux.
   flux = half*( euler_physical_flux(wL) + euler_physical_flux(wR) );


    DEBUG("\n dissipation term \n");
//Add the matrix dissipation term to complete the Roe flux.
    for (j=0; j<3; ++j){
        for (k=0; k<3; ++k){
            flux(j) = flux(j) - half*ws(k)*dV(k)*R(j,k) ;
        }
   }
    return flux;
} //end roe_flux
//--------------------------------------------------------------------------------

//*******************************************************************************
// Physical flux of the Euler equations (inviscid part only).
//
// ------------------------------------------------------------------------------
//  Input:     w(1:5) = primitive variables (rho, u, p, 0, 0)
//
// Output:  flux(1:5) = physical inviscid flux (rho, rho*u, rho*H, 0, 0)
// ------------------------------------------------------------------------------
//
//*******************************************************************************
//function euler_physical_flux(w) result(flux)
Array2D<float> EulerSolver1D::Solver::euler_physical_flux(Array2D<float>& w){

    Array2D<float> flux(3,1); //Output

    //Local parameters
    const float half = 0.5;
    //Local variables
    float rho, u, p;
    float a2;

    rho = w(0);
      u = w(1);
      p = w(2);

    a2 = gamma*p/rho;

    flux(0) = rho*u;
    flux(1) = rho*u*u + p;
    flux(2) = rho*u*( a2/(gamma-one) + half*u*u ); // H = a2/(gamma-one) + half*u*u
    return flux;
}//end function euler_physical_flux
//-------------------------------------------------------------------------------



//********************************************************************************
//* Write output data file
//*
//* ------------------------------------------------------------------------------
//* Output:  Data file "solution.dat" containing for each cell the following:
//*          cell-center coordinate, density, velocity, pressure, entropy
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
    void EulerSolver1D::Solver::output(){

    float entropy;

    ofstream outfile;
    outfile.open ("solution.dat");
    for (int i=1; i<ncells+1; ++i){
        entropy = log( cell[i].w(2)* pow(cell[i].w(0) , (-gamma)) ) / (gamma-one);
        outfile << std::setprecision(16) << cell[i].xc << '\t'
                << std::setprecision(16) << cell[i].w(0) << '\t' 
                << std::setprecision(16) << cell[i].w(1) << '\t'
                << std::setprecision(16) << cell[i].w(2) << '\t'
                << std::setprecision(16) << entropy <<  "\n";
    }
    outfile.close();

}
//--------------------------------------------------------------------------------


void EulerSolver1D::driverEuler1D(){
    Solver solver;
    solver.Euler1D();
    solver.output();
    return;
}