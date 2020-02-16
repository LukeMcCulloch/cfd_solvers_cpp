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
//*   (see "I do like CFD, VOL.1", whose PDF can be downloaded at cfdbooks.com).
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to c++ by Luke McCulloch, PhD.
//*
//*
//* This is Version 0 (June 2019 translation)
//*
//* ------------------------------------------------------------------------------
//* Files: There were 3 files (in the fortran version).
//* This translated version will differ.
//*
//* ------------------------------------------
//* - Main driver program file   : This reads grid and BC files, and dummy NC/CC programs.
//*
//*     edu2d_euler_rk2_main.f90 (This file), which contains a main driver program
//*      -- edu2d_euler_rk2       : Main driver code, which  an Euler solver
//*
//* ------------------------------------------
//* - Basic EDU2D package file   : Arranged for a 2D Euler code
//*
//*     edu2d_basic_package_euler_rk2.f90, which contains the following modules.
//*      -- EulerSolver2D      : Numerical values defined
//*      -- _data_type : Grid data types defined
//*      -- edu2d_main_data      : Main grid data and parameters declared
//*      -- _data      : Read/construct/check grid data
//*
//* ------------------------------------------
//* - Euler solver file   : This computes a solution to the shock diffraction problem.
//*
//*     edu2d_euler_rk2_main.f90, which contains a 2D Euler solver with RK2.
//*      -- edu2d_euler_rk2_solver : Node-centered Explicit Euler solver with RK2
//*
//* ------------------------------------------------------------------------------
//* Notes:
//*
//*  The purpose of this code is to give a beginner an opportunity to learn how to
//*  write an unstructured CFD code. Hence, the focus here is on the simplicity.
//*  The code is not optimized for efficiency.
//*
//*  This code is set up specifi for a shock diffraction problem.
//*  It can be modified easily to solve other problems:
//*
//*   1. Define your own free stream Mach number, M_inf, at the beginning of the main.
//*   2. Modify the subroutine, initial_solution_shock_diffraction(), which is in this
//*      file, to set up an appropriate initial condition for your problem.
//*   3. Delete the special treatment at the corner in euler_solver_main.f90
//*      (Special treatment is done in two places in that file.)
//*
//*  If the code is not simple enough to understand, please send questions to Hiro
//*  at sunmasen(at)hotmail.com. I'll greatly appreciate it and revise the code.
//*
//*  If the code helps you understand how to write your own code that is more
//*  efficient and has more features, it'll have served its purpose.
//*
//* ------------------------------------------------------------------------------
//* Examples of additional features you might want to add.
//*
//*  1. Local time-stepping      (faster convergence for steady problems)
//*  2. Implicit time-stepping   (remove CFL restriction)
//*  3. More boundary conditions (periodic, symmetry, suction, etc.)
//*  4. Other reconstruction     (Van Leer's kappa scheme)
//*  5. Other limiters           (Venkat/Barth limiter,etc.)
//*  6. Other flux functions     (HLL, LDFSS, AUSM, etc.)
//*  7. Local-preconditioning    (low-Mach number accuracy and stiffness removal)
//*  8. CFL ramping              (increase CFL gradually for a stable start-up)
//*  9. More output              (convergence history, etc.)
//* 10. Parallelization          (large-scale problems)
//* 11. Grid adaptation          (h-refinement, steady or unsteady)
//* 12. Adjoint capability       (aerodynamic design, adaptation, etc.)
//* 13. Moving grid              (sliding mesh, overset grid, etc.)
//* 14. Multigrid                (grid-independent convergence)
//* ------------------------------------------------------------------------------
//*
//* Katate Masatsuka, http://www.cfdbooks.com
//********************************************************************************

//********************************************************************************
//* Main program: Node-centered finite-volume Euler code
//*
//* This code computes an unsteady solution of the Euler equations.
//* It is set up to solve a shock diffraction problem.
//* So, it is not really (but almost) a general purpose code.
//*
//* Input -------------------------------------------------------
//*
//*   project.grid  = grid file containing boundary information
//*   project.bcmap = file that specifies boundary conditions
//*
//*   (Note: See the subroutine "read_grid", which is in this file,
//*          for the format of these files.)
//*
//*  Parameters to be specified inside the main program:
//*
//*          M_inf = Upstream Mach number
//*          gamma = Ratio of specific heats (1.4 for air)
//*            CFL = CFL number (global time step)
//*        t_final = Final time to stop the calculation
//*  time_step_max = Max iterations (just a big enough number)
//*  inviscid_flux = Inviscid flux selection (1 = Roe, 2 = Rotated-RHLL)
//*   limiter_type = Limiter selection (1=Van Albada limiter, 0=No limiter)
//*
//*
//* Output ------------------------------------------------------
//*
//*  "project_tecplot.dat" = Tecplot file containing the grid and the solution.
//*
//*
//*  NOTE: The variables are nondimensionalized values (compressible scale),
//*           rho=rho/rho_inf, u=u/a_inf, v=v/a_inf, rho=p/(rho_inf*a_inf^2)
//*
//*  NOTE: Many variables are passed to subroutines via USE statement.
//*        Each module contains specific data, and they are accessed by USE.
//*
//*
//* Katate Masatsuka, http://www.cfdbooks.com
//* Translated to c++ by Luke McCulloch
//*
//*
//*
//********************************************************************************
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
// 2D Eiuler approximate Riemann sovler
#include "../include/EulerUnsteady2D.h"
#include "../include/EulerUnsteady2D_basic_package.h"

//======================================
//using namespace std;

//======================================
//
//namespace EulerSolver2D{



EulerSolver2D::MainData2D::MainData2D(std::string datafile_grid_in, 
                                               std::string datafile_bcmap_in) {


   M_inf             = 0.0;                  // Freestream Mach number to be set in the function
                                             //    -> "initial_solution_shock_diffraction"
                                             //    (Specify M_inf here for other problems.)
   gamma             = 1.4;                  // Ratio of specific heats
   CFL               = 0.95;                 // CFL number
   t_final           = 0.18;                 // Final time to stop the calculation.
   time_step_max     = 5000;                 // Max time steps (just a big enough number)
   inviscid_flux     = "rhll";               // = Rotated-RHLL      , "roe"  = Roe flux
   limiter_type      = "vanalbada";          // = Van Albada limiter, "none" = No limiter
   nq        = 4;                    // The number of equtaions/variables in the target equtaion.
   gradient_type     = "linear";             // or "quadratic2 for a quadratic LSQ.
   gradient_weight   = "none";               // or "inverse_distance"
   gradient_weight_p =  EulerSolver2D::one;  // or any other real value


   // (1) Read grid files
   read_grid(datafile_grid_in, datafile_bcmap_in);

   std::cout << "Allocate arrays" << std::endl;
   std::cout << "there are " << nnodes << " nodes " << std::endl;

   for (size_t i = 0; i < nnodes; i++) {
      std::cout << "i = " << i << " of " << nnodes << std::endl;
      //std::cout << node[i].x << std::endl;

    // declare within the class
      // node[i].u     = new Array2D<real>(E2Ddata.nq,1);
      // node[i].du    = new Array2D<real>(E2Ddata.nq,1);
      // node[i].w     = new Array2D<real>(E2Ddata.nq,1);
      // node[i].gradw = new Array2D<real>(E2Ddata.nq,2); //<- 2: x and y components.
      // node[i].res   = new Array2D<real>(E2Ddata.nq,1);
      //SE2Ddata.node[i].cell
   }

   std::cout << "E2Ddata.nq, = " << nq << std::endl;
// (2) Construct grid data
   construct_grid_data();
}




EulerSolver2D::MainData2D::~MainData2D() {
   delete[] node;
}

EulerSolver2D::Solver::Solver(){}


EulerSolver2D::Solver::~Solver(){
    //delete[] cell;
    }


//********************************************************************************
//********************************************************************************
//********************************************************************************
//* Euler solver: Node-Centered Finite-Volume Method (Edge-Based)
//*
//* - Node-centered finite-volume method for unstructured grids(quad/tri/mixed)
//* - Roe flux with an entropy fix and Rotated-RHLL flux
//* - Reconstruction by unweighted least-squares method (2x2 system for gradients)
//* - Van Albada slope limiter to the primitive variable gradients
//* - 2-Stage Runge-Kutta time-stepping
//*
//********************************************************************************
void EulerSolver2D::Solver::euler_solver_main(){

}
//********************************************************************************
// End of program
//********************************************************************************

//namespace EulerSolver2D{
void program_2D_euler_rk2(){
   // procedural fortran ends up in this function
   int i;

   // euler solver 2D:
   EulerSolver2D::Solver E2Dsolver;

   //Set file names, Inout data files
   std::string  datafile_grid_in  = "project.grid";  //Grid file
   std::string  datafile_bcmap_in = "project.bcmap"; //Boundary condition file
   std::string  datafile_tec      = "project_tecplot.dat";  //Tecplot file for viewing the result.


//--------------------------------------------------------------------------------
// Input Parameters

   // euler main data:
   //typedef EulerSolver2D::EulerSolver2D 2Ddata;
   EulerSolver2D::MainData2D E2Ddata(datafile_grid_in, datafile_bcmap_in);
   //2Ddata = new EulerSolver2D();

   //              E2Ddata.M_inf  = 0.0;         // Freestream Mach number to be set in the function
   //                                  //    -> "initial_solution_shock_diffraction"
   //                                  //    (Specify M_inf here for other problems.)
   //              E2Ddata.gamma = 1.4;         // Ratio of specific heats
   //                E2Ddata.CFL = 0.95;        // CFL number
   //            E2Ddata.t_final = 0.18;        // Final time to stop the calculation.
   //      E2Ddata.time_step_max = 5000;        // Max time steps (just a big enough number)
   //      E2Ddata.inviscid_flux = "rhll";      // = Rotated-RHLL      , "roe"  = Roe flux
   //       E2Ddata.limiter_type = "vanalbada"; // = Van Albada limiter, "none" = No limiter
   //                 E2Ddata.nq = 4;           // The number of equtaions/variables in the target equtaion.
   //  E2Ddata.gradient_type     = "linear";    // or "quadratic2 for a quadratic LSQ.
   //  E2Ddata.gradient_weight   = "none";      // or "inverse_distance"
   //  E2Ddata.gradient_weight_p =  EulerSolver2D::one;        // or any other real value
//--------------------------------------------------------------------------------
// Solve the Euler equations and write the output datafile.
//
// // (1) Read grid files
//    E2Ddata.read_grid(datafile_grid_in, datafile_bcmap_in);

//    std::cout << "Allocate arrays" << std::endl;
//    std::cout << "there are " << E2Ddata.nnodes << " nodes " << std::endl;

//    for (size_t i = 0; i < E2Ddata.nnodes; i++) {
//       std::cout << "i = " << i << " of " << E2Ddata.nnodes << std::endl;
//       //std::cout << E2Ddata.node[i].x << std::endl;

//     // declare within the class
//       // E2Ddata.node[i].u     = new Array2D<real>(E2Ddata.nq,1);
//       // E2Ddata.node[i].du    = new Array2D<real>(E2Ddata.nq,1);
//       // E2Ddata.node[i].w     = new Array2D<real>(E2Ddata.nq,1);
//       // E2Ddata.node[i].gradw = new Array2D<real>(E2Ddata.nq,2); //<- 2: x and y components.
//       // E2Ddata.node[i].res   = new Array2D<real>(E2Ddata.nq,1);
//       //SE2Ddata.node[i].cell
//    }

//    std::cout << "E2Ddata.nq, = " << E2Ddata.nq << std::endl;
// // (2) Construct grid data
//    E2Ddata.construct_grid_data();

// (3) Check the grid data (It is always good to check them before use//)
//       check_grid_data();

// // (4) Prepare LSQ gradients
//       compute_lsq_coeff_nc();
//       check_lsq_coeff_nc();

// // (5) Set initial solution for a shock diffraction problem
// //     (Re-write or replace it by your own subroutine for other problems.)
//       initial_solution_shock_diffraction();

// // (6) Compute the solution (March in time to the final time)
//       euler_solver_main();

// // (7) Write out the tecplot data file (Solutions at nodes)
//       write_tecplot_file(datafile_tec);

}

void EulerSolver2D::driverEuler2D(){
    program_2D_euler_rk2();
    //solver.output();
    return;
}

//}  //end namespace EulerSolver2D