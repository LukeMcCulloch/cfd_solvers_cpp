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
//*        translated by Luke McCulloch
//*
//*
//* This is Version 0 (June 2019 translation)
//*
//* ------------------------------------------------------------------------------
//* Files: There are 3 files.
//*
//* ------------------------------------------
//* - Main driver program file   : This reads grid and BC files, and call dummy NC/CC programs.
//*
//*     edu2d_euler_rk2_main.f90 (This file), which contains a main driver program
//*      -- edu2d_euler_rk2       : Main driver code, which calls an Euler solver
//*
//* ------------------------------------------
//* - Basic EDU2D package file   : Arranged for a 2D Euler code
//*
//*     edu2d_basic_package_euler_rk2.f90, which contains the following modules.
//*      -- edu2d_constants      : Numerical values defined
//*      -- edu2d_grid_data_type : Grid data types defined
//*      -- edu2d_main_data      : Main grid data and parameters declared
//*      -- edu2d_grid_data      : Read/construct/check grid data
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
//*  This code is set up specifically for a shock diffraction problem.
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
//* Translated ti c++ by Luke McCulloch
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
//#include "../include/.h"

//======================================
using namespace std;

//======================================
//fwd declarations
struct cell_data; 
class Solver;
