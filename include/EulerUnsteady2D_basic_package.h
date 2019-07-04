//********************************************************************************
//* Educationally-Designed Unstructured 2D (EDU2D) Code
//*
//*
//*
//*     This module belongs to the inviscid version: EDU2D-Euler-RK2
//*
//*
//*
//* This file contains 4 modules:
//*
//*  1. module edu2d_constants      - Some numerical values, e.g., zero, one, pi, etc.
//*  2. module edu2d_grid_data_type - Grid data types: node, edge, face, element, etc.
//*  3. module edu2d_my_main_data   - Parameters and arrays mainly used in a solver.
//*  4. module edu2d_my_allocation  - Subroutines for dynamic allocation
//*  5. module edu2d_grid_data      - Subroutines for reading/constructing/checking grid data
//*
//* All data in the modules can be accessed by the use statement, e.g., 'use constants'.
//*
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to C++ by Luke McCulloch, PhD.
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015).
//* This F90 code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************

//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//* 1. module edu2d_constants
//*
//* Some useful constants are defined here.
//* They can be accessed by the use statement, 'use constants'.
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to C++ by Luke McCulloch, PhD.
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015). -> C++ June 2019
//* In the spirit of its fortran parent,
//* This C++ code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************
//=================================
// include guard
#ifndef __eulerUnsteady2d_basic_package_INCLUDED__
#define __eulerUnsteady2d_basic_package_INCLUDED__

//======================================
//stl
#include <vector> 
//======================================
// my simple array class template (type)
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"

//#define REAL_IS_DOUBLE true
#ifdef REAL_IS_DOUBLE
  typedef double real;
#else
  typedef float real;
#endif
//
//********************************************************************************
//
struct edu2d_constants 
{
  edu2d_constants() = default; // asks the compiler to generate the default implementation
                         
  real                 zero = 0.0,
                        one = 1.0,
                        two = 2.0,
                      three = 3.0,
                       four = 4.0,
                       five = 5.0,
                        six = 6.0,
                      seven = 7.0,
                      eight = 8.0,
                       nine = 9.0,
                        ten = 10.0,
                     eleven = 11.0,
                       half = 0.5,
                      third = 1.0 / 3.0,
                     fourth = 1.0 / 4.0,
                      fifth = 1.0 / 5.0,
                      sixth = 1.0 / 6.0,
                  two_third = 2.0 / 3.0,
                 four_third = 4.0 / 3.0,
               three_fourth = 3.0 / 4.0,
                    twelfth = 1.0 /12.0,
           one_twentyfourth = 1.0 /24.0;

          real pi = 3.141592653589793238;
};


//********************************************************************************

//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//* 2. module grid_data_type
//*
//* This module defines custom grid data types for unstructured grids
//*
//* NOTE: These data types are designed to make it easier to understand the code.
//*       They may not be the best in terms of efficiency.
//*
//* NOTE: Custom grid data types (derived types) are very useful.
//*       For example, if I declare a variable, "a", by the statemant:
//*           type(node_type), dimension(100) :: a
//*       The variable, a, is a 1D array each component of which contains all data
//*       defined as below. These data can be accessed by %, e.g.,
//*           a(1)%x, a(1)%y, a(1)%nghbr(1:nnghbrs), etc.
//*       In C-programming, this type of data is called "structure", I think.
//*
//*post
//*
//* NOTE: Not all data types below are used in EDU2D-Euler-RK2.
//*
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to C++ by Luke McCulloch
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015).  Translated in June/July 2019
//* This C++ code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************
namespace edu2d_grid_data_type{


//----------------------------------------------------------
// Data type for nodal quantities (used for node-centered schemes)
// Note: Each node has the following data.
//----------------------------------------------------------
  class node_type{

    //  to be read from a grid file
    real x, y;                  //nodal coordinates
    //  to be constructed in the code

    int nnghbrs;   //number of neighbors
    //int,   dimension(:), pointer  :: nghbr     //list of neighbors
    vector<int> nghbr;          //list of neighbors

    int nelms;     //number of elements
    //Array2D<int>*  elm        //list of elements
    vector<int> elm;            //list of elements

    real vol;                   //dual-cell volume
    int bmark;                  //Boundary mark
    int nbmarks;                //# of boundary marks
    //  to be computed in the code
    //Below are arrays always allocated.
    Array2D<real>* u;           //conservative variables
    Array2D<real>* uexact;      //conservative variables
    //Array2D<real>* gradu      //gradient of u
    Array2D<real>* gradu;       //gradient of u (2D) 
    Array2D<real>* res;         //residual (rhs)
    real ar;                    // Control volume aspect ratio
    Array2D<real>* lsq2x2_cx;   //    Linear LSQ coefficient for ux
    Array2D<real>* lsq2x2_cy;   //    Linear LSQ coefficient for uy
    Array2D<real>* lsq5x5_cx;   // Quadratic LSQ coefficient for ux
    Array2D<real>* lsq5x5_cy;   // Quadratic LSQ coefficient for uy
    Array2D<real>* dx, dy;      // Extra data used by Quadratic LSQ
    //Array2D<real>* dw         // Extra data used by Quadratic LSQ
    Array2D<real>* dw;          // Extra data used by Quadratic LSQ (2D)

    //Below are optional: Pointers need to be allocated in the main program if necessary.
    Array2D<real>* du;          //change in conservative variables
    Array2D<real>* w;           //primitive variables(optional)
    //Array2D<real>* gradw      //gradient of w
    Array2D<real>* gradw;       //gradient of w (2D)
    real phi;                   //limiter function (0 <= phi <= 1)
    real dt;                    //local time step
    real wsn;                   //Half the max wave speed at face
    Array2D<real>*  r_temp;     // For GCR implementation
    Array2D<real>*  u_temp;     // For GCR implementation
    Array2D<real>*  w_temp;     // For GCR implementation
};

//----------------------------------------------------------
// Data type for element/cell quantities (used for cell-centered schemes)
// Note: Each element has the following data.
//----------------------------------------------------------
  class elm_type{
    //  to be read from a grid file
    int nvtx;                   //number of vertices
    Array2D<int>*  vtx;         //list of vertices
    //  to be constructed in the code
    int nnghbrs;                //number of neighbors
    Array2D<int>*  nghbr;       //list of neighbors
    real x, y;                  //cell center coordinates
    real vol;                   //cell volume

    Array2D<int>*  edge;        //list of edges
    Array2D<real>* u;            //conservative variables
    Array2D<real>* uexact;       //conservative variables
    //NotUsed   Array2D<real>* du       //change in conservative variables
    Array2D<real>* gradu;        //gradient of u
    Array2D<real>* res;          //residual (rhs)
    real dt;                     //local time step
    real wsn;                    //??
    int bmark;                   //Boundary mark
    int nvnghbrs;                //number of vertex neighbors
    Array2D<int>* vnghbr;        //list of vertex neighbors
    real ar;                     //Element volume aspect ratio
    Array2D<real>* lsq2x2_cx;    //Linear LSQ coefficient for ux
    Array2D<real>* lsq2x2_cy;    //Linear LSQ coefficient for uy

};

//----------------------------------------------------------
// Data type for edge quantities (used for node-centered scheemes)
// Note: Each edge has the following data.
//----------------------------------------------------------
  class edge_type{
//  to be constructed in the code
   int n1, n2;            //associated nodes
   int e1, e2;            //associated elements
   Array2D<real>* dav;    //unit directed-area vector
   real           da;     //magnitude of the directed-area vector
   Array2D<real>* ev;     //unit edge vector
   real            e;     //magnitude of the edge vector
   int kth_nghbr_of_1;    //neighbor index
   int kth_nghbr_of_2;    //neighbor index

  };

//----------------------------------------------------------
// Data type for boundary quantities (for both node/cell-centered schemes)
// Note: Each boundary segment has the following data.
//----------------------------------------------------------
  class bgrid_type{
    //  to be read from a boundary grid file
    char bc_type[80];     //type of boundary condition
    int nbnodes; //# of boundary nodes
    Array2D<int>* bnode;  //list of boundary nodes
    //  to be constructed in the code
    int nbfaces; //# of boundary faces
    Array2D<real>* bfnx;  //x-component of the face outward normal
    Array2D<real>* bfny;  //y-component of the face outward normal
    Array2D<real>* bfn;   //magnitude of the face normal vector
    Array2D<real>* bnx;   //x-component of the outward normal
    Array2D<real>* bny;   //y-component of the outward normal
    Array2D<real>* bn;    //magnitude of the normal vector
    Array2D<int>*  belm;  //list of elm adjacent to boundary face
    Array2D<int>*  kth_nghbr_of_1;
    Array2D<int>*  kth_nghbr_of_2;

  };

//----------------------------------------------------------
// Data type for face quantities (used for cell-centered schemes)
//
// A face is defined by a line segment connecting two nodes.
// The directed area is defined as a normal vector to the face,
// pointing in the direction from e1 to e2.
//
//      n2
//       o------------o
//     .  \         .
//    .    \   e2  .
//   .  e1  \    .
//  .        \ .         Directed area is positive: n1 -> n2
// o----------o         e1: left element
//             n1       e2: right element (e2 > e1 or e2 = 0)
//
// Note: Each face has the following data.
//----------------------------------------------------------
  class face_type{
// to be constructed in the code (NB: boundary faces are excluded.)
   int n1, n2;    //associated nodes
   int e1, e2;    //associated elements
   real da;       //magnitude of the directed-area vector
   Array2D<real> dav = Array2D<real>(2,1);   //unit directed-area vector

};

}

#endif