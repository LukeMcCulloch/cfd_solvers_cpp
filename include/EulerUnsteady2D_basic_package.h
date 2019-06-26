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
//
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
                         
     real                   zero = 0.0,
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

#endif