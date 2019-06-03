/*
Specialization of Array2D to vector
*/

//=================================
#ifndef Special_VECTOR_H
#define Special_VECTOR_H

//=================================
#include <cassert>
#include <iostream>
#include <limits>


//=================================
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"

template <class T>
class Vector{

public:
    // CONSTRUCTORS ----------------------------------------------------------
    //Vector();                        // default init to zero
    Vector( void );                        // initializes all components to zero
    Vector( T x, T y, T z); // initializes with specified components
    Vector(const Vector& v );             // initializes from existing vector
   

    // destructor
    ~Vector ();

    // STORAGE ---------------------------------------------------------------
    T x, y, z; // components

};


// CONSTRUCTORS ----------------------------------------------------------------

// template <class T>
// Vector<T> :: Vector( ){
//     x=0.;
//     y=0.;
//     z=0.;
//     printf("0");
// }

template <class T>
Vector<T>::Vector( void )
// initializes all components to zero
: x( 0. ),
  y( 0. ),
  z( 0. )
{
    printf("\n1");
};


template <class T>
Vector<T>::Vector( T x0,
                   T y0,
                   T z0 )
// initializes with specified components
: x( x0 ),
  y( y0 ),
  z( z0 )
{
    printf("\n2");
};


template <class T>
Vector<T>::Vector(const  Vector& v )
// initializes from existing vector
: x( v.x ),
  y( v.y ),
  z( v.z )
{
    printf("\n3\n");
};


template<class T>
Vector<T>::~Vector(){
}

#endif