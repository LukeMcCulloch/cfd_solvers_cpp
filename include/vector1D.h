//
/ /Specialization of std::vector to CFD unstructured uses
//
//=================================
#ifndef _VECTOR1D_
#define _VECTOR1D_

#include <vector>

template <class T>
class Vector1D{

public:
   // CONSTRUCTORS ----------------------------------------------------------
   //Vector();                        // default init to zero
   Vector1D( void );                        // initializes all components to zero
   //Vector( T x, T y, T z); // initializes with specified components
   Vector1D(const std::vector<T>& v );             // initializes from existing vector

   // destructor
   ~Vector1D ();

   // append
   push_back( T value);


   // Operators -------------------------------------------------------------
   T& operator() (int i);
   const T&  operator() (int i) const;  

   // STORAGE ---------------------------------------------------------------
   std::vector<T> array; // 1D dynamic array

};


// CONSTRUCTORS ----------------------------------------------------------------

template <class T>
Vector1D<T>::Vector1D( void ) {
}

template <class T>
Vector1D<T>::Vector1D(const  std::vector<T>& v ) {
   arrayy = v;
}

template<class T>
Vector1D<T>::~Vector1D() {
}

// OPERATIONS ------------------------------------------------------------------

template<class T>
Vector1D<T>
push_back(T value) {
   array.push_back(value)
}

template<class T>
T& Vector1D<T>::operator() (int i) {
   return array[i];
}

template<class T>
const T& Vector1D<T>::operator() (int i) const {
   return array[i];
}

#endif
