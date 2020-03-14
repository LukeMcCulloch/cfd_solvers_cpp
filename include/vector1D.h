//
// Specialization of std::vector to CFD unstructured uses
//
//=================================
#ifndef _VECTOR1D_
#define _VECTOR1D_

#include <vector>




//bracket_proxy is parameterized with Vector1D
template <typename ArrayClass, typename Result>
class bracket_proxy_vec{
    public:
        bracket_proxy_vec(ArrayClass& A, int r): A(A), r(r){}

        Result& operator[](int c){return A(r,c); }
    private:
        ArrayClass& A;
        int r;
};



template <class T>
class Vector1D{

public:



   // STORAGE --------------------------------------------------------------- 
   typedef std::vector<T> ContainerType;
   typedef typename ContainerType::iterator Iterator;
   typedef typename ContainerType::const_iterator ConstIterator;
   ContainerType array; // 1D dynamic array



   // CONSTRUCTORS ----------------------------------------------------------
   //Vector();                        // default init to zero
   Vector1D( void );                        // initializes all components to zero
   //Vector( T x, T y, T z); // initializes with specified components
   Vector1D(const std::vector<T>& v );             // initializes from existing vector

   // destructor
   ~Vector1D ();

   // append
   void append( T value);


   // Operators -------------------------------------------------------------
   T& operator() (int i);
   const T&  operator() (int i) const;  


   // array[][] access:
   bracket_proxy_vec<Vector1D, T> operator[](int r){
      return bracket_proxy_vec<Vector1D, T>(*this, r);
   }


};


// CONSTRUCTORS ----------------------------------------------------------------

template <class T>
Vector1D<T>::Vector1D( void ) {
}

template <class T>
Vector1D<T>::Vector1D(const  std::vector<T>& v ) {
   array = v;
}

template<class T>
Vector1D<T>::~Vector1D() {
}

// OPERATIONS ------------------------------------------------------------------

template<class T>
void Vector1D<T>::append(T value) {
   array.push_back(value);
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
