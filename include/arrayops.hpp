

//=================================
// include guard
#ifndef __ARRAYOPS_TEMPLATE_INCLUDED__
#define __ARRAYOPS_TEMPLATE_INCLUDED__

// template <class T>
// class Array2D; //fwd declaration

// #include "array_template.hpp"


// addition of two Array2Ds
template <class T>
Array2D<T> 
operator+(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.nrows == b.nrows);
	assert(a.ncols == b.ncols);

    Array2D<T> result(a.nrows,a.ncols);

    int size = a.storage_size;
    for(size_t i=0; i < size; i++) {
    	result.array[i] = a.array[i] + b.array[i];
    }
    return result;
}



// multiplication of scalar and Array2D
template<typename T>
Array2D<T> 
operator*(T const& s, Array2D<T> const& a)
{
    Array2D<T> result(a.nrows,a.ncols);
    for (size_t k = 0; k<a.size(); ++k) {
        result.array[k] = s*a.array[k];
    }
    return result;
}


// multiplication of two Array2Ds
template <class T>
Array2D<T> 
operator*(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.nrows == b.nrows);
	assert(a.ncols == b.ncols);

    Array2D<T> result(a.nrows,a.ncols);

    int size = a.storage_size;
    for(size_t i=0; i < size; i++) {
    	result.array[i] = a.array[i] * b.array[i];
    }
    return result;
}



// matmul of two Array2D into a third
/**
    this matmul creates a new Array2D, c.
*/
template <class T>
Array2D<T> 
matmul(const Array2D<T>& a, const Array2D<T>& b) {

	assert(a.ncols == b.nrows);

    Array2D<T> c(a.nrows,b.ncols);
    c = 0.;

    for( size_t i=0; i<a.nrows; ++i ) {
        for( size_t k=0; k<b.ncols; ++k ) {
            c(i,k) = a(i,0) * b(0,k);
        }
        for( size_t j=1; j<a.ncols; ++j ) {
            for( size_t  k=0; k<b.ncols; ++k ) {
                c(i,k) += a(i,j) * b(j,k);
            }
        }
    }


    // for (size_t i = 0; i < a.nrows; i++) {
    //     for (size_t k = 0; k < a.ncols; k++) {
    //         for (size_t j = 0; j < b.ncols; j++) {
    //             c(i, j) += a(i, k) * b(k, j);
    //         }
    //     }
    // }
    return c;
}



// matmul of two Array2D into a third
/**
    this matmul over-writes Array2D c
*/
template <class T>
Array2D<T> 
matmul(const Array2D<T>& a, const Array2D<T>& b, Array2D<T>& c) {
    
	assert(a.ncols == b.nrows);
    assert(a.nrows == c.nrows);
	assert(b.ncols == c.ncols);

    c = 0.;

    for( size_t i=0; i<a.nrows; ++i ) {
        for( size_t k=0; k<b.ncols; ++k ) {
            c(i,k) = a(i,0) * b(0,k);
        }
        for( size_t j=1; j<a.ncols; ++j ) {
            for( size_t  k=0; k<b.ncols; ++k ) {
                c(i,k) += a(i,j) * b(j,k);
            }
        }
    }

    // for (size_t i = 0; i < a.nrows; i++) {
    //     for (size_t k = 0; k < a.ncols; k++) {
    //         for (size_t j = 0; j < b.ncols; j++) {
    //             c(i, j) += a(i, k) * b(k, j);
    //         }
    //     }
    // }
    return c;
}


// multiplication of Array2D and scalar
// addition of scalar and Array2D
// addition of Array2D and scalar
//...



// Dot of two Array2Ds
// this matmul creates a new Array2D, c.
template <class T>
Array2D<T> 
Dot(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.ncols == b.nrows);
    Array2D<T> c(a.nrows,b.ncols);
    c = 0.;
    for( size_t i=0; i<a.nrows; ++i ) {
        for( size_t k=0; k<b.ncols; ++k ) {
            c(i,k) = a(i,0) * b(0,k);
        }
        for( size_t j=1; j<a.ncols; ++j ) {
            for( size_t  k=0; k<b.ncols; ++k ) {
                c(i,k) += a(i,j) * b(j,k);
            }
        }
    }
    return c;
}

// Cross two 1D Array2Ds
// this matmul creates a new 1D Array2D, c.
template <class T>
Array2D<T> 
Cross(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.ncols == b.nrows);
	assert(a.ncols == 1);
	assert(b.nrows == 1);

    Array2D<T> c(a.nrows,1);

    c(0,0) =   a(1,0)*b(2,0) - a(2,0)*b(1,0);
    c(1,0) = -(a(0,0)*b(2,0) - a(2,0)*b(0,0));
    c(2,0) =   a(0,0)*b(1,0) - a(1,0)*b(0,0);

    return c;
}



#endif //__ARRAYOPS_TEMPLATE_INCLUDED__