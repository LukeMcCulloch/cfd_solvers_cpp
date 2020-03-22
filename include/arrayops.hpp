

//=================================
// include guard
#ifndef __ARRAYOPS_TEMPLATE_INCLUDED__
#define __ARRAYOPS_TEMPLATE_INCLUDED__

// template <class T>
// class Array2D; //fwd declaration

// #include "array_template.hpp"

#include <cmath>

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



// subtraction of two Array2Ds
template <class T>
Array2D<T> 
operator-(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.nrows == b.nrows);
	assert(a.ncols == b.ncols);

    Array2D<T> result(a.nrows,a.ncols);

    int size = a.storage_size;
    for(size_t i=0; i < size; i++) {
    	result.array[i] = a.array[i] - b.array[i];
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


// division of an Array2D by a scalar
template<typename T>
Array2D<T> 
operator/(Array2D<T> const& a , T const& s )
{
    Array2D<T> result(a.nrows,a.ncols);
    for (size_t k = 0; k<a.size(); ++k) {
        result.array[k] = a.array[k]/s;
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


// Makes a m x m matrix with all diagonal elements to 1, and other elements
// to 0.
// template <class T>
// Array2D<T> 
// makeIdentity(size_t m) {
//     Array2D<T> c(m,m);
//     c = 0.;
//     for(size_t i = 0; i < c.nrows; ++i ) {
//         c(i,i) = 1.0;
//     }
//     return c;
// }




template <class T>
T Dotscalar(const Array2D<T>& a, const Array2D<T>& b) {
	assert(a.ncols == 1);
	assert(b.ncols == 1);
	assert(a.nrows == b.nrows);

    T c;
    c = 0.;
    for( size_t i=0; i<a.nrows; ++i ) {
        for( size_t k=0; k<b.nrows; ++k ) {
            c += a(i,0) * b(k),0;
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


//
//***********************************************************
// misc functions
template <class T>
T MaxColVal(Array2D<T> A, int col) {
  T max_ = 0.0;
  for (int i = 0; i < A.nrows; i++) {
    max_ = max(max_,A(i,col));
  }
  return max_;
}


template <class T>
Array2D<T>
abs( Array2D<T>& A ) {
    Array2D<T> result(A.nrows,A.ncols);
    for (int i = 0; i < A.storage_size; i++) {
        result.array[i] = abs( A.array[i] );
    }
    return result;
}



// inverse


// Gauss Siedel Solve Am = b
/**
*    a m = b
*    n = inv(a)*b
* 
*/
template <class T>
Array2D<T> 
GaussSeidelInv(const Array2D<T>& a, 
                        Array2D<T>& m, 
                        const Array2D<T>& b) {
    
	assert(a.ncols == b.nrows);
    assert(a.nrows == m.nrows);
	assert(b.ncols == m.ncols);

    int p = a.nrows; // a = pxp matrix


    Array2D<T> n = m; // bad sizing.  The ans size = size(a)
    n = 0;
    Array2D<T> np1 = n;

    int q   = 100;
    T convg = 10.;
    T tol   = 1.0;


    while (q > 0 & convg > tol) {
        for (size_t i = 0; i < p; i++) {
            //np1 = m;
            n(i) = b(i) / a(i,i);
            for (size_t j = 0; j < p; j++) {
                if (j == i){
                    continue;
                }
                np1(i) = -m(i);
                  n(i) =  n(i) - ( ( a(i,j) / a(i,i) ) * m(j) );
                  m(i) =  n(i);
            }
            //cout<<"x"<<i + 1 << "="<<n[i]<<" ";
        }
        //cout << "\n";
        convg = Dotscalar(n,np1);
        q--;
    }
    if (q==0) n.istat = 1;
    return n;
}


template <typename T>
Array2D<T> Array2D<T>::invert() {
    //JET_ASSERT(isSquare());
    if (not isSquare()) {
        cout << "ERROR, trying to directly invert nonsquare matrix " << endl;
        std::exit(0);
    }

    // Computes inverse matrix using Gaussian elimination method.
    // https://martin-thoma.com/solving-linear-equations-with-gaussian-elimination/
    size_t n = getnrows();
    Array2D& a = *this;
    //Array2D<T> rhs = makeIdentity(n);
    // Array2D<T> rhs(n,n); 
    // rhs = makeIdentity(n);
    Array2D<T> rhs(true, n, n);

    for (size_t i = 0; i < n; ++i) {
        // Search for maximum in this column
        T maxEl = std::fabs(a(i, i));
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::fabs(a(k, i)) > maxEl) {
                maxEl = std::fabs(a(k, i));
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        if (maxRow != i) {
            for (size_t k = i; k < n; ++k) {
                std::swap(a(maxRow, k), a(i, k));
                std::swap(rhs(maxRow, k), rhs(i, k));
            }
        }

        // Make all rows except this one 0 in current column
        for (size_t k = 0; k < n; ++k) {
            if (k == i) {
                continue;
            }
            T c = -a(k, i) / a(i, i);
            for (size_t j = 0; j < n; ++j) {
                rhs(k, j) += c * rhs(i, j);
                if (i == j) {
                    a(k, j) = 0;
                } else if (i < j) {
                    a(k, j) += c * a(i, j);
                }
            }
        }

        // Scale
        for (size_t k = 0; k < n; ++k) {
            T c = 1 / a(k, k);
            for (size_t j = 0; j < n; ++j) {
                a(k, j) *= c;
                rhs(k, j) *= c;
            }
        }
    }

    //set(rhs);
    return rhs;
}


#endif //__ARRAYOPS_TEMPLATE_INCLUDED__