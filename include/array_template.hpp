//=================================
// include guard
#ifndef __ARRAY_TEMPLATE_INCLUDED__
#define __ARRAY_TEMPLATE_INCLUDED__



//=================================
// array class
#include <cassert>
#include <iostream>
#include <limits>
#include <string.h>




//=================================
//using namespace std; // do not do this in headers
using std::cout;
using std::cin;
using std::endl;
using std::numeric_limits;
using std::ofstream;
using std::max;





template <class T>
class Array2D; //fwd declaration



// undefine to disable range checking
#define RANGE_CHECK


class overdetermined : public std::domain_error{

   public:
     overdetermined() 
        : std::domain_error("solution is over-determined")
     {}
};

class underdetermined : public std::domain_error{

   public:
      underdetermined() 
         : std::domain_error("solution is under-determined")
      {}
};



//bracket_proxy is parameterized with Array2D
template <typename ArrayClass, typename Result>
class bracket_proxy{
    public:
        bracket_proxy(ArrayClass& A, int r): A(A), r(r){}

        Result& operator[](int c){return A(r,c); }
        //const Result& operator[](int c){return const A(r,c); }
    private:
        ArrayClass& A;
        int r;
};






template <class T>
class Array2D{

    void check_indices(int i, int j) const{
        assert(i >= 0 && i < nrows);
        assert(j >= 0 && j < ncols);
    }

    void check_size(int other_nrows,int other_ncols) const { 
        assert(nrows == other_nrows);
        assert(ncols == other_ncols); 
    }

    void check_index(int i) const { assert(i >= 0 && i < nrows); }

    void check_size(int i) const {
        //printf("\nchecl size i = %d, storage_size = %d\n",i,storage_size);
        assert(i >= 0 && i < storage_size);
    }

    private:
        int nBytes;

    public:
        // set up data size of matrix
        int nrows, ncols;
        int storage_size;
        int tracked_index;
        int istat = 0;
        bool allocated;
        bool makeidentidy;

        // malloc host memory
        T* array;
        
        // explicit constructor declaring size nrow,ncol:
        explicit Array2D(int numrows, int numcols): 
                        nrows(numrows), ncols(numcols){
            //cout << "building \n" << endl;
            //buildWithParameters(numrows, numcols);
            tracked_index = 0;
            build();
            allocated = true;
            //cout << "built \n" << endl;
            //initialize();
            //cout << "initialized \n" << endl;
        }
            
        Array2D();
        Array2D(bool makeidentidy, size_t m, size_t n);
        
        //copy constructor 
        Array2D(const Array2D& A);


        // destructor
        ~Array2D ();


        // resize! constructor
        Array2D(const Array2D& A, int nrows, int ncols);

        

        //methods:
        //
        //setting
        void build();
        void buildWithParameters(int, int);
        //void initialize();
        void setonce( T data);
        //
        // inversion
        Array2D<T>  invert();
        //
        //printing
        void print();
        std::ostringstream& print(std::ostringstream&);

       
        void cache() const {}

        // Array2D& val() const {
        //     return *this;
        // }
        Array2D<T> val() const {
            return *this;
        }


        // Returns true if this matrix is a square matrix.
        bool isSquare() const;

        //
        size_t getnrows() const {
            return nrows;
        }
        size_t getncols() const {
            return ncols;
        }
        // return size
        size_t size() const {
            return storage_size;
        }
        
        // Operators:
        T& operator() (int row, int col);        // Subscript operators in pairs
        const T&  operator() (int row, int col) const;  
        
        T& operator() (int i);        // Subscript operator
        const T&  operator() (int i) const;  
        
        
        Array2D operator = (const Array2D&);
        Array2D operator = (const T a);


        // for linear algebra operations see: arrayops.hpp
        
        // array[][] access:
        bracket_proxy<Array2D, T> operator[](int r){
            return bracket_proxy<Array2D, T>(*this, r);
        }

        Array2D inverse() const;
};


template<class T>
Array2D<T>::~Array2D(){
    delete[] array;
}

template <class T>
Array2D<T>::Array2D(){
    //nrows = 4;
    //ncols = 4;
    //cout << "build \n" << endl;
    //build();
    // cout << "initialize\n" << endl;

    //memset(array, 0, nBytes);
}

//
// initialize as an identity matrix
template <class T>
Array2D<T>::Array2D(bool makeidentidy,  size_t m, size_t n){
    nrows = n;
    ncols = n;
    build();
    if (makeidentidy) {
        assert(m == n);
        for(size_t i = 0; i < nrows; i++) {
            array[i*ncols + i] = 1.0;
        }
    }
}



// new array:
template <class T>
void Array2D<T>::build(){
    storage_size = nrows*ncols;
    nBytes = storage_size * sizeof(T);
    array = new T[storage_size];
    //printf("\n build  \n");

    int i = 0;
    for(i=0; i < storage_size; i++) {
        array[i] = 0.;
    }

}


template <class T>
void Array2D<T>::buildWithParameters(int numrows, int numcols){
    nrows = numrows;
    ncols = numcols;
    //printf("\n build with parameters \n");
    std::cout << " build with parameters  " << std::endl;
    build();
}


// copy constructor:
template <class T>
Array2D<T>::Array2D(const Array2D& other)
    : nrows(other.nrows), ncols(other.ncols){
    storage_size = nrows*ncols;
    nBytes = storage_size * sizeof(T);

    //cout << "copy constructor" << endl;

    array = new T[storage_size];

    //printf("\narray copy constructor\n");
    int i = 0;
    for(i=0; i < storage_size; i++) {
       array[i] = other.array[i];
    }
}


// resize copy constructor:
template <class T>
Array2D<T>::Array2D(const Array2D& other, int nrows, int ncols)
    : nrows(other.nrows+nrows), ncols(other.ncols+ncols){
    storage_size = nrows*ncols;
    nBytes = storage_size * sizeof(T);

    array = new T[storage_size];

    //printf("\narray copy constructor\n");
    int i = 0;
    for(i = 0; i < other.storage_size; i++) {
       array[i] = other.array[i];
    }
}



template <class T>
T& Array2D<T>::operator()(int i, int j) {
	check_indices(i,j);
    return array[i*ncols + j];
}
template <class T>
const T& Array2D<T>::operator()(int i, int j) const {
	check_indices(i,j);
    return array[i*ncols + j];
}


template <class T>
T& Array2D<T>::operator()(int i) {
	check_size(i);
    return array[i];
}
template <class T>
const T& Array2D<T>::operator()(int i) const {
	check_size(i);
    return array[i];
}



template <class T>
Array2D<T> Array2D<T>::operator=(const Array2D& that) {
    /*not using this check  allows delayed allocation 
    of this array type when using it within another class*/
	assert(that.nrows == nrows);
	assert(that.ncols == ncols);

    int i;
    for(i=0; i < storage_size; i++) {
    	array[i] = that.array[i];
    }
    return *this;
}
template <class T>
Array2D<T> Array2D<T>::operator=(const T a) {
    int i;
    for(i=0; i < storage_size; i++) {
    	array[i] = a;
    }
    return *this;
}



template <typename T>
bool Array2D<T>::isSquare() const {
    return getnrows() == getncols();
}




template <typename T>
void print (T const& c){

    int ncol = c.getncols();
    int nrow = c.getnrows();
    for (int i=0; i<nrow; ++i) {
        std::cout << '\n';
        for (int j=0; j<ncol; ++j) {
            std::cout << c(i,j) << ' ';
        }
    }
    std::cout << " \n" << std::endl;
}


/**/  //TODO, use this: 
template <class T>
void Array2D<T>::setonce(T data){
    int i = 0;
    int j = 0;

    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++) {
            array[i*ncols + j] = data;
            data++;
        }
    }
}
/**/





template <class T>
void Array2D<T>::print() {
    int i,j;
    int oldp = cout.precision(numeric_limits<T>::digits10 + 1);
    for( i=0; i<nrows; i++) {
        std::cout << "\n";
        for(j=0; j<ncols; j++) {
            //cout << i << " " << j << " " << &array[i*ncols + j] << endl;
            cout << array[i*ncols + j] << " ";
        }
    }
    cout << "\n" << endl;
    cout.precision(oldp);
}



// template <typename CharT, typename traits, typename T>
// std::basic_ostringstream<CharT,traits>&& operator<<(std::basic_ostringstream<CharT,traits>&& out, const T& t) {
//   static_cast<std::basic_ostream<CharT,traits>&>(out) << t;
//   return std::move(out);
// }

template <class T>
std::ostringstream& Array2D<T>::print(
                        std::ostringstream& out) {
    int i,j;
    //int oldp = cout.precision(numeric_limits<T>::digits10 + 1);
    for( i=0; i<nrows; i++) {
        //out << "\n";
        //out << " " << endl;
        for(j=0; j<ncols; j++) {
            if (j==ncols-1){
                out << array[i*ncols + j] << " \n";
            }
            else {
                out << array[i*ncols + j] << " ";
            }
        }
    }
    //out << "\n";
    //out << " " << endl;
    return out;
}



template <class T>
T sum (T a, T b){
  T result;
  result = a + b;
  return result;
}





#endif //__ARRAY_TEMPLATE_INCLUDED__
