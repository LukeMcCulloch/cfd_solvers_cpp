//=================================
// include guard
#ifndef __eulershock1d_INCLUDED__
#define __eulershock1d_INCLUDED__


//======================================
// my simple array class template (type)
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"


namespace EulerSolver1D
{




// use an array of structs (may be inefficient//)
struct cell_data{
    float xc;  // Cell-center coordinate
    Array2D<float> u  = Array2D<float>(3,1);  // Conservative variables = [rho, rho*u, rho*E]
    Array2D<float> u0 = Array2D<float>(3,1);  // Conservative variables at the previous time step
    Array2D<float> w  = Array2D<float>(3,1);  // Primitive variables = [rho, u, p]
    Array2D<float> dw = Array2D<float>(3,1);  // Slope (difference) of primitive variables
    Array2D<float> res= Array2D<float>(3,1);  // Residual = f_{j+1/2) - f_{j-1/2)
};



class Solver{

public:

    //constructor
    Solver();
    // destructor
    ~Solver();


    void Euler1D();
    void initialize( int ncells, 
                float dx, float xmin, const float gamma);
    float timestep(float cfl, float dx, float gamma, int ncells);
    
    // limiter:
    float minmod(float a, float b);
    
    // transforms:
    void w2u_efficient( Array2D<float>& w, Array2D<float>& u );
    void u2w_efficient( Array2D<float>& u, Array2D<float>& w );
    Array2D<float> u2w( Array2D<float>& u);
    Array2D<float> w2u( Array2D<float>& w);
    
    //flux:
    Array2D<float> roe_flux(Array2D<float>&  wL, Array2D<float>&  wR);
    Array2D<float> euler_physical_flux(Array2D<float>& w);

    //print;
    void output();

    struct constants{
        const float  zero = 0.0;
        const float   one = 1.0;
        const float  half = 0.5;
        const float gamma = 1.4;  //Ratio of specific heats for air
    };

    //Numeric parameters: [Note: no Fortran-like way to handle precision?]
    //const int p2 = 10;
    const float  zero = 0.0;
    const float   one = 1.0;
    const float  half = 0.5;
    const float gamma = 1.4;  //Ratio of specific heats for air

    float xmin, xmax; //Left and right ends of the domain
    float dx;         //Cell spacing (uniform grid)
    float t, tf;      //Current time and final time
    float cfl, dt;    //CFL number and global time step
    int   ncells;     //Total number of cells
    int   nsteps;     //Number of time steps
    int   itime;      //Index for time stepping
    int   istage;     //Index for Runge-Kutta stages
    int   i, j;

    //Local variables used for computing numerical fluxes.
    // init arrays here
    Array2D<float>  dwl = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  dwr = Array2D<float>(3,1);  //Slopes between j and j-1, j and j+1
    Array2D<float>  wL = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  wR = Array2D<float>(3,1);   //Extrapolated states at a face
    Array2D<float>  flux = Array2D<float>(3,1); //Numerical flux

    cell_data* cell;

};


//=================================
// the driver function
void driverEuler1D();

} //namespace Euler1D



#endif 