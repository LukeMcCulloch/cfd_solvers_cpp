//=================================
// include guard
#ifndef __gridGen2D_INCLUDED__
#define __gridGen2D_INCLUDED__



namespace Grid2D{



// use an array of structs (may be inefficient//)
struct cell_data{
    float xc;  // Cell-center coordinate
    Array2D<float> u  = Array2D<float>(3,1);  // Conservative variables = [rho, rho*u, rho*E]
    Array2D<float> u0 = Array2D<float>(3,1);  // Conservative variables at the previous time step
    Array2D<float> w  = Array2D<float>(3,1);  // Primitive variables = [rho, u, p]
    Array2D<float> dw = Array2D<float>(3,1);  // Slope (difference) of primitive variables
    Array2D<float> res= Array2D<float>(3,1);  // Residual = f_{j+1/2) - f_{j-1/2)
};

class gridGen2D{

public:


    //constructor
    gridGen2D();
    //gridGen2D(int nxi, int nyi);
    // explicit constructor declaring size nrow,ncol:
    explicit gridGen2D(int nxi, int nyi):
                    nx(nxi), ny(nyi){
        //  Define the domain: here we define a unit square.
        xmin = zero;
        xmax = one;

        ymin = zero;
        ymax = one;
        build();
    }
    //destructor
    ~gridGen2D();

    void build();
    void generate_tria_grid();
    void generate_quad_grid();

    //output
    void write_tecplot_file(const std::string& datafile);
    void write_grid_file(const std::string& datafile);


    //Input  - domain size and grid dimensions
    float xmin, xmax;           //Minimum x and Max x.
    float ymin, ymax;           //Minimum y and Max y

    const float  zero = 0.0;    // minimum x and max x
    const float   one = 1.0;    // minimum y and max y
    int nx;                     // number of nodes in the x-direction
    int ny;                     // number of nodes in the y-direction


    //Output - grid files
    std::string  datafile_tria_tec = "tria_grid_tecplot.dat";
    std::string  datafile_quad_tec = "quad_grid_tecplot.dat";
    std::string  datafile_tria = "tria.grid";
    std::string  datafile_quad = "quad.grid";
    std::string  datafile_bcmap = "project.dat";


    // structured grid data
    Array2D<float>*  xs;  //Slopes between j and j-1, j and j+1
    Array2D<float>*  ys;  //Slopes between j and j-1, j and j+1


    //Local variables
    int nnodes; //Total number of nodes
    int  ntria; //Total number of triangles
    int  nquad; //Total number of quadrilaterals
    int  inode; //Local variables used in the 1D nodal array


    Array2D<int>* tria;      //Triangle connectivity data
    Array2D<int>* quad;      //Quad connectivity data
    Array2D<float>*  x;   //Nodal x coordinates, 1D array
    Array2D<float>*  y;   //Nodal y coordinates, 1D array

    float dx; //Uniform grid spacing in x-direction = (xmax-xmin)/nx
    float dy; //Uniform grid spacing in y-direction = (ymax-ymin)/ny
    int i, j, os;
};
//

//=================================
// the actual function
void driverGrid2D();

}

#endif 