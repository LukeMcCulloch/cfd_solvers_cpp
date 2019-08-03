//*******************************************************************************
// This program generates 2D quad and triangular grids in a rectangular domain.
//
// Note that this is separate conceputally and practically from the 
// Euler 2D code to read in a grid.
// (and we want to keep it that way!)
//
//      Boundary information is set up for a shock-diffraction problem
//
//                                 Wall
//                         --------------------
//     Post-shock (inflow) |                  |
//                         |->Shock           |            o: Corner node
//                         |  Mach=5.09       |
//                  .......o                  |Outflow
//                    Wall |                  |
//                         |                  |
//                         |                  |
//                         --------------------
//                               Outflow
//
// 3-Step Generation:
//
// 1. Generate a temporary structured grid data for nodes: xs(i,j) and ys(i,j)
// 2. Generate a 1D node array: x(1:nnodes), y(1:nnodes)
// 3. Generate element connectivity data: tria(1:ntria,3), quad(1:nquad,4)
//
//
// Inuput:
//        xmin, xmax = x-coordinates of the left and right ends
//        ymin, ymax = y-coordinates of the bottom and top ends
//                nx = number of nodes in x-direction
//                ny = number of nodes in y-direction
//
//        (NOTE: All input parameters are defined inside the program.)
//
// Output:
//        tria_grid_tecplot.dat = tecplot file of the triangular grid
//        quad_grid_tecplot.dat = tecplot file of the quadrilateral grid
//                    tria.grid = triangular-grid file for EDU2D-Euler solver
//                    quad.grid = quadrilateral-grid file for EDU2D-Euler solver
//                project.bcmap = file that contains boundary condition info
//
//        (NOTE: The grid file format is specific to the EDU2D-Euler solver.)
//
//        (NOTE: The grid file contains 5 boundary segments which is designed
//               for the shock-diffraction problem. See write_grid_file() on the
//               bottom of this file. Boudnary condition files project.bcmap
//               specify the boundary condition to be applied to each segment.)
//
//
//
//        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//
// the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//
// This is Version 0 (December 2011).
// This F90 code is written and made available for an educational purpose as well
// as for generating grids for the EDU2D-Euler code.
// This file may be updated in future.
//
// Katate Masatsuka, January 2012. http://www.cfdbooks.com
// translated by Luke McCulloch June 2019 LukeMcCulloch@github.com
//*******************************************************************************

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