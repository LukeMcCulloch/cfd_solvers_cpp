//*******************************************************************************
// This program generates 2D quad and triangular grids in a rectangular domain.
//
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

//#define CHECKPT {printf("Checkpoint: .s, line .d\n",__FILE__,__LINE__);\
//fflush(stdout);}
#ifdef DEBUG_BUILD
#  define DEBUG(x) fprintf(stderr, x)
#else
#  define DEBUG(x) do {} while (0)
#endif

//=================================
#include <iostream>     // std::cout, std::fixed
#include <fstream>      // write to file
#include <iomanip>    // std::setprecision - only works for output :(
#include <math.h>       // sqrt
//=================================
#include <cstring> //needed for memset
#include <string.h>

//======================================
// my simple array class template (type)
#include "../include/tests_array.hpp"
#include "../include/array_template.hpp"
#include "../include/arrayops.hpp"

//======================================
// my simple vector class template
#include "../include/vector.h"

//======================================
using namespace std;

//======================================
//fwd declarations
struct cell_data;
class gridGen2D;



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
    float xmin, xmax;            //Minimum x and Max x.
    float ymin, ymax;            //Minimum y and Max y

    const float  zero = 0.0;    // minimum x and max x
    const float   one = 1.0;    // minimum y and max y
    int nx;                     // number of nodes in the x-direction
    int ny;                     // number of nodes in the y-direction


    //Output - grid files
    // tria_grid_tecplot.dat
    // quad_grid_tecplot.dat
    // tria.dat
    // quad.dat
    // project.dat


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
gridGen2D::~gridGen2D(){
    printf("destruct");
    delete xs;
    delete ys;
    delete tria;
    delete quad;
    delete x;
    delete y;

}

//
// Default Constructor
gridGen2D::gridGen2D(){

    //  Define the domain: here we define a unit square.
    xmin = zero;
    xmax = one;

    ymin = zero;
    ymax = one;

    //  Define the grid size: the number of nodes in each direction.
    //  NOTE: "ny" is better to be an odd number to place a node at the midpoint
    //        on the left boundary, which will be a corner node in shock diffraction problem.
    nx = 401;
    ny = 401;
    build();
}

// Default Constructor
void gridGen2D::build(){//build

    // structured grid data
    xs = new Array2D<float>(nx,ny);  //Slopes between j and j-1, j and j+1
    ys = new Array2D<float>(nx,ny);  //Slopes between j and j-1, j and j+1

    //--------------------------------------------------------------------------------
    // 1. Generate a structured 2D grid data, (i,j) data: go up in y-direction//
    //
    // j=5 o--------o--------o--------o--------o
    //     |        |        |        |        |
    //     |        |        |        |        |   On the left is an example:
    //     |        |        |        |        |         nx = 5
    // j=4 o--------o--------o--------o--------o         ny = 5
    //     |        |        |        |        |
    //     |        |        |        |        |    +y
    //     |        |        |        |        |    |
    // j=3 o--------o--------o--------o--------o    |
    //     |        |        |        |        |    |____ +x
    //     |        |        |        |        |
    //     |        |        |        |        |
    // j=2 o--------o--------o--------o--------o
    //     |        |        |        |        |
    //     |        |        |        |        |
    //     |        |        |        |        |
    // j=1 o--------o--------o--------o--------o
    //     i=1      i=2      i=3      i=4      i=5

    printf("\nGenerating structured data...\n");

    //  Compute the grid spacing in x-direction
    dx = (xmax-xmin)/float(nx-1);

    // //  Compute the grid spacing in y-direction
    dy = (ymax-ymin)/float(ny-1);

    //  Generate nodes in the domain.


    for (int j=0; j<ny; ++j) {       // Go up in y-direction.
        for (int i=0; i<nx; ++i) {   // Go to the right in x-direction.
            //printf("\ni = %d, j = %d",i,j);
            (*xs)(i,j) = xmin + dx*float(i);
            (*ys)(i,j) = ymin + dy*float(j);
        }
    }

//--------------------------------------------------------------------------------
// 2. Generate unstructured data: 1D array to store the node information.
//
//    - The so-called lexcographic ordering -
//
//   21       22       23       24       25
//     o--------o--------o--------o--------o
//     |        |        |        |        |
//     |        |        |        |        |   On the left is an example:
//   16|      17|      18|      19|      20|           nx = 5
//     o--------o--------o--------o--------o           ny = 5
//     |        |        |        |        |
//     |        |        |        |        |   nnodes = 5x5 = 25
//   11|      12|      13|      14|      15|
//     o--------o--------o--------o--------o
//     |        |        |        |        |
//     |        |        |        |        |
//    6|       7|       8|       9|      10|
//     o--------o--------o--------o--------o
//     |        |        |        |        |
//     |        |        |        |        |
//    1|       2|       3|       4|       5|
//     o--------o--------o--------o--------o
//
   printf("\nGenerating 1D node array for unstructured grid data...\n");

//  Total number of nodes
   nnodes = nx*ny;

//  Allocate the arrays
   x = new Array2D<float>(nnodes,1);
   y = new Array2D<float>(nnodes,1);

// Node data: the nodes are ordered in 1D array.

    for (int j=0; j<ny; ++j) {  //Go up in y-direction.
        for (int i=0; i<nx; ++i) { //Go to the right in x-direction.

            inode = i + (j)*nx;   //<- Node number in the lexcographic ordering
            //(*x)(inode) = (*xs)(inode);
            //(*y)(inode) = (*ys)(inode);
            (*x)(inode) = (*xs)(i,j);
            (*y)(inode) = (*ys)(i,j);
        }
    }
    printf("\n");
    printf( "\n Nodes have been generated:");
    printf( "\n       nx  = %d", nx);
    printf( "\n       ny  = %d", ny);
    printf( "\n    nx*ny  = %d", nx*ny);
    printf( "\n    nnodes = %d", nnodes);
    printf("\n");
    printf( "\n Now, generate elements...");
    printf("\n");




//--------------------------------------------------------------------------------
// 3. Generate unstructured element data:
//
//    We generate both quadrilateral and triangular grids.
//    Both grids are constructed in the unstructured (finite-element) data.
//

// Allocate arrays of triangular and quad connectivity data.

//   Number of quadrilaterals = (nx-1)(ny-1)
    quad = new Array2D<int>( (nx-1)*(ny-1) , 4 );

// //   Number of triangles = 2*(nx-1)*(ny-1)
    tria = new Array2D<int>( 2*(nx-1)*(ny-1) ,  3 );

// (1)Genearte a triangular grid

    printf( "\nGenerating triangular grid...");
    generate_tria_grid();
    printf("\n");
    printf( "\n Number of triangles = %d", ntria);
    printf("\n");
    printf( "\nWriting a tecplot file for the triangular grid...");
    write_tecplot_file("tria_grid_tecplot.dat");//datafile_tria_tec)
    printf( "\n --> File generated:  tria_grid_tecplot.dat");//%d", tria_grid_tecplot);

    printf( "\nWriting a grid file for the triangular grid...");
    write_grid_file("tria.dat");//(datafile_tria)
    printf( "\n --> File generated: tria.dat");//%d", tria);

    printf("\n");

// (2)Generate a quadrilateral grid

    printf( "\nGenerating quad grid...");
    generate_quad_grid();
    printf("\n");
    printf( "\n Number of quads =  %d", nquad);
    printf("\n");
    printf( "\nWriting a tecplot file for the quadrilateral grid...");
    write_tecplot_file("quad_grid_tecplot.dat");//datafile_quad_tec)
    printf( "\n --> File generated:  quad_grid_tecplot.dat");//%d", datafile_quad_tec);

    printf( "\nWriting a grid file for the quadrilateral grid...");
    write_grid_file("quad.dat");//datafile_quad)
    printf( "\n --> File generated:  quad.dat");//%d", datafile_quad);

// (3)Generate a mixed grid. (not implemented. I'll leave it to you// You can do it//)
//

// (4)Write a boundary condition file: to be read by EDU2D-Euler code
    printf( "\nGenerating bcmap file...");
    //open(unit=1, file=datafile_bcmap, status="unknown", iostat=os)
    printf("\nBoundary Segment  Boundary Condition");
    printf("\n               1          freestream");
    printf("\n               2           slip_wall");
    printf("\n               3  outflow_supersonic");
    printf("\n               4  outflow_supersonic");
    printf("\n               5           slip_wall");
    //close(1)

//--------------------------------------------------------------------------------

    printf("\n");
    printf( "\nSuccessfully completed. Stop.");


}


//********************************************************************************
// This subroutine generates triangles by constructing the connectivity data.
//********************************************************************************
void gridGen2D::generate_tria_grid(){
    //Local variables
    int inode, i1, i2, i3, i4;

    // No quads
    nquad = 0;
    //quad = 0;

// Trianguler grid with right-up diagonals (i.e., / ).
//
//  inode+nx   inode+nx+1     i4      i3
//       o--------o           o--------o
//       |     .  |           |     .  |
//       |   .    |     or    |   .    |
//       | .      |           | .      |
//       o--------o           o--------o
//    inode    inode+1        i1      i2
//
// Triangle is defined by the counterclockwise ordering of nodes.
// normals thus "point" out of a closed surface of such triangles

    ntria = -1;

    for (int j=0; j<ny-1; ++j) {
        for (int i=0; i<nx-1; ++i) {

            inode = i + (j)*nx;

//     Define the local numbers (see figure above)
            i1 = inode;
            i2 = inode + 1;
            i3 = inode + nx + 1;
            i4 = inode + nx;

            ntria = ntria + 1;
            (*tria)(ntria,0) = i1;
            (*tria)(ntria,1) = i2;
            (*tria)(ntria,2) = i3;

            ntria = ntria + 1;
            (*tria)(ntria,0) = i1;
            (*tria)(ntria,1) = i3;
            (*tria)(ntria,2) = i4;
        }
    }
    ntria += 1;

 }
//********************************************************************************


//********************************************************************************
// This subroutine generates quads by constructing the connectivity data.
//********************************************************************************
void gridGen2D::generate_quad_grid(){
    //Local variables
    int inode, i1, i2, i3, i4;

    // No triangles
    ntria = 0;
    //tria = 0;

//
//  inode+nx   inode+nx+1     i4      i3
//       o--------o           o--------o
//       |        |           |        |
//       |        |     or    |        |
//       |        |           |        |
//       o--------o           o--------o
//     inode   inode+1        i1      i2
//
// Quad is defined by the counterclockwise ordering of nodes.

    // Quadrilateral grid

    nquad = -1;

    for (int j=0; j<ny-1; ++j) {
        for (int i=0; i<nx-1; ++i) {

            inode = i + (j)*nx;
            //Define the local numbers (see figure above)
            i1 = inode;
            i2 = inode + 1;
            i3 = inode + nx + 1;
            i4 = inode + nx;

            //Order the quad counterclockwise:
            nquad = nquad + 1;
            (*quad)(nquad,0) = i1;
            (*quad)(nquad,1) = i2;
            (*quad)(nquad,2) = i3;
            (*quad)(nquad,3) = i4;
        }
    }
    nquad += 1;
}
//********************************************************************************




//********************************************************************************
//* Write output data file
//*
//* ------------------------------------------------------------------------------
//* Output:  Data file "solution.dat" containing for each cell the following:
//*          cell-center coordinate, density, velocity, pressure, entropy
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
    //void gridGen2D::write_tecplot_file(char* datafile){
    void gridGen2D::write_tecplot_file(const std::string& datafile){

        int os;
        // float entropy;

        ofstream outfile;
        //outfile.open ("tria_grid_tecplot.dat");
        outfile.open (datafile);
        outfile << "title = \"grid\" \n";
        outfile << "variables = \"x\" \"y\" \n";
        outfile << "zone N="  << nnodes << ",E= " << ntria+nquad << ",ET=quadrilateral,F=FEPOINT \n";
        //--------------------------------------------------------------------------------
        for (int i=0; i<nnodes; ++i) {
            outfile << (*x)(i) << '\t' 
                    << (*y)(i) << "\n"; 
        }
        //--------------------------------------------------------------------------------
        //Triangles
        if (ntria > 0) {
            for (int i=0; i<ntria; ++i) {
                outfile << (*tria)(i,0) << '\t' 
                        << (*tria)(i,1) << '\t' 
                        << (*tria)(i,2) << '\t' 
                        << (*tria)(i,2) <<  "\n"; //The last one is a dummy.
            }
        }

        //Quadrilaterals
        if (nquad > 0) {
            for (int i=0; i<nquad; ++i) {
                outfile << (*quad)(i,0) << '\t' 
                        << (*quad)(i,1) << '\t' 
                        << (*quad)(i,2) << '\t' 
                        << (*quad)(i,3) <<  "\n";
            }
        }
        //--------------------------------------------------------------------------------
        outfile.close();



        // ofstream outfile;
        // outfile.open ("tria_grid_tecplot.dat");
        // for (int i=1; i<ncells+1; ++i){
        //     outfile << std::setprecision(16) << cell[i].xc << '\t'
        //             << std::setprecision(16) << cell[i].w(0) << '\t'
        //             << std::setprecision(16) << cell[i].w(1) << '\t'
        //             << std::setprecision(16) << cell[i].w(2) << '\t'
        //             << std::setprecision(16) << entropy <<  "\n";
        // }
        // outfile.close();

}
//--------------------------------------------------------------------------------




//********************************************************************************
// This subroutine writes a grid file to be read by a solver.
// NOTE: Unlike the tecplot file, this files contains boundary info.
//********************************************************************************
void gridGen2D::write_grid_file(const std::string& datafile) {//(char* datafile)
    int i,j,os;
//--------------------------------------------------------------------------------
    ofstream outfile;
    outfile.open (datafile);

//--------------------------------------------------------------------------------
// Grid size: # of nodes, # of triangles, # of quadrilaterals
    outfile <<  nnodes << ntria << nquad << "\n";

//--------------------------------------------------------------------------------
// Node data
    for (int i=0; i<nnodes; ++i) {
        outfile <<  (*x)(i) <<(*y)(i) << "\n";
;    }

//--------------------------------------------------------------------------------
// Triangle connectivity
    if (ntria > 0) {
        for (int i=0; i<ntria; ++i) {
            outfile <<  (*tria)(i,0) << (*tria)(i,1) << (*tria)(i,2) <<  "\n";
        }
    }

// Quad connectivity
    if (nquad > 0) {
        for (int i=0; i<nquad; ++i) {
            outfile <<  (*quad)(i,0) << (*quad)(i,1) << (*quad)(i,2) << (*quad)(i,3) <<  "\n";
        }
    }

// Boundary data:
// NOTE: These boundary data are specific to the shock diffraction problem.
//
//  Example: nx=ny=7
//
//   in = inflow
//    w = wall
//    e = outflow
//    o = interior nodes
//  inw = this node belongs to both inflow and wall boundaries.
//   we = this node belongs to both wall and outflow boundaries.
//
//   inw----w----w----w----w----w----we
//     |    |    |    |    |    |    |
//    in----o----o----o----o----o----e
//     |    |    |    |    |    |    |
//    in----o----o----o----o----o----e
//     |    |    |    |    |    |    |
//   inw----o----o----o----o----o----e
//     |    |    |    |    |    |    |
//     w----o----o----o----o----o----e
//     |    |    |    |    |    |    |
//     w----o----o----o----o----o----e
//     |    |    |    |    |    |    |
//    we----e----e----e----e----e----e
//

// Number of boundary segments
    outfile <<  5 << "\n";

    outfile <<  (ny-1)/2+1  << "\n"; //Inflow
    outfile <<  (ny-1)/2+1  << "\n"; //Left Wall
    outfile <<   nx << "\n";         //Bottom Outflow
    outfile <<   ny  << "\n";        //Right  Outflow
    outfile <<   nx  << "\n";        //Top Wall

    outfile << "\n";

// Inflow boundary
    //do j = ny, (ny-1)/2+1, -1
    for (int j=ny-1; j<(ny-1)/2+1; ++j) {
        i = 1;
        outfile <<  i + (j-1)*nx  << "\n";  
    }

// // Left wall boundary
    //do j = (ny-1)/2+1, 1, -1
    for (int j=(ny-1)/2+1; j>=1; --j) {
        i = 1;
        outfile <<  i + (j-1)*nx  << "\n";  
    }

// // Bottom outflow boundary
//     do i = 1, nx
    for (int i=nx-1; i<nx; ++i) {
        j = 1;
        outfile <<  i + (j-1)*nx  << "\n";  
    }

// // Right outflow boundary
//     do j = 1, ny
//     i = nx
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

// // Top wall boundary
//     do i = nx, 1, -1
//     j = ny
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

//--------------------------------------------------------------------------------
    outfile.close();
}
//********************************************************************************



void driverGrid2D(){
    gridGen2D Grid;
    printf("\nGridding Complete\n");
    return;
}