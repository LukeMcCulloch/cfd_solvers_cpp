//********************************************************************************
//* Educationally-Designed Unstructured 2D (EDU2D) Code
//*
//*
//*
//*     This module belongs to the inviscid version: EDU2D-Euler-RK2
//*
//*
//*
//* This file contains 4 modules:
//*
//*  1. module EulerSolver2D      - Some numerical values, e.g., zero, one, pi, etc.
//*  2. module EulerSolver2D_type - Grid data types: node, edge, face, element, etc.
//*  3. module EulerSolver2D   - Parameters and arrays mainly used in a solver.
//*  4. module EulerSolver2D  - Subroutines for dynamic allocation
//*  5. module EulerSolver2D      - Subroutines for reading/constructing/checking grid data
//*
//* All data in the modules can be accessed by the use statement, e.g., 'use constants'.
//*
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to C++ by Luke McCulloch, PhD.
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015).
//* This F90 code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************

//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//* 1. module EulerSolver2D
//*
//* Some useful constants are defined here.
//* They can be accessed by the use statement, 'use constants'.
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*        translated to C++ by Luke McCulloch, PhD.
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015). -> C++ June 2019
//* In the spirit of its fortran parent,
//* This C++ code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************

//======================================
// 2D Euler approximate Riemann sovler
#include "../include/EulerUnsteady2D_basic_package.h"
//======================================
// i/o
#include <iostream>     // cout, std::fixed
#include <fstream>      // write to file
//======================================
// line based parsing, including streams
#include <sstream>
#include <string>

//======================================
// mesh-y math-y functions
#include "../include/MathGeometry.h" 


using std::cout;
using std::endl;

// constructors and destrutors that do nothing
//EulerSolver2D::MainData2D::MainData2D() {}
//EulerSolver2D::MainData2D::~MainData2D() {
//   if ()
//}



//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//* 5. module EulerSolver2D
//*
//* This module contians subroutines used for reading a grid, constructing
//* additional grid data, and check the grid data.
//*
//*  - my_alloc_int_ptr       : Allocate/reallocate an integer 1D array
//*  - my_alloc_p2_ptr        : Allcoate/reallocate a real 1D array
//*  - my_alloc_p2_matrix_ptr : Allcoate/reallocate a real 2D array
//*
//*
//* Public subroutines:
//*
//*  - read_grid            : Read a grid file, allocate necessary data.
//*  - construct_grid_data  : Construct additional data, allocate more data.
//*  - check_grid_data      : Check the whole grid data.
//*
//* Private functions and subroutines:
//*
//*  - tri_area             : Computed a triangle area
//*  - check_skewness_nc    : Compute the skewness (e*n).
//*  - compute_ar           : Compute aspect ratio at node and element.
//*
//*
//*
//*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
//*
//* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
//*
//* This is Version 0 (July 2015).
//* This F90 code is written and made available for an educational purpose.
//* This file may be updated in future.
//*
//********************************************************************************
//namespace EulerSolver2D{



//********************************************************************************
//* Read the grid and the exact solution.
//* ------------------------------------------------------------------------------
//*  Input: datafile_grid_in  = filename of the grid file
//*         datafile_bcmap_in = filename of the bc file
//*
//* Output: nnodes, ncells, node(:), elm(:), bound(:) = data used in the solver
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
//* 1. "datafile_grid_in" is assumed to have been written in the following format:
//*
//*   -----------------------------------------------------------------------
//*    write(*,*) nnodes, ntria, nquad //Numbers of nodes, triangles and quads
//*
//*   do i = 1, nnodes
//*    write(*,*) x(i), y(i) //(x,y) coordinates of each node
//*   end do
//*
//*   do i = 1, ntria        //Nodes of triangles ordered counterclockwise
//*    write(*,*) node_1(i), node_2(i), node_3(i)
//*   end do
//*
//*   do i = 1, nquad        //Nodes of quadrilaterals ordered counterclockwise
//*    write(*,*) node_1(i), node_2(i), node_3(i), node_4(i)
//*   end do
//* 
//*    write(*,*) nbound     //Number of boundary segments
//*
//*   do i = 1, nbound
//*    write(*,*) nbnodes(i) //Number of nodes on each segment
//*   end do
//*
//*   do i = 1, nbound
//*    do j = 1, nbnodes(i)
//*     write(*,*) bnode(j)  //Node number of each node j in segment i
//*    end do
//*   end do
//*   -----------------------------------------------------------------------
//*
//*   NOTE: Add the first node to the end if the segment is closed
//*         (e.g., airfoil) The number of nodes will be the actual number + 1
//*         in that case.
//*
//*   NOTE: Boundary nodes must be ordered such that the domain is on the left.
//*
//********************************************************************************
//*
//* 2. "datafile_bcmap_in" is assumed have been written in the following format:
//*
//*   -----------------------------------------------------------------------
//*    write(*,*) "Boundary Segment              Boundary Condition"
//*   do i = 1, nbound
//*    write(*,*) i, bc_name
//*   end do
//*   -----------------------------------------------------------------------
//*
//*   NOTE: bc_name is the name of the boundary condition, e.g.,
//*
//*         1. "freestream"
//*             Roe flux with freestream condition on the right state.
//*
//*         2. "slip_wall"
//*             Solid wall condition. Mass flux through the boundary is set zero.
//*
//*         3. "outflow_supersonic"
//*             Just compute the boundary flux by the physical Euler flux
//*             (equivalent to the interior-extrapolation condition.)
//*
//*         4. "outflow_back_pressure"
//*             Fix the back pressure. This should work for subsonic flows in a
//*             large enough domain.
//*
//*         Something like the above needs to be implemented in a solver.
//*
//********************************************************************************
//* Data to be read and stored:
//*
//* 1. Some numbers
//*    nnodes        = Number of nodes
//*    ntria         = Number of triangular elements
//*    nquad         = Number of quadrilateral elements
//*    nelms         = Total number of elements (=ntria+nquad)
//*
//* 2. Element data:
//*    elm(1:nelms).nvtx   =  Number of vertices of each element
//*    elm(1:nelms).vtx(:) = Pointer to vertices of each element
//*
//* 3. Node data: nodes are stored in a 1D array
//*    node(1:nnodes).x     = x-coordinate of the nodes
//*    node(1:nnodes).y     = y-coordinate of the nodes
//*
//* 4. Boundary Data:
//*    nbound                   = Number of boundary segments
//*    bound(1:nbound).nbnodes  = Number of nodes in each segment
//*    bound(1:nbound).bnode(:) = List of node numbers for each segment
//*    bound(1:nbound).bc_type  = Boundary condition name for each segment
//*    bound(1:nbound).bc_type  = Boundary condition name for each segment
//*
//********************************************************************************


void EulerSolver2D::MainData2D::read_grid(std::string datafile_grid_in, 
                                               std::string datafile_bcmap_in)
{
   /*
   *   This is inside out from the C++ perspective.
   *   'Somebody' needs to read the grid and figure 
   *   out how many nodes we have.
   *   Then create a MainData2D instance and allocate node arrays
   *   within the constructor.
   */
   //use EulerSolver2D, only : nnodes, node, ntria, nquad, nelms, elm, nbound, bound

   //Local variables
   int i, j, os, dummy_int;

   //--------------------------------------------------------------------------------
   // 1. Read grid file>: datafile_grid_in

   cout << "Reading the grid file...." << datafile_grid_in << endl;

   //  Open the input file.
   std::ifstream infile;
   infile.open(datafile_grid_in);

   // line based parsing, using string streams:
   std::string line;

   // READ: Get the size of the grid.
   std::getline(infile, line);
   std::istringstream iss(line);
   iss >> nnodes >> ntria >> nquad;
   cout <<  "Found... " <<
            " nnodes =  " << nnodes << 
            "   ntria = " <<  ntria <<  
            "   nquad = " <<  nquad << endl;
   nelms = ntria + nquad;

   // //  Allocate node and element arrays.
   std::cout << "Allocating node_type" << std::endl;
   std::cout << "    for " << nnodes << " nodes " << std::endl;
   node = new node_type[nnodes];


   std::cout << " Allocating elm_type" << std::endl;
   std::cout << "    for " << nelms << " elements " << std::endl;
   elm = new elm_type[nelms];


   // // READ: Read the nodal coordinates
   cout << "reading nodal coords" << endl;
   for (size_t i = 0; i < nnodes; i++) {
      std::getline(infile, line);
      std::istringstream iss(line);
      iss >> node[i].x >> node[i].y ;
      // could declare here:
      // node[i].u     = new Array2D<real>(nq,1);
      // node[i].du    = new Array2D<real>(nq,1);
      // node[i].w     = new Array2D<real>(nq,1);
      // node[i].gradw = new Array2D<real>(nq,2); //<- 2: x and y components.
      // node[i].res   = new Array2D<real>(nq,1);

      //node[i].nelms = 0;
      //std::cout << "i = " << i << " of " << nnodes-1 << std::endl;

   }
   cout << "done reading nodal coords" << endl;

      
   // Read element-connectivity information

   // Triangles: assumed that the vertices are ordered counterclockwise
   //
   //         v3
   //         /\
   //        /  \
   //       /    \
   //      /      \
   //     /        \
   //    /__________\
   //   v1           v2

   // READ: read connectivity info for triangles
   if (ntria > 0) {
      for (size_t i = 0; i < ntria; i++) {
         std::getline(infile, line);
         std::istringstream in(line);
         elm[i].nvtx = 3;
         elm[i].vtx = new Array2D<int>(3,1) ;

         std::string type;
         in >> type;                  //and read the first whitespace-separated token


         int x, y, z;
         in >> x >> y >> z;       //now read the whitespace-separated floats
         (*elm[i].vtx)(0,0) = x;
         (*elm[i].vtx)(1,0) = y;
         (*elm[i].vtx)(2,0) = z;
      }
   }
   // // Quads: assumed that the vertices are ordered counterclockwise
   // //
   // //        v4________v3
   // //         /        |
   // //        /         |
   // //       /          |
   // //      /           |
   // //     /            |
   // //    /_____________|
   // //   v1             v2

   // // READ: read connectivity info for quadrilaterals
   if (nquad > 0) {
      for (size_t i = 0; i < nquad; i++) {
         std::getline(infile, line);
         std::istringstream in(line);
         elm[ntria+i].nvtx = 4;
         elm[ntria+i].vtx = new Array2D<int>(4,1);

         int x1,x2,x3,x4;
         in >> x1 >> x2 >> x3 >> x4;       //now read the whitespace-separated ...ints
         (*elm[ntria+i].vtx)(0,0) = x1;
         (*elm[ntria+i].vtx)(1,0) = x2;
         (*elm[ntria+i].vtx)(2,0) = x3;
         (*elm[ntria+i].vtx)(3,0) = x4;
      }
   }

   //  Write out the grid data.
   cout << " Total numbers:" << endl;
   cout << "       nodes = " << nnodes << endl;
   cout << "   triangles = " << ntria << endl;
   cout << "       nquad = " << nquad << endl;
   cout << "       nelms = " << nelms << endl;
   

   // Read the boundary grid data

   // READ: Number of boundary condition types
   std::getline(infile, line);
   std::istringstream in(line);
   in >> nbound;
   bound = new bgrid_type[nbound];


   // // READ: Number of Boundary nodes (including the starting one at the end if
   // // it is closed such as an airfoil.)
   for (size_t i = 0; i < nbound; i++) {
      std::getline(infile, line);
      std::istringstream in(line);
      in >> bound[i].nbnodes;
      bound[i].bnode = new Array2D<int>(bound[i].nbnodes,1);
   }
   cout << "       nbnodes = " << bound[i].nbnodes << endl;

   // // READ: Read boundary nodes
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {
         std::getline(infile, line);
         std::istringstream in(line);
         int init;
         in >> init;
         (*bound[i].bnode)[j][0] = init;
         //(*bound[i].bnode)(j,0) = init;

         // testing array access:
         //int some = (*bound[i].bnode)[j][0];
         //cout << "get some " << some << endl;
      }
   }

   //  Print the boundary grid data.
   std::cout << " Boundary nodes:" << std::endl;
   std::cout << "    segments = " << nbound << std::endl;
      for (size_t i = 0; i < 2; i++) {
         std::cout <<  " boundary = " << i << 
                     "   bnodes = " <<  bound[i].nbnodes <<  
                     "   bfaces = " <<  bound[i].nbnodes << std::endl;
      }
      for (size_t i = nbound-2; i < nbound; i++) {
         std::cout <<  " boundary = " << i << 
                     "   bnodes = " <<  bound[i].nbnodes <<  
                     "   bfaces = " <<  bound[i].nbnodes << std::endl;
      }
   

   infile.close(); // close datafile_grid_in
   // End of Read grid file>: datafile_grid_in
   //--------------------------------------------------------------------------------

   //--------------------------------------------------------------------------------
   // 2. Read the boundary condition data file

   std::cout << "" << std::endl;
   std::cout << "Reading the boundary condition file...." << datafile_bcmap_in << std::endl;
   std::cout << "" << std::endl;

   // // Open the input file.
   //std::ofstream outfile;
   std::ifstream outfile;
   outfile.open (datafile_bcmap_in);

   std::getline(outfile, line);

   // READ: Read the boundary condition type
   for (size_t i = 0; i < nbound; i++) {
      std::getline(outfile, line);
      std::istringstream in(line);
      in >> dummy_int, bound[i].bc_type;
   }

   //  Print the data
   std::cout << " Boundary conditions:" << std::endl;
   for (size_t i = 0; i < 2; i++) {
      std::cout << " boundary" << i << "  bc_type = " << bound[i].bc_type << std::endl;
   }
   for (size_t i = nbound-2; i < nbound; i++) {
      std::cout << " boundary" << i << "  bc_type = " << bound[i].bc_type << std::endl;
   }

   std::cout << "" << std::endl;

   // close(2)
   outfile.close(); // close datafile_bcmap_in

   // End of Read the boundary condition data file
   //--------------------------------------------------------------------------------
   return;

 } // end function read_grid



//********************************************************************************
//* Construct the grid data:
//*
//* The following data, needed for NCFV method, will be constructed based on the
//* data read from the grid file.
//*
//* 1. Element data:
//*    elm(:).nnghbrs  = Number of element neighbors of each element
//*    elm(:).nghbr(:) = List of element neighbors of each element
//*    elm(:).x        = x-coordinate of the centroid
//*    elm(:).y        = y-coordinate of the centroid
//*    elm(:).vol      = Volume of the element
//*
//*
//* 2. Node data:
//*    node(:).nnghbrs = Number of node neighbors of each node
//*    node(:).nghbr(:)= List of node neighbors of each node
//*    node(:).nelms   = Number of adjacent elements of each node
//*    node(:).elm     = List of adjacent elements of each node
//*    node(:).vol     = Volume of the dual volume around each node
//*
//* 3. Edge data:
//*    edge(:).n1, n2  = End nodes of each edge (edge points n1 -> n2)
//*    edge(:).e1, e2  = Left and right elements of each edge
//*    edge(:).dav     = Unit directed area vector of each edge
//*    edge(:).da      = Magnitude of the directed area vector for each edge
//*    edge(:).ev      = Unit edge vector of each edge (vector n1 -> n2)
//*    edge(:).e       = Magnitude of the edge vector for each edge
//*
//*
//* 4. Boudnary data
//*    bound(:).bnx    = Outward normal at boundary nodes (x-component of unit vector)
//*    bound(:).bny    = Outward normal at boundary nodes (y-component of unit vector)
//*    bound(:).bn     = Magnitude of (bnx,bny)
//*    NOTE: In this code, the above normal vector at boundary nodes is computed by
//*          a quadratic fit. It is sufficiently accuarte for 3rd-order schemes.
//*          See http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2015v281pp518-555_preprint.pdf
//*          for details on the quadratic approximation for computing more accurate normals.
//*    bound(:).bfnx   = Outward normal at boundary nodes (x-component of unit vector)
//*    bound(:).bfny   = Outward normal at boundary nodes (y-component of unit vector)
//*    bound(:).bfn    = Magnitude of (bfnx,bfny)
//*    bound(:).belm   = Element to which the boundary face belongs
//*
//********************************************************************************

void EulerSolver2D::MainData2D::construct_grid_data(){

//  use EulerSolver2D , only : nnodes, node, nelms, elm, nedges, edge, nbound, bound, face, nfaces
//  use EulerSolver2D    , only : p2, zero, half, third
//  use EulerSolver2D, only : my_alloc_int_ptr, my_alloc_p2_ptr, my_alloc_p2_matrix_ptr

// //Local variables
int i, j, k, ii, in, im, jelm, v1, v2, v3, v4;
real x1, x2, x3, x4, y1, y2, y3, y4, xm, ym, xc, yc;
real xj, yj, xm1, ym1, xm2, ym2, dsL,dsR,dx,dy;
bool found;
int vL, vR, n1, n2, e1, e2;
int vt1, vt2, ielm;
int ave_nghbr, min_nghbr, max_nghbr, imin, imax;
int iedge;

//  real(p2)                          :: ds
real ds;

// // Some initialization
 v2 = 0;
 vL = 0;
 im = 0;
 jelm = 0;

//   write(*,*) "Constructing grid data...."
cout << "construct grid data...." << endl;

// // Initializations
for (size_t i = 0; i < nnodes; i++) {
   node[i].nelms = 0;
} 
nedges = 0;

// //--------------------------------------------------------------------------------
// // Loop over elements and construct the fololowing data.
// //
// // 1. Surrounding elements: node(:).nelms, node(:).elm(:)
// //
// //    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
// //             node[i].nelms = 4
// //             node[i].elm(1) = 23
// //             node[i].elm(2) = 13
// //             node[i].elm(3) = 41
// //             node[i].elm(4) = 101
// //
// //        o-------o-------------o
// //       /        |   .         |
// //      /    23   |      41     |
// //     o----------o-------------o
// //      \        i \            |
// //       \   101    \     13    |
// //        \          \          | 
// //         o----------o---------o
// //
// // 2. Element quantities  : elm(:).x,elm(:).y,elm(:).vol
// //
// //  o-----------o            
// //   \          |            o
// //    \    (x,y)|           / \
// //     \   .    |          /   \
// //      \       |         /  .  \    (x,y): centroid coordinates
// //       \      |        / (x,y) \     vol: volume of element
// //        o-----o       o---------o

//   elements : do i = 1, nelms


for ( int i = 0; i < nelms; ++i ) {

   v1 = (*elm[i].vtx)(0,0);
   v2 = (*elm[i].vtx)(0,0);
   v3 = (*elm[i].vtx)(0,0);

   x1 = node[v1].x;
   x2 = node[v2].x;
   x3 = node[v3].x;

   y1 = node[v1].y;
   y2 = node[v2].y;
   y3 = node[v3].y;

// // Distribute the element index to nodes.

   /*
   * DESIGN CHOICE HERE -- use std<vector?> or similar?
   * I don't think we need a true dynamic (i.e. resizable) array
   */
   node[v1].nelms = node[v1].nelms + 1;
   node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
   node[v1].elm[node[v1].nelms] = i;

   node[v2].nelms = node[v2].nelms + 1;
   node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
   node[v2].elm[node[v2].nelms] = i;

   node[v3].nelms = node[v3].nelms + 1;
   node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
   node[v3].elm[node[v3].nelms] = i;

// // Compute the cell center and cell volume.
   //tri_or_quad : if (elm(i).nvtx==3) then
   if (elm[i].nvtx==3) {

//   Triangle centroid and volume
    elm[i].x   = third*(x1+x2+x3);
    elm[i].y   = third*(y1+y2+y3);
    elm[i].vol = tri_area(x1,x2,x3,y1,y2,y3);
   }
   else if (elm[i].nvtx==4) {

//   this is a quad. Get the 4th vertex.
      v4 = (*elm[i].vtx)(4,0);
      x4 = node[v4].x;
      y4 = node[v4].y;
//   Centroid: median dual
//   (Note: There is an alternative. See Appendix B in Nishikawa AIAA2010-5093.)
      xm1 = half*(x1+x2);
      ym1 = half*(y1+y2);
      xm2 = half*(x3+x4);
      ym2 = half*(y3+y4);
      elm[i].x   = half*(xm1+xm2);
      elm[i].y   = half*(ym1+ym2);
//   Volume is computed as a sum of two triangles: 1-2-3 and 1-3-4.
      elm[i].vol = tri_area(x1,x2,x3,y1,y2,y3) + \
                  tri_area(x1,x3,x4,y1,y3,y4);

      xc = elm[i].x;
      yc = elm[i].y;
      if (tri_area(x1,x2,xc,y1,y2,yc)<zero) {
         cout << " Centroid outside the quad element 12c: i=" << i << endl;
         cout << "  (x1,y1)=" << x1 << y1  << endl;
         cout << "  (x2,y2)=" << x2 << y2  << endl;
         cout << "  (x3,y3)=" << x3 << y3  << endl;
         cout << "  (x4,y4)=" << x4 << y4  << endl;
         cout << "  (xc,yc)=" << xc << yc  << endl;
         //stop
      }

      if (tri_area(x2,x3,xc,y2,y3,yc)<zero) {
         cout << " Centroid outside the quad element 23c: i=" << i << endl;
         cout << "  (x1,y1)=" << x1 << y1  << endl;
         cout << "  (x2,y2)=" << x2 << y2  << endl;
         cout << "  (x3,y3)=" << x3 << y3  << endl;
         cout << "  (x4,y4)=" << x4 << y4  << endl;
         cout << "  (xc,yc)=" << xc << yc  << endl;
         //stop
      }

      if (tri_area(x3,x4,xc,y3,y4,yc)<zero) {
         cout << " Centroid outside the quad element 34c: i=" << i << endl;
         cout << "  (x1,y1)=" << x1 << y1 << endl;
         cout << "  (x2,y2)=" << x2 << y2 << endl;
         cout << "  (x3,y3)=" << x3 << y3 << endl;
         cout << "  (x4,y4)=" << x4 << y4 << endl;
         cout << "  (xc,yc)=" << xc << yc << endl;
         //stop
      }

      if (tri_area(x4,x1,xc,y4,y1,yc)<zero) {
         cout << " Centroid outside the quad element 41c: i=" << i << endl;
         cout << "  (x1,y1)=" << x1 << y1  << endl;
         cout << "  (x2,y2)=" << x2 << y2  << endl;
         cout << "  (x3,y3)=" << x3 << y3  << endl;
         cout << "  (x4,y4)=" << x4 << y4  << endl;
         cout << "  (xc,yc)=" << xc << yc  << endl;
         //stop
      }

//  Distribution of element number to the 4th node of the quadrilateral
   node[v4].nelms = node[v4].nelms + 1;
   node[v4].elm = new Array2D<int>(node[v4].nelms, 0);
   (*node[v4].elm)[node[v4].nelms][0] = i;

   }//    endif tri_or_quad

}//   end do elements (i loop)

// // Median dual volume

for (int i = 0; i < nnodes; i++) {
   node[i].vol = zero;
}

//   elementsv : do i = 1, nelms
for ( int i = 0; i < nelms; ++i ) {
   

   v1 = (*elm[i].vtx)(0,0);
   v2 = (*elm[i].vtx)(1,0);
   v3 = (*elm[i].vtx)(2,0);

//    tri_or_quadv : 
   if (elm[i].nvtx==3) {
//   Dual volume is exactly 1/3 of the volume of the triangle.
      node[v1].vol = node[v1].vol + third*elm[i].vol;
      node[v2].vol = node[v2].vol + third*elm[i].vol;
      node[v3].vol = node[v3].vol + third*elm[i].vol;

   }  else if (elm[i].nvtx==4) {
         v4 = (*elm[i].vtx)(4,0);

         x1 = node[v1].x;
         x2 = node[v2].x;
         x3 = node[v3].x;
         x4 = node[v4].x;
         xc = elm[i].x;

         y1 = node[v1].y;
         y2 = node[v2].y;
         y3 = node[v3].y;
         y4 = node[v4].y;
         yc = elm[i].y;

// - Vertex 1
         xj = node[v1].x;
         yj = node[v1].y;
         xm1 = half*(xj+x2);
         ym1 = half*(yj+y2);
         xm2 = half*(xj+x4);
         ym2 = half*(yj+y4);

//   Median volume is computed as a sum of two triangles.
         node[v1].vol = node[v1].vol + \
                        tri_area(xj,xm1,xc,yj,ym1,yc) + \
                        tri_area(xj,xc,xm2,yj,yc,ym2);

// - Vertex 2
         xj = node[v2].x;
         yj = node[v2].y;
         xm1 = half*(xj+x3);
         ym1 = half*(yj+y3);
         xm2 = half*(xj+x1);
         ym2 = half*(yj+y1);

//   Median volume is computed as a sum of two triangles.
         node[v2].vol = node[v2].vol + \
                        tri_area(xj,xm1,xc,yj,ym1,yc) + \
                        tri_area(xj,xc,xm2,yj,yc,ym2);

// - Vertex 3
         xj = node[v3].x;
         yj = node[v3].y;
         xm1 = half*(xj+x4);
         ym1 = half*(yj+y4);
         xm2 = half*(xj+x2);
         ym2 = half*(yj+y2);

//   Median volume is computed as a sum of two triangles.
         node[v3].vol = node[v3].vol + \
                        tri_area(xj,xm1,xc,yj,ym1,yc) + \
                        tri_area(xj,xc,xm2,yj,yc,ym2);

   // - Vertex 4
         xj = node[v4].x;
         yj = node[v4].y;
         xm1 = half*(xj+x1);
         ym1 = half*(yj+y1);
         xm2 = half*(xj+x3);
         ym2 = half*(yj+y3);

//   Median volume is computed as a sum of two triangles.
         node[v4].vol = node[v4].vol + \
                        tri_area(xj,xm1,xc,yj,ym1,yc) + \
                        tri_area(xj,xc,xm2,yj,yc,ym2);
 
   }//    endif tri_or_quadv

}//   end do elementsv

// //--------------------------------------------------------------------------------
// // Loop over elements 2
// //
// //  Allocate elm[:].nghbr[:] : elm[:].nnghrs, elm[:].nghr[:]
// //  Construct element nghbr data: elm(:).nghbr(:)
// //  Order of neighbor elements [e1,e2,e3,..] are closely related to
// //  the order of vertices [v1,v2,v3,..] (see below).
// //
// //          o------o
// //          |      |                
// //        v4|  e1  |v3                     v3
// //    o-----o------o------o      o---------o------------o
// //    |     |      |      |       .      .   .        .
// //    | e2  |      |  e4  |        . e2 .     . e1  .
// //    o-----o------o------o         .  .       .  .
// //       v1 |     .v2              v1 o---------o v2   
// //          | e3 .                     .   e3  .
// //          |   .                        .    .
// //          |  .                           . .
// //          | .                             o
// //          o
// //

   // Allocate the neighbor array

   for (int i = 0; i < nelms; i++) {

   //  3 neighbors for triangle
      if (elm[i].nvtx==3) {

         elm[i].nnghbrs = 3;
         elm[i].nghbr = new Array2D<int>(3,1);
      }
   //  4 neighbors for quadrilateral
      else if (elm[i].nvtx==4) {

         elm[i].nnghbrs = 4;
         elm[i].nghbr = new Array2D<int>(4,1);

      }

   }

// Begin constructing the element-neighbor data

   //elements2 : do i = 1, nelms
   for (int i = 0; i < nelms; i++) {

      //elm_vertex : do k = 1, elm(i).nvtx
      for (int k = 0; k < elm[i].nvtx; k++) {
      //   Get the face of the element i:
      //
      //             vL      vR
      //              o------o
      //             /       |
      //            /        |
      //           o---------o
      //
      //if (k  < elm[i].nvtx)  vL = (*elm[i].vtx)(k+1,0); //TLM warning:  apparent step out of bounds 
         if (k  < elm[i].nvtx)  vL = (*elm[i].vtx)(k); //TLM dumb fix
         if (k == elm[i].nvtx) vL = (*elm[i].vtx)(0); 
         vR = (*elm[i].vtx)(k);
      //   Loop over the surrounding elements of the node vR,
      //   and find the element neighbor from them.
         found = false;
      //elms_around_vR : do j = 1, node(vR).nelms;
         for (int j = 0; j < node[vR].nelms; j++) {
            jelm = (*node[vR].elm)(j); 
         }}}
//          //edge_matching : do ii = 1, elm(jelm).nvtx;
//          // I just remembered that 
//          //  the "2D" array works just as well as 1D
//          for (int ii = 0; ii < elm[jelm].nvtx; ii++) {
//                          v1 = (*elm[jelm].vtx)(ii);
//             if (ii  > 1) { v2 = (*elm[jelm].vtx)(ii-1); }
//             if (ii == 1) { v2 = (*elm[jelm].vtx)(elm[jelm].nvtx); }

//             if (v1==vR and v2==vL) {
//                found = true;
//                im = ii+1;
//                if (im > elm[jelm].nvtx) { im = im - elm[jelm].nvtx; }
//                break; //exit edge_matching
//             } //endif
//          } //end do edge_matching

//      //if (found) exit elms_around_vR
//      if (found) break;

//       } //end do elms_around_vR

//       in = k + 2;
//       if (in > elm[i].nvtx) { in = in - elm[i].nvtx; }

//       if (found) {
//          (*elm[   i].nghbr)(in) = jelm;
//          (*elm[jelm].nghbr)(im) = i;
//       }
//       else {
//          (*elm[   i].nghbr)(in) = 0;
//       }

//       }//    end do elm_vertex

//    }//   end do elements2

// //--------------------------------------------------------------------------------
// // Edge-data for node-centered (edge-based) scheme.
// //
// // Loop over elements 3
// // Construct edge data: edge(:).n1, n2, e1, e2.
// // Edge points from node n1 to node n2.
// //
// //      n2
// //       o------------o
// //     .  \         .
// //    .    \   e2  .
// //   .  e1  \    .
// //  .        \ .         Directed area is positive: n1 -> n2
// // o----------o         e1: left element
// //             n1       e2: right element (e2 > e1 or e2 = 0)

// // First count the number of edges.
// //
// // NOTE: Count edges only if the neighbor element number is
// //       greater than the current element (i) to avoid double
// //       count. Zero element number indicates that it is outside
// //       the domain (boundary face).

//   elements0 : do i = 1, nelms

//    v1 = elm[i].vtx(1)
//    v2 = elm[i].vtx(2)
//    v3 = elm[i].vtx(3)

//    tri_quad0 : if (elm[i].nvtx==3) then

//     if ( elm[i].nghbr(3) > i  .or. elm[i].nghbr(3)==0 ) then
//      nedges = nedges + 1
//     endif

//     if ( elm[i].nghbr(1) > i .or. elm[i].nghbr(1)==0 ) then
//      nedges = nedges + 1
//     endif

//     if ( elm[i].nghbr(2) > i .or. elm[i].nghbr(2)==0 ) then
//      nedges = nedges + 1
//     endif

//    elseif (elm[i].nvtx==4) then

//     v4 = elm[i].vtx(4)

//     if ( elm[i].nghbr(3) > i .or. elm[i].nghbr(3) ==0 ) then
//      nedges = nedges + 1
//     endif

//     if ( elm[i].nghbr(4) > i .or. elm[i].nghbr(4) ==0 ) then
//      nedges = nedges + 1
//     endif

//     if ( elm[i].nghbr(1) > i .or. elm[i].nghbr(1) ==0 ) then
//      nedges = nedges + 1
//     endif

//     if ( elm[i].nghbr(2) > i .or. elm[i].nghbr(2) ==0 ) then
//      nedges = nedges + 1
//     endif

//    endif tri_quad0

//   end do elements0

// // Allocate the edge array.
//   allocate(edge(nedges))
//   nedges = 0
//   edge(:).e1 = 0
//   edge(:).e2 = 0

// // Construct the edge data:
// //  two end nodes (n1, n2), and left and right elements (e1, e2)

//   elements3 : do i = 1, nelms

//    v1 = elm[i].vtx(1)
//    v2 = elm[i].vtx(2)
//    v3 = elm[i].vtx(3)

// // Triangular element
//    tri_quad2 : if (elm[i].nvtx==3) then

//     if ( elm[i].nghbr(3) > i  .or. elm[i].nghbr(3)==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v1
//      edge(nedges).n2 = v2
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(3)
//     endif

//     if ( elm[i].nghbr(1) > i .or. elm[i].nghbr(1)==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v2
//      edge(nedges).n2 = v3
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(1)
//     endif

//     if ( elm[i].nghbr(2) > i .or. elm[i].nghbr(2)==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v3
//      edge(nedges).n2 = v1
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(2)
//     endif

// //  Quadrilateral element
//    elseif (elm[i].nvtx==4) then

//     v4 = elm[i].vtx(4)

//     if ( elm[i].nghbr(3) > i .or. elm[i].nghbr(3) ==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v1
//      edge(nedges).n2 = v2
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(3)
//     endif

//     if ( elm[i].nghbr(4) > i .or. elm[i].nghbr(4) ==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v2
//      edge(nedges).n2 = v3
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(4)
//     endif

//     if ( elm[i].nghbr(1) > i .or. elm[i].nghbr(1) ==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v3
//      edge(nedges).n2 = v4
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(1)
//     endif

//     if ( elm[i].nghbr(2) > i .or. elm[i].nghbr(2) ==0 ) then
//      nedges = nedges + 1
//      edge(nedges).n1 = v4
//      edge(nedges).n2 = v1
//      edge(nedges).e1 = i
//      edge(nedges).e2 = elm[i].nghbr(2)
//     endif

//    endif tri_quad2

//   end do elements3

// // Loop over edges
// // Construct edge vector and directed area vector.
// //
// // Edge vector is a simple vector pointing froom n1 to n2.
// // For each edge, add the directed area vector (dav) from
// // the left and right elements.
// //
// //              n2
// //   o-----------o-----------o
// //   |     dav   |  dav      |
// //   |       ^   |   ^       |
// //   |       |   |   |       |
// //   |   c - - - m - - -c    |
// //   |           |           |
// //   |           |           |    m: edge midpoint
// //   |           |           |    c: element centroid
// //   o-----------o-----------o
// //                n1
// //
//   edges : do i = 1, nedges

//    n1 = edge[i].n1
//    n2 = edge[i].n2
//    e1 = edge[i].e1
//    e2 = edge[i].e2
//    xm = half*( node(n1).x + node(n2).x )
//    ym = half*( node(n1).y + node(n2).y )

//    edge[i].dav = zero

// // Contribution from the left element
//   if (e1 > 0) then
//    xc = elm(e1).x
//    yc = elm(e1).y
//    edge[i].dav(1) = -(ym-yc)
//    edge[i].dav(2) =   xm-xc
//   endif

// // Contribution from the right element
//   if (e2 > 0) then
//    xc = elm(e2).x
//    yc = elm(e2).y
//    edge[i].dav(1) = edge[i].dav(1) -(yc-ym)
//    edge[i].dav(2) = edge[i].dav(2) + xc-xm
//   endif

//   if (e1 < 0 .and. e2 < 0) then
//    write(*,*) "////////// e1 and e2 are both negative... No way..."
//   endif

// // Magnitude and unit vector
//    edge[i].da  = sqrt( edge[i].dav(1)**2 + edge[i].dav(2)**2 )
//    edge[i].dav = edge[i].dav / edge[i].da

// // Edge vector

//   edge[i].ev(1) = node(n2).x - node(n1).x
//   edge[i].ev(2) = node(n2).y - node(n1).y
//   edge[i].e     = sqrt( edge[i].ev(1)**2 + edge[i].ev(2)**2 )
//   edge[i].ev    = edge[i].ev / edge[i].e

//   end do edges

// //--------------------------------------------------------------------------------
// // Construct node neighbor data:
// //  pointers to the neighbor nodes(o)
// //
// //        o     o
// //         \   / 
// //          \ /
// //     o-----*-----o
// //          /|
// //         / |
// //        /  o        *: node in interest
// //       o            o: neighbors (edge-connected nghbrs)
// //

//   do i = 1, nnodes
//    node[i].nnghbrs = 0
//   end do

// // Loop over edges and distribute the node numbers:

//   edges4 : do i = 1, nedges

//    n1 = edge[i].n1
//    n2 = edge[i].n2

// // (1) Add node1 to the neighbor list of n2
//    node(n1).nnghbrs = node(n1).nnghbrs + 1
//    call my_alloc_int_ptr(node(n1).nghbr, node(n1).nnghbrs)
//    node(n1).nghbr(node(n1).nnghbrs) = n2

// // (2) Add node2 to the neighbor list of n1
//    node(n2).nnghbrs = node(n2).nnghbrs + 1
//    call my_alloc_int_ptr(node(n2).nghbr, node(n2).nnghbrs)
//    node(n2).nghbr(node(n2).nnghbrs) = n1

//   end do edges4

// //--------------------------------------------------------------------------------
// // Boundary normal at nodes constructed by accumulating the contribution
// // from each boundary face normal. This vector will be used to enforce
// // the tangency condition, for example.
// //
// //
// //        Interior domain      /
// //                            o
// //                  .        /
// //                  .       /
// // --o-------o-------------o
// //           j   |  .  |   j+1
// //               v  .  v
// //
// //        Left half added to the node j, and
// //       right half added to the node j+1.
// //

// // Allocate and initialize the normal vector arrays
//   do i = 1, nbound

//    allocate(bound[i].bnx(bound[i].nbnodes))
//    allocate(bound[i].bny(bound[i].nbnodes))
//    allocate(bound[i].bn( bound[i].nbnodes))

//    do j = 1, bound[i].nbnodes
//     bound[i].bnx(j) = zero
//     bound[i].bny(j) = zero
//     bound[i].bn( j) = zero
//    end do

//   end do

// // Normal vector at boundary nodes
// // Note: Below it describes normals of linear approximation.
// //       We will overwrite it by a quadratic approximation.
// //
// // Linear approximation:
// //
// // Step 1. Compute the outward normals
//   do i = 1, nbound
//    do j = 1, bound[i].nbnodes-1

//     x1 = node(bound[i].bnode(j  )).x
//     y1 = node(bound[i].bnode(j  )).y

//     x2 = node(bound[i].bnode(j+1)).x
//     y2 = node(bound[i].bnode(j+1)).y

// //   Normal vector pointing into the domain at this point.
//     bound[i].bnx(j) = bound[i].bnx(j) + half*( -(y2-y1) )
//     bound[i].bny(j) = bound[i].bny(j) + half*(   x2-x1  )

//     bound[i].bnx(j+1) = bound[i].bnx(j+1) + half*( -(y2-y1) )
//     bound[i].bny(j+1) = bound[i].bny(j+1) + half*(   x2-x1  )

//    end do
//   end do

// // Step 2. Compute the magnitude and turn (bnx,bny) into a unit vector
//   do i = 1, nbound
//    do j = 1, bound[i].nbnodes

//     bound[i].bn(j)  = sqrt( bound[i].bnx(j)**2 + bound[i].bny(j)**2 )
// //   Minus sign to make it pont out towards the outside of the domain.
//     bound[i].bnx(j) =  - bound[i].bnx(j) / bound[i].bn(j)
//     bound[i].bny(j) =  - bound[i].bny(j) / bound[i].bn(j)

//    end do
//   end do

// // Now, ignore the linear approximation, and let us construct
// // more accurate surfae normal vectors and replace the linear ones.
// // So, we will overwrite the unit normal vectors: bnx, bny.
// // Note: We keep the magnitude of the normal vector.
// //
// // Quadratic approximation:
// // See http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2015v281pp518-555_preprint.pdf
// // for details on the quadratic approximation for computing more accurate normals.
// // 
//   boundary_type0 : do i = 1, nbound
//    boundary_nodes0 : do j = 1, bound[i].nbnodes

//      if (j==1) then
//       v1 = bound[i].bnode(j  )
//       v2 = bound[i].bnode(j+1)
//       v3 = bound[i].bnode(j+2)
//      elseif (j==bound[i].nbnodes) then
//       v1 = bound[i].bnode(j-2)
//       v2 = bound[i].bnode(j-1)
//       v3 = bound[i].bnode(j  )
//      else
//       v1 = bound[i].bnode(j-1)
//       v2 = bound[i].bnode(j)
//       v3 = bound[i].bnode(j+1)
//      endif

//      x1 = node(v1).x
//      x2 = node(v2).x
//      x3 = node(v3).x

//      y1 = node(v1).y
//      y2 = node(v2).y
//      y3 = node(v3).y

// //----------------------------------------------------------------------
// //   Fit a quadratic over 3 nodes

// //    Skip the last one if the boundary segment is a closed boundary 
// //    in which case the last node is the same as the first one.
//      if (j==bound[i].nbnodes .and. bound[i].bnode(j)==bound[i].bnode(1) ) then
//       bound[i].bn(j)  = bound[i].bn(1)
//       bound[i].bnx(j) = bound[i].bnx(1)
//       bound[i].bny(j) = bound[i].bny(1)
//       cycle
//      endif

//      dsL = sqrt( (x2-x1)**2 + (y2-y1)**2 )
//      dsR = sqrt( (x3-x2)**2 + (y3-y2)**2 )
//       dx = dsR*x1/(dsL*(-dsL-dsR))-x2/dsR+x2/dsL+dsL*x3/((dsR+dsL)*dsR)
//       dy = dsR*y1/(dsL*(-dsL-dsR))-y2/dsR+y2/dsL+dsL*y3/((dsR+dsL)*dsR)

//      ds  = sqrt( dx**2 + dy**2 )
//      bound[i].bnx(j) = -( -dy / ds )
//      bound[i].bny(j) = -(  dx / ds )

//    end do boundary_nodes0
//   end do boundary_type0

// //--------------------------------------------------------------------------------
// // Construct neighbor index over edges
// //
// //  Example:
// //
// //        o     o
// //         \   / 
// //          \j/       k-th neighbor
// //     o-----*----------o
// //          /|  edge i
// //         / |
// //        /  o        Note: k-th neighbor is given by "node(j).nghbr(k)"
// //       o
// //
// //  Consider the edge i
// //
// //   node j        k-th neighbor
// //       *----------o
// //      n1  edge i  n2
// //
// //   We store "k" in the edge data structure as
// //
// //    edge[i].kth_nghbr_of_1: n2 is the "edge[i].kth_nghbr_of_1"-th neighbor of n1
// //    edge[i].kth_nghbr_of_2: n1 is the "edge[i].kth_nghbr_of_3"-th neighbor of n2
// //
// //   That is,  we have
// //
// //    n2 = node(n1).nghbr(edge[i].kth_nghbr_of_1)
// //    n1 = node(n2).nghbr(edge[i].kth_nghbr_of_2)
// //
// //   We make use of this data structure to access off-diagonal entries in Jacobian matrix.
// //

// // Loop over edges

//   edges5 : do i = 1, nedges

//    n1 = edge[i].n1
//    n2 = edge[i].n2

//    do k = 1, node(n2).nnghbrs

//     if ( n1 == node(n2).nghbr(k) ) then
//      edge[i].kth_nghbr_of_2 = k
//     endif

//    end do

//    do k = 1, node(n1).nnghbrs

//     if ( n2 == node(n1).nghbr(k) ) then
//      edge[i].kth_nghbr_of_1 = k
//     endif

//    end do

//   end do edges5

// // Boundary mark: It should be an array actually because some nodes are associated with
// //                more than one boundaries.
//   do i = 1, nnodes
//    node[i].bmark   = 0
//    node[i].nbmarks = 0
//   end do

//   do i = 1, nbound
//    do j = 1, bound[i].nbnodes
//     node( bound[i].bnode(j) ).bmark   = i
//     node( bound[i].bnode(j) ).nbmarks = node( bound[i].bnode(j) ).nbmarks + 1
//    end do
//   end do

// //--------------------------------------------------------------------------------
// // Boundary face data
// //
// //      |     Domain      |
// //      |                 |
// //      o--o--o--o--o--o--o  <- Boundary segment
// //   j= 1  2  3  4  5  6  7
// //
// //   In the above case, nbnodes = 7, nbfaces = 6
// //

//   do i = 1, nbound
//    bound[i].nbfaces = bound[i].nbnodes-1
//    allocate(bound[i].bfnx(    bound[i].nbfaces   ))
//    allocate(bound[i].bfny(    bound[i].nbfaces   ))
//    allocate(bound[i].bfn(     bound[i].nbfaces   ))
//    allocate(bound[i].belm(    bound[i].nbfaces   ))
//    allocate(bound[i].kth_nghbr_of_1(    bound[i].nbfaces   ))
//    allocate(bound[i].kth_nghbr_of_2(    bound[i].nbfaces   ))
//   end do

// // Boundary face vector: outward normal
//   do i = 1, nbound
//    do j = 1, bound[i].nbfaces

//     x1 = node(bound[i].bnode(j  )).x
//     y1 = node(bound[i].bnode(j  )).y
//     x2 = node(bound[i].bnode(j+1)).x
//     y2 = node(bound[i].bnode(j+1)).y

//     bound[i].bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
//     bound[i].bfnx(j) = -(y1-y2) / bound[i].bfn(j)
//     bound[i].bfny(j) =  (x1-x2) / bound[i].bfn(j)

//    end do
//   end do

// // Boundary normal vector at nodes: outward normal
//   do i = 1, nbound
//    do j = 1, bound[i].nbfaces

//     x1 = node(bound[i].bnode(j  )).x
//     y1 = node(bound[i].bnode(j  )).y
//     x2 = node(bound[i].bnode(j+1)).x
//     y2 = node(bound[i].bnode(j+1)).y

//     bound[i].bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
//     bound[i].bfnx(j) = -(y1-y2) / bound[i].bfn(j)
//     bound[i].bfny(j) =  (x1-x2) / bound[i].bfn(j)

//    end do
//   end do

// // Neighbor index over boundary edges (faces)

//   do i = 1, nbound
//    do j = 1, bound[i].nbfaces

//     n1 = bound[i].bnode(j  )  //Left node
//     n2 = bound[i].bnode(j+1)  //Right node

//     do k = 1, node(n2).nnghbrs
//      if ( n1 == node(n2).nghbr(k) ) then
//       bound[i].kth_nghbr_of_2(j) = k
//      endif
//     end do

//     do k = 1, node(n1).nnghbrs
//      if ( n2 == node(n1).nghbr(k) ) then
//       bound[i].kth_nghbr_of_1(j) = k
//      endif
//     end do

//    end do
//   end do

// // Find element adjacent to the face: belm
// //
// //  NOTE: This is useful to figure out what element
// //        each boundary face belongs to. Boundary flux needs
// //        special weighting depending on the element.
// //
// //      |_________|_________|________|
// //      |         |         |        | 
// //      |         |         |        | 
// //      |_________|_________|________|
// //      |         |         |        |     <- Grid (e.g., quads)
// //      |         | elmb(j) |        |
// //   ---o---------o---------o--------o---  <- Boundary segment
// //                 j-th face
// //
// // elmb(j) is the element number of the element having the j-th boundary face.
// //

//   do i = 1, nbound
//    do j = 1, bound[i].nbfaces

// //   bface is defined by the nodes v1 and v2.
//     v1 = bound[i].bnode(j  )
//     v2 = bound[i].bnode(j+1)

//     found = .false.
// //   Find the element having the bface from the elements
// //   around the node v1.
//     do k = 1, node(v1).nelms
//      ielm = node(v1).elm(k)
//      do ii = 1, elm(ielm).nvtx
//       in = ii
//       im = ii+1
//       if (im > elm(ielm).nvtx) im = im - elm(ielm).nvtx //return to 1
//       vt1 = elm(ielm).vtx(in)
//       vt2 = elm(ielm).vtx(im)
//        if (vt1 == v1 .and. vt2 == v2) then
//         found = .true.
//         exit
//        endif
//      end do
//      if (found) exit
//     end do

//     if (found) then
//      bound[i].belm(j) = ielm
//     else
//      write(*,*) " Boundary-adjacent element not found. Error..."
//      stop
//     endif

//    end do
//   end do

// //--------------------------------------------------------------------------------
// // Construct least-squares matrix for node-centered schemes.
// //
// //        o     o
// //         \   / 
// //          \ /
// //     o-----*-----o
// //          /|
// //         / |
// //        /  o        *: node in interest
// //       o            o: neighbors (edge-connected nghbrs)
// //

// // Check the number of neighbor nodes (must have at least 2 neighbors)
//   write(*,*) " --- Node neighbor data:"

//   ave_nghbr = node(1).nnghbrs
//   min_nghbr = node(1).nnghbrs
//   max_nghbr = node(1).nnghbrs
//        imin = 1
//        imax = 1
//    if (node(1).nnghbrs==2) then
//     write(*,*) "--- 2 neighbors for the node = ", 1
//    endif

//   do i = 2, nnodes
//    ave_nghbr = ave_nghbr + node[i].nnghbrs
//    if (node[i].nnghbrs < min_nghbr) imin = i
//    if (node[i].nnghbrs > max_nghbr) imax = i
//    min_nghbr = min(min_nghbr, node[i].nnghbrs)
//    max_nghbr = max(max_nghbr, node[i].nnghbrs)
//    if (node[i].nnghbrs==2) then
//     write(*,*) "--- 2 neighbors for the node = ", i
//    endif
//   end do

//   write(*,*) "      ave_nghbr = ", ave_nghbr/nnodes
//   write(*,*) "      min_nghbr = ", min_nghbr, " at node ", imin
//   write(*,*) "      max_nghbr = ", max_nghbr, " at node ", imax
//   write(*,*)

// //--------------------------------------------------------------------------------
// // Cell centered scheme data
// //--------------------------------------------------------------------------------

//   write(*,*) "Generating CC scheme data......"

//   do i = 1, nelms   
//    allocate(elm(i).edge( elm(i).nnghbrs ) )
//   end do

//   edges3 : do i = 1, nedges

//    e1 = edge[i].e1
//    e2 = edge[i].e2

// // Left element
//   if (e1 > 0) then
//    do k = 1, elm(e1).nnghbrs
//     if (elm(e1).nghbr(k)==e2) elm(e1).edge(k) = i
//    end do
//   endif

// // Right element
//   if (e2 > 0) then
//    do k = 1, elm(e2).nnghbrs
//     if (elm(e2).nghbr(k)==e1) elm(e2).edge(k) = i
//    end do
//   endif

//   end do edges3

// // Face-data for cell-centered (edge-based) scheme.
// //
// // Loop over elements 4
// // Construct face data:
// // face is an edge across elements pointing
// // element e1 to element e2 (e2 > e1):
// //
// //       e2
// //        \    
// //         \ face: e1 -> e2 
// //          \
// //  n1 o--------------o n2 <-- face
// //            \
// //             \          n1, n2: end nodes of the face
// //              \         e1: element 1
// //              e1        e2: element 2  (e2 > e1)
// //
// // Note: Face data is dual to the edge data.
// //       It can be trivially constructed from the edge data, but
// //       here the face data is constructed by using the element
// //       neighbor data just for an educational purpose.

//   nfaces = 0
//   elements4 : do i = 1, nelms
//    do k = 1, elm(i).nnghbrs
//    jelm = elm(i).nghbr(k)
//     if (jelm > i) then
//      nfaces = nfaces + 1
//     endif
//    end do
//   end do elements4

//   allocate(face(nfaces))

//   nfaces = 0

//   elements5 : do i = 1, nelms
//    do k = 1, elm(i).nnghbrs
//    jelm = elm(i).nghbr(k)

//     if (jelm > i) then

//      nfaces = nfaces + 1

//      face(nfaces).e1 = i
//      face(nfaces).e2 = jelm

//      iedge = elm(i).edge(k)
//      v1 = edge(iedge).n1
//      v2 = edge(iedge).n2

//      if (edge(iedge).e1 == jelm) then
//       face(nfaces).n1 = v1
//       face(nfaces).n2 = v2
//      else
//       face(nfaces).n1 = v2
//       face(nfaces).n2 = v1
//      endif
   
//     elseif (jelm == 0) then
// //    Skip boundary faces.
//     endif

//    end do
//   end do elements5

// // Loop over faces
// // Construct directed area vector.

//   faces : do i = 1, nfaces

//    n1 = face(i).n1
//    n2 = face(i).n2
//    e1 = face(i).e1
//    e2 = face(i).e2

// // Face vector
//   face(i).dav(1) = -( node(n2).y - node(n1).y )
//   face(i).dav(2) =    node(n2).x - node(n1).x
//   face(i).da     = sqrt( face(i).dav(1)**2 + face(i).dav(2)**2 )
//   face(i).dav    = face(i).dav / face(i).da

//   end do faces

// // Construct vertex-neighbor data for cell-centered scheme.
// //
// // For each element, i, collect all elements sharing the nodes
// // of the element, i, including face-neighors.
// //
// //      ___________
// //     |     |     |
// //     |  o  |  o  |
// //     |_____|_____|
// //    /\    / \    \
// //   / o\ o/ i \  o \
// //  /____\/_____\____\
// //  \    /      /\    \
// //   \o /  o   / o\ o  \
// //    \/______/____\____\
// //
// //          i: Element of interest
// //          o: Vertex neighbors (k = 1,2,...,9)

//   write(*,*) " --- Vertex-neighbor data:"

//   do i = 1, nelms
//    elm(i).nvnghbrs = 1
//    call my_alloc_int_ptr(elm(i).vnghbr, 1)
//   end do

//   ave_nghbr = 0
//   min_nghbr = 10000
//   max_nghbr =-10000
//        imin = 1
//        imax = 1

// // Initialization
//   elements6 : do i = 1, nelms
//    elm(i).nvnghbrs = 0
//   end do elements6

// // Collect vertex-neighbors
//   elements7 : do i = 1, nelms

// // (1)Add face-neighbors
//    do k = 1, elm(i).nnghbrs
//     if ( elm(i).nghbr(k) > 0 ) then
//      elm(i).nvnghbrs = elm(i).nvnghbrs + 1
//      call my_alloc_int_ptr(elm(i).vnghbr, elm(i).nvnghbrs)
//      elm(i).vnghbr(elm(i).nvnghbrs) = elm(i).nghbr(k)
//     endif
//    end do

// // (2)Add vertex-neighbors
//    do k = 1, elm(i).nvtx
//     v1 = elm(i).vtx(k)

//     velms : do j = 1, node(v1).nelms
//      e1 = node(v1).elm(j)
//      if (e1 == i) cycle velms

// //    Check if the element is already added.
//        found = .false.
//      do ii = 1, elm(i).nvnghbrs
//       if ( e1 == elm(i).vnghbr(ii) ) then
//        found = .true.
//        exit
//       endif
//      end do

// //    Add the element, e1, if not added yet.
//      if (.not.found) then
//       elm(i).nvnghbrs = elm(i).nvnghbrs + 1
//       call my_alloc_int_ptr(elm(i).vnghbr, elm(i).nvnghbrs)
//       elm(i).vnghbr(elm(i).nvnghbrs) = e1
//      endif
//     end do velms

//    end do

//    ave_nghbr = ave_nghbr + elm(i).nvnghbrs
//    if (elm(i).nvnghbrs < min_nghbr) imin = i
//    if (elm(i).nvnghbrs > max_nghbr) imax = i
//    min_nghbr = min(min_nghbr, elm(i).nvnghbrs)
//    max_nghbr = max(max_nghbr, elm(i).nvnghbrs)
//    if (elm(i).nvnghbrs < 3) then
//     write(*,*) "--- Not enough neighbors: elm = ", i, &
//                "elm(i).nvnghbrs=",elm(i).nvnghbrs
//    endif

//   end do elements7

//   write(*,*) "      ave_nghbr = ", ave_nghbr/nelms
//   write(*,*) "      min_nghbr = ", min_nghbr, " elm = ", imin
//   write(*,*) "      max_nghbr = ", max_nghbr, " elm = ", imax
//   write(*,*)


//   do i = 1, nelms
//    elm(i).bmark = 0
//   end do

//   bc_loop : do i = 1, nbound
//    if (trim(bound[i].bc_type) == "dirichlet") then
//     do j = 1, bound[i].nbfaces
//      elm( bound[i].belm(j) ).bmark = 1
//     end do
//    endif
//   end do bc_loop

// //--------------------------------------------------------------------------------

   return;
};
//  end subroutine construct_grid_data

//********************************************************************************

