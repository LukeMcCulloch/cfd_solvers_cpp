/*
*
* THE MESHING MODULE
*
*/
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

//======================================
// string trimfunctions
#include "StringOps.h" 

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
//*    cout <<  nnodes, ntria, nquad //Numbers of nodes, triangles and quads
//*
//*   do i = 1, nnodes
//*    cout <<  x(i), y(i) //(x,y) coordinates of each node
//*   end do
//*
//*   do i = 1, ntria        //Nodes of triangles ordered counterclockwise
//*    cout <<  node_1(i), node_2(i), node_3(i)
//*   end do
//*
//*   do i = 1, nquad        //Nodes of quadrilaterals ordered counterclockwise
//*    cout <<  node_1(i), node_2(i), node_3(i), node_4(i)
//*   end do
//* 
//*    cout <<  nbound     //Number of boundary segments
//*
//*   do i = 1, nbound
//*    cout <<  nbnodes(i) //Number of nodes on each segment
//*   end do
//*
//*   do i = 1, nbound
//*    do j = 1, nbnodes(i)
//*     cout <<  bnode(j)  //Node number of each node j in segment i
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
//*    cout <<  "Boundary Segment              Boundary Condition"
//*   do i = 1, nbound
//*    cout <<  i, bc_name
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
   //std::cout << "Allocating node_type" << std::endl;
   //std::cout << "    for " << nnodes << " nodes " << std::endl;
   node = new node_type[nnodes];


   //std::cout << " Allocating elm_type" << std::endl;
   //std::cout << "    for " << nelms << " elements " << std::endl;
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
      //if (i< 200) std::cout  << i << "  " << node[i].x << " " << node[i].y << std::endl;

   }
   //cout << "done reading nodal coords" << endl;

      
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
         elm[i].vtx = new Array2D<int>(3,1);

         //std::string type;
         //in >> type;                  //and read the first whitespace-separated token


         int x, y, z;
         in >> x >> y >> z;       //now read the whitespace-separated ints
         // Fix indices for 0 indexed code//
         (*elm[i].vtx)(0) = x-1;
         (*elm[i].vtx)(1) = y-1;
         (*elm[i].vtx)(2) = z-1;
         //cout << (*elm[i].vtx).ncols << endl;
         
         // if (i<20) cout << "\nInput = " << (*elm[i].vtx)(0) <<
         //                                 "  " << (*elm[i].vtx)(1) << 
         //                                 "  " << (*elm[i].vtx)(2) ;
         // if (i>ntria-20) cout << "\nInput = " << (*elm[i].vtx)(0) <<
         //                                 "  " << (*elm[i].vtx)(1) << 
         //                                 "  " << (*elm[i].vtx)(2) ;
               
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
         // Fix indices for 0 indexed code//
         (*elm[ntria+i-1].vtx)(0) = x1-1;
         (*elm[ntria+i-1].vtx)(1) = x2-1;
         (*elm[ntria+i-1].vtx)(2) = x3-1;
         (*elm[ntria+i-1].vtx)(3) = x4-1;
         // if (i<20) cout << "\nx, y, z = " << x1 <<"  " << x2 << "  " << x3 << x4;
         // if (i<20) cout << "\nelm x, elm y, elm z = " << (*elm[i].vtx)(0,0) <<"  " << (*elm[i].vtx)(1,0) << "  " << (*elm[i].vtx)(2,0)<< "  " << (*elm[i].vtx)(3,0);
         // if (i<20) cout << "\nelm x, elm y, elm z = " << (*elm[i].vtx)(0) <<"  " << (*elm[i].vtx)(1) << "  " << (*elm[i].vtx)(2) << "  " << (*elm[i].vtx)(3);
      }
   }
   else{
      cout << "\nNo Quads in this Mesh\n" << endl;
   }

   //  Write out the grid data.
   cout << " " << endl;
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
   cout << "\n nbound = " << nbound << "\n" <<endl;


   // // READ: Number of Boundary nodes (including the starting one at the end if
   // // it is closed such as an airfoil.)
   for (size_t i = 0; i < nbound; i++) {
      std::getline(infile, line);
      std::istringstream in(line);
      in >> bound[i].nbnodes;
      //cout << "Got in bnodes = " << bound[i].nbnodes << endl;
      bound[i].bnode = new Array2D<int>(bound[i].nbnodes , 1);
   }

   // // READ: Read boundary nodes
   cout << "Reading boundary nodes" << endl;
   std::getline(infile, line); //TLM: need to skip line here
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {
         std::getline(infile, line);
         std::istringstream in(line);
         int init;
         in >> init;
         //(*bound[i].bnode)[j][0] = init;
         //TLM bnode issue : we are reading in indices that 
         //  were starting from 1, but our code starts form 0//
         (*bound[i].bnode)(j,0) = init-1;
         
         
         //cout << i << " " << j << " " << (*bound[i].bnode)(j,0)  << endl;

         // if ( j>bound[i].nbnodes-21 ) {
         //    cout << (*bound[i].bnode)(j,0) << endl;
         // }

         //cout << (*bound[i].bnode)(j,0) << "   " << init << endl;
         // testing array access:
         //int some = (*bound[i].bnode)[j][0];
         //cout << "get some " << some << endl;
      }
   }
   cout << "Done Reading boundary nodes" << endl;

   //  Print the boundary grid data.
   std::cout << " Boundary nodes:" << std::endl;
   std::cout << "    segments = " << nbound << std::endl;
      for (size_t i = 0; i < nbound; i++) {
         std::cout <<  " boundary = " << i << 
                     "   bnodes = " <<  bound[i].nbnodes <<  
                     "   bfaces = " <<  bound[i].nbnodes-1 << std::endl;
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
   std::ifstream outfile;
   outfile.open (datafile_bcmap_in);

   std::getline(outfile, line);

   // READ: Read the boundary condition type
   for (size_t i = 0; i < nbound; i++) {
      std::getline(outfile, line);
      std::istringstream in(line);
      in >> dummy_int >> bound[i].bc_type;
   }

   //  Print the data
   std::cout << " Boundary conditions:" << std::endl;
   for (size_t i = 0; i < nbound; i++) {
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
   vL = -1;
   im = 0;
   jelm = 0;

   cout << " " << endl;
   cout << "construct grid data...." << endl;

   // // Initializations
   for (size_t i = 0; i < nnodes; i++) {
      node[i].nelms = 0;
   } 
   nedges = 0;
   cout << "\nnnodes = " << nnodes << endl;
   cout << "nelms = " << nelms << "\n" << endl;

//--------------------------------------------------------------------------------
// Loop over elements and construct the fololowing data.
//
// 1. Surrounding elements: node(:).nelms, node(:).elm(:)
//
//    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
//             node[i].nelms = 4
//             node[i].elm(1) = 23
//             node[i].elm(2) = 13
//             node[i].elm(3) = 41
//             node[i].elm(4) = 101
//
//        o-------o-------------o
//       /        |   .         |
//      /    23   |      41     |
//     o----------o-------------o
//      \        i \            |
//       \   101    \     13    |
//        \          \          | 
//         o----------o---------o
//
// 2. Element quantities  : elm(:).x,elm(:).y,elm(:).vol
//
//  o-----------o            
//   \          |            o
//    \    (x,y)|           / \
//     \   .    |          /   \
//      \       |         /  .  \    (x,y): centroid coordinates
//       \      |        / (x,y) \     vol: volume of element
//        o-----o       o---------o

//   elements : do i = 1, nelms

   // for ( int i = 0; i < nelms; ++i ) {
   //    v1 = (*elm[i].vtx)(0,0);
   //    v2 = (*elm[i].vtx)(1,0);
   //    v3 = (*elm[i].vtx)(2,0);
   //    node[v1].nelms = node[v1].nelms + 1;
   //    node[v2].nelms = node[v2].nelms + 1;
   //    node[v3].nelms = node[v3].nelms + 1;

   // }
   // for ( int i = 0; i < nelms; ++i ) {
   //    v1 = (*elm[i].vtx)(0,0);
   //    v2 = (*elm[i].vtx)(1,0);
   //    v3 = (*elm[i].vtx)(2,0);
   //    node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
   //    node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
   //    node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
   // }

   //dummy allocations:
   // node[v1].elm = new Array2D<int>(1, 1);
   // node[v2].elm = new Array2D<int>(1, 1);
   // node[v3].elm = new Array2D<int>(1, 1);
   for ( int i = 0; i < nelms; ++i ) {
      v1 = (*elm[i].vtx)(0,0);
      v2 = (*elm[i].vtx)(1,0);
      v3 = (*elm[i].vtx)(2,0);
      // node[v1].elm = new Array2D<int>(1, 1);
      // node[v2].elm = new Array2D<int>(1, 1);
      // node[v3].elm = new Array2D<int>(1, 1);
   }

   for ( int i = 0; i < nelms; ++i ) {

      v1 = (*elm[i].vtx)(0,0);
      v2 = (*elm[i].vtx)(1,0);
      v3 = (*elm[i].vtx)(2,0);

      x1 = node[v1].x;
      x2 = node[v2].x;
      x3 = node[v3].x;

      y1 = node[v1].y;
      y2 = node[v2].y;
      y3 = node[v3].y;

   // Distribute the element index to nodes.

      /*
      * DESIGN CHOICE HERE -- use std<vector> or similar?
      * should use containers
      */



      // save and reallocate:
      node[v1].nelms = node[v1].nelms + 1;
      if (node[v1].nelms == 1) {
         node[v1].elm = new Array2D<int>(1, 1);
      }
      else {
         Array2D<int> save1 =  (*node[v1].elm); // copy constructor
         node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
         for (size_t ra = 0; ra < save1.storage_size; ra++){
            (*node[v1].elm).array[ra] = save1.array[ra];
         }
      }
      (*node[v1].elm)(node[v1].nelms-1, 0) = i;
      node[v1].elmV.append(i);



      // save and reallocate:
      node[v2].nelms = node[v2].nelms + 1;
      if (node[v2].nelms == 1) {
         node[v2].elm = new Array2D<int>(1, 1);
      }
      else {
         Array2D<int> save2 =  (*node[v2].elm); // copy constructor
         node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
         for (size_t ra = 0; ra < save2.storage_size; ra++){
            (*node[v2].elm).array[ra] = save2.array[ra];
         }
      }
      (*node[v2].elm)(node[v2].nelms-1, 0) = i;
      node[v2].elmV.append(i);



      // // save and reallocate:
      node[v3].nelms = node[v3].nelms + 1;
      if (node[v3].nelms == 1) {
         node[v3].elm = new Array2D<int>(1, 1);
      }
      else {
         Array2D<int> save3 =  (*node[v3].elm); // copy constructor
         node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
         for (size_t ra = 0; ra < save3.storage_size; ra++){
            (*node[v3].elm).array[ra] = save3.array[ra];
         }
      }
      (*node[v3].elm)(node[v3].nelms-1, 0) = i;
      node[v3].elmV.append(i);
      // end save and reallocate:


      // Compute the cell center and cell volume.
      //tri_or_quad : if (elm(i).nvtx==3) then
      if (elm[i].nvtx==3) {

         // Triangle centroid and volume
         elm[i].x   = third*(x1+x2+x3);
         elm[i].y   = third*(y1+y2+y3);
         //cout << " tri area -1" << endl;
         elm[i].vol = tri_area(x1,x2,x3,y1,y2,y3);
      }
      else if (elm[i].nvtx==4) {

         cout << "need to fix reallocator for 4th node \n";
         std::exit(0);
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
         //cout << " tri area 0" << endl;
         elm[i].vol = tri_area(x1,x2,x3,y1,y2,y3) + \
                     tri_area(x1,x3,x4,y1,y3,y4);

         xc = elm[i].x;
         yc = elm[i].y;
         //cout << " tri area 1" << endl;
         if (tri_area(x1,x2,xc,y1,y2,yc)<=zero) {
            cout << " Centroid outside the quad element 12c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            std::exit(0);
         }

         //cout << " tri area 2" << endl;
         if (tri_area(x2,x3,xc,y2,y3,yc)<=zero) {
            cout << " Centroid outside the quad element 23c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            std::exit(0);
         }

         //cout << " tri area 3" << endl;
         if (tri_area(x3,x4,xc,y3,y4,yc)<=zero) {
            cout << " Centroid outside the quad element 34c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1 << endl;
            cout << "  (x2,y2)=" << x2 << y2 << endl;
            cout << "  (x3,y3)=" << x3 << y3 << endl;
            cout << "  (x4,y4)=" << x4 << y4 << endl;
            cout << "  (xc,yc)=" << xc << yc << endl;
            std::exit(0);
         }

         //cout << " tri area 4" << endl;
         if (tri_area(x4,x1,xc,y4,y1,yc)<=zero) {
            cout << " Centroid outside the quad element 41c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            std::exit(0);
         }

      //  Distribution of element number to the 4th node of the quadrilateral
         // node[v4].nelms = node[v4].nelms + 1;
         // node[v4].elm = new Array2D<int>(node[v4].nelms, 1);
         // (*node[v4].elm)[ node[v4].nelms-1 ][0] = i;

         // easier way TODO: fixme for quads
         (*node[v4].elm)((*node[v4].elm).tracked_index, 0) = i;
         (*node[v4].elm).tracked_index +=1;

         // node[v4].elm = new Array2D<int>(node[v4].nelms+1, 1);
         // (*node[v4].elm)( node[v4].nelms , 0) = i;
         // node[v4].nelms = node[v4].nelms + 1;

      }//    endif tri_or_quad
      else {
         cout << "ERROR: not a tri or quad" << endl;
         std:exit(0);
      }

   }//   end do elements (i loop)

// Median dual volume

   for (size_t i = 0; i < nnodes; i++) {
      node[i].vol = zero;
   }

//   elementsv : do i = 1, nelms
   for ( int i = 0; i < nelms; ++i ) {
      
//TLM here 2/23/2020 6::11
      v1 = (*elm[i].vtx)(0);
      v2 = (*elm[i].vtx)(1);
      v3 = (*elm[i].vtx)(2);

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

//--------------------------------------------------------------------------------
// Loop over elements 2
//
//  Allocate elm[:].nghbr[:] : elm[:].nnghrs, elm[:].nghr[:]
//  Construct element nghbr data: elm(:).nghbr(:)
//  Order of neighbor elements [e1,e2,e3,..] are closely related to
//  the order of vertices [v1,v2,v3,..] (see below).
//
//          o------o
//          |      |                
//        v4|  e1  |v3                     v3
//    o-----o------o------o      o---------o------------o
//    |     |      |      |       .      .   .        .
//    | e2  |      |  e4  |        . e2 .     . e1  .
//    o-----o------o------o         .  .       .  .
//       v1 |     .v2              v1 o---------o v2   
//          | e3 .                     .   e3  .
//          |   .                        .    .
//          |  .                           . .
//          | .                             o
//          o
//

   // Allocate the neighbor array
   // narrow stencil

   for (size_t i = 0; i < nelms; i++) {

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

      (*elm[   i].nghbr) = -1;

   }
int nbrprint = 2;
// Begin constructing the element-neighbor data
   cout << "Begin constructing the element-neighbor data \n" << endl;
   //elements2 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {

      //elm_vertex : do k = 1, elm(i).nvtx
      for (size_t k = 0; k < elm[i].nvtx; k++) {
         //   Get the face of the element i:
         //
         //             vL      vR
         //              o------o
         //             /       |
         //            /        |
         //           o---------o
         //
         //TLM warning:  fixed step out of bounds
         if (k  < elm[i].nvtx-1) vL = (*elm[i].vtx)(k+1); //k+1 - 1.. nope, K goes from 0
         if (k == elm[i].nvtx-1) vL = (*elm[i].vtx)(0); //1-1
         vR = (*elm[i].vtx)(k);

         //   Loop over the surrounding elements of the node vR,
         //   and find the element neighbor from them.
         found = false;
         //elms_around_vR
         for (size_t j = 0; j < node[vR].nelms; j++) {
            jelm = (*node[vR].elm)(j);
            
            // if (i < nbrprint) cout <<  " i = " << i <<  " j = " << j << endl;
            // if (i < nbrprint) cout << "vR = " << vR << " " << "   jelm = " << (*node[vR].elm)(j) << endl;
            //if (i < nbrprint) cout << "vR , jelm = "<< vR << "   " << jelm << endl;

            //edge_matching
            for (size_t ii = 0; ii < elm[jelm].nvtx; ii++) {
               
               v1 = (*elm[jelm].vtx)(ii);
               //cout << ii << endl;
               if (ii  > 0) { 
                  v2 = (*elm[jelm].vtx)(ii-1); 
               }
               if (ii == 0) { 
                  v2 = (*elm[jelm].vtx)(elm[jelm].nvtx-1); 
               } //TLM fix: array bounds overrun fixed here
               
               // if (i < nbrprint)  cout << " v = " << vR 
               //                         << "  " << v1 
               //                         << "  " << vL 
               //                         << "  " << v2 << endl;
               
               if (v1==vR and v2==vL) {
                  found = true;
                  
                  im = ii+1;
                  if (im > (elm[jelm].nvtx-1)) { 
                     im = im - (elm[jelm].nvtx-0); 
                  }
                  // if (i < nbrprint)  cout << "found v1==VR, v2==VL " << v1 << " " << vR << "   " << v2 << " " <<  vL << endl;
                  break; //exit edge_matching  |
               } //endif       

            } //end do     edge_matching   <---V

      //if (found) exit elms_around_vR
      if (found) {break;}

      } //end do elms_around_vR

      // Q: why is this k+2 when we already loop all the way to nvtx?
      in = k + 2; 
      if (in > elm[i].nvtx-1) { in = in - elm[i].nvtx-0; } // A: simple fix here: in > elm[i].nvtx had to be ammended and not [0,1,2](3) => 2 len=3
      // i.e. if n > 2, then c = 3; so subtract 3 (i.e. nvtx) to get back to zero
      if (found) {
         (*elm[   i].nghbr)(in) = jelm;
         (*elm[jelm].nghbr)(im) = i;
      }
      else {
         (*elm[   i].nghbr)(in) = -1; //boundary
      }


      // if (i < nbrprint) {
      //    cout << " elm nghbrs...\n";
      //    print((*elm[   i].nghbr));
      //    cout << " ------------------\n";
      // }

      }//    end do elm_vertex

   }//   end do elements2

cout << "DONE constructing the element-neighbor data " << endl;

//--------------------------------------------------------------------------------
// Edge-data for node-centered (edge-based) scheme.
//
// Loop over elements 3
// Construct edge data: edge(:).n1, n2, e1, e2.
// Edge points from node n1 to node n2.
//
//      n2
//       o------------o
//     .  \         .
//    .    \   e2  .
//   .  e1  \    .
//  .        \ .         Directed area is positive: n1 -> n2
// o----------o         e1: left element
//             n1       e2: right element (e2 > e1 or e2 = 0)

// First count the number of edges.
//
// NOTE: Count edges only if the neighbor element number is
//       greater than the current element (i) to avoid double
//       count. -1 element number indicates that it is outside
//       the domain (boundary face).

   //   elements0 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {

      v1 = (*elm[i].vtx)(0);
      v2 = (*elm[i].vtx)(1);
      v3 = (*elm[i].vtx)(2);

//    tri_quad0 : if (elm[i].nvtx==3) then
      if (elm[i].nvtx==3) {

         if ( (*elm[i].nghbr)(2) > i  or (*elm[i].nghbr)(2) == -1 ) {
            nedges = nedges + 1;
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0) == -1 ) {
            nedges = nedges + 1;
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1) == -1 ) {
            nedges = nedges + 1;
         }
      }
      
      else if (elm[i].nvtx==4) {

      v4 = (*elm[i].vtx)(3);

      if ( (*elm[i].nghbr)(2) > i or (*elm[i].nghbr)(2) == -1 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(3) > i or (*elm[i].nghbr)(3) == -1 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0) == -1 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1) == -1 ) {
       nedges = nedges + 1;
      }

      }//    endif tri_quad0

   }//   end do elements0

// Allocate the edge array.
   edge = new edge_type[nedges];
   for (size_t i = 0; i < nedges; i++) {
      edge[i].e1 = -1;
      edge[i].e2 = -1;
   }
   nedges = -1; //TLM fence post fix//

// Construct the edge data:
//  two end nodes (n1, n2), and left and right elements (e1, e2)

   int maxprint = 1;

   //elements3 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {

      v1 = (*elm[i].vtx)(0);
      v2 = (*elm[i].vtx)(1);
      v3 = (*elm[i].vtx)(2);

   
   // Triangular element
      //tri_quad2 : 
      if (elm[i].nvtx==3) {

         // if (i<maxprint) {
         //    cout << "printing edge vars to be set \n";
         //    cout << "v1 = " << v1 << "\n";
         //    cout << "v2 = " << v2 << "\n";
         //    cout << "(*elm[i].nghbr)(0) = " << (*elm[i].nghbr)(0) << "\n";
         //    cout << "(*elm[i].nghbr)(1) = " << (*elm[i].nghbr)(1) << "\n";
         //    cout << "(*elm[i].nghbr)(2) = " << (*elm[i].nghbr)(2) << "\n";
         // }

         if ( (*elm[i].nghbr)(2) > i  or (*elm[i].nghbr)(2)==-1 ) {

            // if (i<maxprint) cout << "set edge 1\n";
            nedges = nedges + 1;
            edge[nedges].n1 = v1;
            edge[nedges].n2 = v2;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(2);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -1");
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0)==-1 ) {
            // if (i<maxprint) cout << "set edge 2\n";
            nedges = nedges + 1;
            edge[nedges].n1 = v2;
            edge[nedges].n2 = v3;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(0);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -2" );
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1)==-1 ) {
            // if (i<maxprint) cout << "set edge 3\n";
            nedges = nedges + 1;
            edge[nedges].n1 = v3;
            edge[nedges].n2 = v1;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(1);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -3");
         }

         // else {
         //    cout << "ERROR: missed edge case! \n";
         // }

      }
   //  Quadrilateral element
      else if (elm[i].nvtx==4) {

         v4 = (*elm[i].vtx)(3);

         if ( (*elm[i].nghbr)(2) > i or (*elm[i].nghbr)(2) == -1 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v1;
            edge[nedges].n2 = v2;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(2);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -q1");
         }

         if ( (*elm[i].nghbr)(3) > i or (*elm[i].nghbr)(3) == -1 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v2;
            edge[nedges].n2 = v3;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(3);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -q2");
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0) == -1 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v3;
            edge[nedges].n2 = v4;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(0);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -q3");
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1) == -1 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v4;
            edge[nedges].n2 = v1;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(1);
            // assert(v1 //= v2 && "v1 should not be equal to v2 -q4");
         }

      }//    endif tri_quad2

   }//   end do elements3
   nedges += 1;

// Loop over edges
// Construct edge vector and directed area vector.
//
// Edge vector is a simple vector pointing froom n1 to n2.
// For each edge, add the directed area vector (dav) from
// the left and right elements.
//
//              n2
//   o-----------o-----------o
//   |     dav   |  dav      |
//   |       ^   |   ^       |
//   |       |   |   |       |
//   |   c - - - m - - -c    |
//   |           |           |
//   |           |           |    m: edge midpoint
//   |           |           |    c: element centroid
//   o-----------o-----------o
//                n1
//
//int pi = 0;
//   edges : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      n1 = edge[i].n1;
      n2 = edge[i].n2;
      e1 = edge[i].e1;
      e2 = edge[i].e2;
      // if (i < maxprint) {
      //    cout << "edge[i] -> i, e1, e2 " << i << " " << e1 << " " << e2 << "\n";
      // }
      // edge centroids:
      xm = half*( node[n1].x + node[n2].x );
      ym = half*( node[n1].y + node[n2].y );

      
      edge[i].dav = zero;


      // Contribution from the left element
      if (e1 > -1) {
         xc = elm[e1].x;
         yc = elm[e1].y;
         edge[i].dav(0) = -(ym-yc);
         edge[i].dav(1) =   xm-xc;
         // if (i<maxprint) {
         //    cout << "elm(e1) = " <<  elm[e1].x << " " << elm[e1].y << "\n";
         // }
      }

      // Contribution from the right element
      if (e2 > -1) {
         xc = elm[e2].x;
         yc = elm[e2].y;
         edge[i].dav(0) = edge[i].dav(0) -(yc-ym);
         edge[i].dav(1) = edge[i].dav(1) + xc-xm;
         // if (i<maxprint) {
         //    cout << "elm(e2) = " <<  elm[e2].x << " " << elm[e2].y << "\n";
         // }
      }

      if (e1 < 0 and e2 < 0) {
         cout << "ERROR: e1 and e2 are both negative... " << endl;
         cout << "n1 = " << n1 << "\n";
         cout << "n2 = " << n2 << "\n";
         cout << "e1 = " << e1 << "\n";
         cout << "e2 = " << e2 << "\n";
      }

      // Magnitude and unit vector
      edge[i].da  = std::sqrt( edge[i].dav(0) * edge[i].dav(0) + \
                                 edge[i].dav(1) * edge[i].dav(1) );
      //cout << " edge dav before division = " << edge[i].dav(0) <<  " " << edge[i].dav(1) << endl;
      edge[i].dav = edge[i].dav / edge[i].da;
      
      if (i<maxprint) {
         cout << "printing edge[i].dav \n";
         print(edge[i].dav);
         cout << " edge dav after division = " <<  edge[i].dav(0) <<  " " << edge[i].dav(1) << endl;
         if (edge[i].da < 1.e-5) {
            cout << "ERROR: collapsed edge" << endl;
            //std::exit(0);
         }
         //pi += 1;
      }
      // Edge vector
      edge[i].ev(0) = node[n2].x - node[n1].x;
      edge[i].ev(1) = node[n2].y - node[n1].y;
      edge[i].e     = std::sqrt( edge[i].ev(0) * edge[i].ev(0) + \
                                  edge[i].ev(1) * edge[i].ev(1) );
      edge[i].ev    = edge[i].ev / edge[i].e;

      
   }//   end do edges

//--------------------------------------------------------------------------------
// Construct node neighbor data:
//  pointers to the neighbor nodes(o)
//
//        o     o
//         \   / 
//          \ /
//     o-----*-----o
//          /|
//         / |
//        /  o        *: node in interest
//       o            o: neighbors (edge-connected nghbrs)
//

   cout << " --- Node-neighbor (edge connected vertex) data:" << endl;
   for (size_t i = 0; i < nnodes; i++) {
      node[i].nnghbrs = 0;
   }

   //dummy allocations are bad!:
   for (size_t i = 0; i < nedges; i++) {
      n1 = edge[i].n1;
      n2 = edge[i].n2;
      //node[n1].nghbr = new Array2D<int>(1, 1);
      //node[n2].nghbr = new Array2D<int>(1, 1);
   }

// Loop over edges and distribute the node numbers:

   //edges4 : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      n1 = edge[i].n1;
      n2 = edge[i].n2;
      

      // (1) Add n1 to the neighbor list of n2
      node[n1].nnghbrs = node[n1].nnghbrs + 1;
      if (node[n1].nnghbrs == 1) {
         node[n1].nghbr = new Array2D<int>(1, 1);
      }
      else{
         Array2D<int> save7 = (*node[n1].nghbr);
         node[n1].nghbr = new Array2D<int>(node[n1].nnghbrs, 1);
            for (size_t ra = 0; ra < save7.storage_size; ra++) {
               (*node[n1].nghbr).array[ra] = save7.array[ra];
            }
      }
      (*node[n1].nghbr)( node[n1].nnghbrs - 1 ) = n2; // more fence post trickery


      // (2) Add n2 to the neighbor list of n1
      node[n2].nnghbrs = node[n2].nnghbrs + 1;
      if (node[n2].nnghbrs == 1) {
         node[n2].nghbr = new Array2D<int>(1, 1);
      }
      else{
         Array2D<int> save8 = (*node[n2].nghbr);
         node[n2].nghbr = new Array2D<int>(node[n2].nnghbrs, 1);
            for (size_t ra = 0; ra < save8.storage_size; ra++) {
               (*node[n2].nghbr).array[ra] = save8.array[ra];
            }
      }
      (*node[n2].nghbr)( node[n2].nnghbrs - 1 ) = n1;


   } //end do edges4

//--------------------------------------------------------------------------------
// Boundary normal at nodes constructed by accumulating the contribution
// from each boundary face normal. This vector will be used to enforce
// the tangency condition, for example.
//
//
//        Interior domain      /
//                            o
//                  .        /
//                  .       /
// --o-------o-------------o
//           j   |  .  |   j+1
//               v  .  v
//
//        Left half added to the node j, and
//       right half added to the node j+1.
//

// Allocate and initialize the normal vector arrays
   for (size_t i = 0; i < nbound; i++) {

      bound[i].bnx = new Array2D<real>( bound[i].nbnodes, 1 );
      bound[i].bny = new Array2D<real>( bound[i].nbnodes, 1 );
      bound[i].bn  = new Array2D<real>( bound[i].nbnodes, 1 );

      for (size_t j = 0; j < bound[i].nbnodes; j++) {
         (*bound[i].bnx)(j) = zero;
         (*bound[i].bny)(j) = zero;
         (*bound[i].bn)( j) = zero;
      }

   }

// Normal vector at boundary nodes
// Note: Below it describes normals of linear approximation.
//       We will overwrite it by a quadratic approximation.
//
// Linear approximation:
//
// Step 1. Compute the outward normals
  for (size_t i = 0; i < nbound; i++) {
      //do j = 1, bound[i].nbnodes-1
      for (size_t j = 0; j < bound[i].nbnodes-1; j++) {

         // quick expansion for understanding:
         // int num = (*bound[i].bnode)(j);
         // x1 = node[num].x;
         x1 = node[ (*bound[i].bnode)(j) ].x;
         y1 = node[ (*bound[i].bnode)(j) ].y;

         x2 = node[ (*bound[i].bnode)(j+1) ].x;
         y2 = node[ (*bound[i].bnode)(j+1) ].y;

      //   Normal vector pointing into the domain at this point.
         (*bound[i].bnx)(j) = (*bound[i].bnx)(j) + half*( -(y2-y1) );
         (*bound[i].bny)(j) = (*bound[i].bny)(j) + half*(   x2-x1  );

         (*bound[i].bnx)(j+1) = (*bound[i].bnx)(j+1) + half*( -(y2-y1) );
         (*bound[i].bny)(j+1) = (*bound[i].bny)(j+1) + half*(   x2-x1  );

      }
   }

// Step 2. Compute the magnitude and turn (bnx,bny) into a unit vector
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {

         (*bound[i].bn)(j)  = std::sqrt( (*bound[i].bnx)(j) * (*bound[i].bnx)(j) + \
                                     (*bound[i].bny)(j) * (*bound[i].bny)(j) );
         //   Minus sign for outward pointing normal
         (*bound[i].bnx)(j) =  - (*bound[i].bnx)(j) / (*bound[i].bn)(j);
         (*bound[i].bny)(j) =  - (*bound[i].bny)(j) / (*bound[i].bn)(j);

      }
   }

// Now, ignore the linear approximation, and let us construct
// more accurate surfae normal vectors and replace the linear ones.
// So, we will overwrite the unit normal vectors: bnx, bny.
// Note: We keep the magnitude of the normal vector.
//
// Quadratic approximation:
// See http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2015v281pp518-555_preprint.pdf
// for details on the bnode(quadratic approximation for computing more accurate normals.
// 

//   boundary_type0 : do i = 1, nbound
//    boundary_nodes0 : do j = 1, bound[i].nbnodes
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {

         if (j==0) {
            v1 = (*bound[i].bnode)[j  ][0];
            v2 = (*bound[i].bnode)[j+1][0];
            v3 = (*bound[i].bnode)[j+2][0];
         }
         else if (j==bound[i].nbnodes-1) {
            v1 = (*bound[i].bnode)[j-2][0];
            v2 = (*bound[i].bnode)[j-1][0];
            v3 = (*bound[i].bnode)[j  ][0];
         }
         else {
            v1 = (*bound[i].bnode)[j-1][0];
            v2 = (*bound[i].bnode)[j  ][0];
            v3 = (*bound[i].bnode)[j+1][0];
         }

         x1 = node[v1].x;
         x2 = node[v2].x;
         x3 = node[v3].x;

         y1 = node[v1].y;
         y2 = node[v2].y;
         y3 = node[v3].y;

// //----------------------------------------------------------------------
// //   Fit a quadratic over 3 nodes

//    Skip the last one if the boundary segment is a closed boundary 
//    in which case the last node is the same as the first one.
      //if (j==bound[i].nbnodes and (*bound[i].bnode)(j)==(*bound[i].bnode)(0) ) {
      if (j==bound[i].nbnodes-1 and (*bound[i].bnode)[j][0]==(*bound[i].bnode)[0][0] ) { //TLM TODO FIX? better?
         (*bound[i].bn)(j)  = (*bound[i].bn)(0);
         (*bound[i].bnx)(j) = (*bound[i].bnx)(0);
         (*bound[i].bny)(j) = (*bound[i].bny)(0);
         continue; // cycle
      }

     dsL = std::sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
     dsR = std::sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) );
      dx = dsR*x1/(dsL*(-dsL-dsR))-x2/dsR+x2/dsL+dsL*x3/((dsR+dsL)*dsR);
      dy = dsR*y1/(dsL*(-dsL-dsR))-y2/dsR+y2/dsL+dsL*y3/((dsR+dsL)*dsR);

     ds  = std::sqrt( dx*dx + dy*dy );
     (*bound[i].bnx)(j) = -( -dy / ds );
     (*bound[i].bny)(j) = -(  dx / ds );

   }//    end do boundary_nodes0
}  //   end do boundary_type0

//--------------------------------------------------------------------------------
// Construct neighbor index over edges
//
//  Example:
//
//        o     o
//         \   / 
//          \j/       k-th neighbor
//     o-----*----------o
//          /|  edge i
//         / |
//        /  o        Note: k-th neighbor is given by "(*node(j).nghbr)(k)"
//       o
//
//  Consider the edge i
//
//   node j        k-th neighbor
//       *----------o
//      n1  edge i  n2
//
//   We store "k" in the edge data structure as
//
//    edge[i].kth_nghbr_of_1: n2 is the "edge[i].kth_nghbr_of_1"-th neighbor of n1
//    edge[i].kth_nghbr_of_2: n1 is the "edge[i].kth_nghbr_of_3"-th neighbor of n2
//
//   That is,  we have
//
//    n2 = (*node[n1].nghbr)(edge[i].kth_nghbr_of_1)
//    n1 = (*node[n2].nghbr)(edge[i].kth_nghbr_of_2)
//
//   We make use of this data structure to access off-diagonal entries in Jacobian matrix.
//

// Loop over edges

//   edges5 : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      n1 = edge[i].n1;
      n2 = edge[i].n2;

      //do k = 1, node[n2].nnghbrs
      for (size_t k = 0; k < node[n2].nnghbrs; k++) {

         if ( n1 == (*node[n2].nghbr)(k) ) {
         edge[i].kth_nghbr_of_2 = k;
         }

      }//   end do

      //do k = 1, node[n1].nnghbrs
      for (size_t k = 0; k < node[n1].nnghbrs; k++) {

         if ( n2 == (*node[n1].nghbr)(k) ) {
         edge[i].kth_nghbr_of_1 = k;
         }

      }//end do

   }//end do edges5

// Boundary mark: It should be an array actually because some nodes are associated with
//                more than one boundaries.
   for (size_t i = 0; i < nnodes; i++) {
      node[i].bmark   = -1;
      node[i].nbmarks = 0;
   }

   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {

         node[ (*bound[i].bnode)[j][0] ].bmark   = i;
         node[ (*bound[i].bnode)[j][0] ].nbmarks = node[ (*bound[i].bnode)[j][0] ].nbmarks + 1;

      }
   }

//--------------------------------------------------------------------------------
// Boundary face data
//
//      |     Domain      |
//      |                 |
//      o--o--o--o--o--o--o  <- Boundary segment
//   j= 1  2  3  4  5  6  7
//
//   In the above case, nbnodes = 7, nbfaces = 6
//

   for (size_t i = 0; i < nbound; i++) {
      bound[i].nbfaces = bound[i].nbnodes-1;
      cout << "nbfaces = " << bound[i].nbfaces << endl;

      bound[i].bfnx = new Array2D<real>( bound[i].nbfaces , 1 );
      bound[i].bfny = new Array2D<real>( bound[i].nbfaces , 1 );
      bound[i].bfn  = new Array2D<real>( bound[i].nbfaces , 1 );
      bound[i].belm = new Array2D<int>(  bound[i].nbfaces , 1 );
      bound[i].kth_nghbr_of_1 = new Array2D<int>( bound[i].nbfaces , 1 );
      bound[i].kth_nghbr_of_2 = new Array2D<int>( bound[i].nbfaces , 1 );


   }

// Boundary face vector: outward normal
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {

         x1 = node[ (*bound[i].bnode)(j  ,0) ].x;
         y1 = node[ (*bound[i].bnode)(j  ,0) ].y;
         x2 = node[ (*bound[i].bnode)(j+1,0) ].x;
         y2 = node[ (*bound[i].bnode)(j+1,0) ].y;


         if (j==0) {
            // good now
            // cout << "weird i, j  = " << i << "  "  << j << endl;
            // cout << "weird x1, y1 = " << x1 << "  "  << y1 << " x2, y2 = " << x2 << "  "  << y2 << endl;
            // cout << "bound[i].nbfaces = " << bound[i].nbfaces << endl;
            // cout << "bound[i].nbnodes = " << bound[i].nbnodes << endl;
         }

         if (j==bound[i].nbfaces-1) {
            // good now
            // cout << "weird i, j  = " << i << "  "  << j << endl;
            // cout << "weird x1, y1 = " << x1 << "  "  << y1 << " x2, y2 = " << x2 << "  "  << y2 << endl;
            // cout << "bound[i].nbfaces = " << bound[i].nbfaces << endl;
            // cout << "bound[i].nbnodes = " << bound[i].nbnodes << endl;
         }

         (*bound[i].bfn)(j,0)  =  std::sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
         (*bound[i].bfnx)(j,0) = -(y1-y2) / (*bound[i].bfn)(j);
         (*bound[i].bfny)(j) =  (x1-x2) / (*bound[i].bfn)(j);

      }
   }

// Boundary normal vector at nodes: outward normal
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {

         x1 = node[(*bound[i].bnode)[j  ][0] ].x;
         y1 = node[(*bound[i].bnode)[j  ][0] ].y;
         x2 = node[(*bound[i].bnode)[j+1][0] ].x;
         y2 = node[(*bound[i].bnode)[j+1][0] ].y;

         // if (j==bound[i].nbfaces-1){ 
         //    cout << "weird i, j  = " << i << "  "  << j << endl;
         //    cout << "weird x2, y2 = " << x2 << "  " << y2 << endl;
         //    cout << "bound[i].nbfaces = " << bound[i].nbfaces << endl;
         // }

         (*bound[i].bfn)(j)  =  std::sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
         (*bound[i].bfnx)(j) = -(y1-y2) / (*bound[i].bfn)(j);
         (*bound[i].bfny)(j) =  (x1-x2) / (*bound[i].bfn)(j);

      }
   }

// Neighbor index over boundary edges (faces)
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {

         n1 = (*bound[i].bnode)[j  ][0];  //Left node
         n2 = (*bound[i].bnode)[j+1][0];  //Right node


         for (size_t k = 0; k < node[n2].nnghbrs; k++) {
            if ( n1 == (*node[n2].nghbr)(k) ) {
               (*bound[i].kth_nghbr_of_2)(j) = k;
            }
         }

         for (size_t k = 0; k < node[n1].nnghbrs; k++) {
            if ( n2 == (*node[n1].nghbr)(k) ) {
               (*bound[i].kth_nghbr_of_1)(j) = k;
            }
         }

      }
   }

// Find element adjacent to the face: belm
//
//  NOTE: This is useful to figure out what element
//        each boundary face belongs to. Boundary flux needs
//        special weighting depending on the element.
//
//      |_________|_________|________|
//      |         |         |        | 
//      |         |         |        | 
//      |_________|_________|________|
//      |         |         |        |     <- Grid (e.g., quads)
//      |         | elmb(j) |        |
//   ---o---------o---------o--------o---  <- Boundary segment
//                 j-th face
//
// elmb(j) is the element number of the element having the j-th boundary face.
//

//   do i = 1, nbound
//    do j = 1, bound[i].nbfaces
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {

         //   bface is defined by the nodes v1 and v2.
         v1 = (*bound[i].bnode)(j) ;
         v2 = (*bound[i].bnode)(j+1);

         found = false;

      //   Find the element having the bface from the elements
      //   around the node v1.

         //do k = 1, node[v1).nelms
         for (size_t k = 0; k < node[v1].nelms; k ++) {
            //k = node[v1].nelms-1;

            ielm = (*node[v1].elm)(k);


            //cout << "v1, k, ielm  " << v1 << "    " << k << "    " << ielm << endl;
            //do ii = 1, elm[ielm].nvtx;
            for (size_t ii = 0; ii < elm[ielm].nvtx; ii++) {


               in = ii;
               im = ii+1;
               if (im > elm[ielm].nvtx-1 ) { im = im - (elm[ielm].nvtx-0); }//return to 0? (cannot use im = 0; }//)
              
               vt1 = (*elm[ielm].vtx)(in); //(in); //TLM these are bad
               vt2 = (*elm[ielm].vtx)(im); //TLM these are bad


               // vt1 = (*elm[ielm].vtx)(ii);
               // //cout << ii << endl;
               // if (ii  < elm[ielm].nvtx-1) { 
               //    vt2 = (*elm[ielm].vtx)(ii+1); 
               // }
               // else { 
               //    vt2 = (*elm[ielm].vtx)(0); 
               // } //TLM fix: array bounds overrun fixed here
               


               // if (j < 2) cout << " v = " << vt1 << "  " << v1 << "  " << vt2 << "  " << v2 << endl;
               // if (j < 2) cout << "    " << ielm << "    " << im << "    " << in << endl;
               if (vt1 == v1 and vt2 == v2) {
                  found = true;
                  //if (j < 2) cout << "found// " << endl;;//" " << in << " " << im << endl;
                  //if (j < 2) cout << "    " << ielm << "    " << im << "    " << in << endl;
                  // if (j < 2) cout << " v = " << vt1 << "  " << v1 << "  " << vt2 << "  " << v2 << endl;
                  //cout << "break 1" << endl;
                  break; //continue; //exit
               }
               // cout << "break 2" << endl;
               //if (found) {break;} //exit  //extra break needed to account for exit behavior//
            } //end do
            //cout << "break 3" << endl;
            if (found) {break;} //exit
         }//end do
            // cout << "break 3" << endl;
            // if (found) {break;} //exit

         if (found) {
            //cout << " GOOD: Boundary-adjacent element found" << endl;
            (*bound[i].belm)(j) = ielm;
         }
         else {
            cout << " Boundary-adjacent element not found. Error..." << endl;
            std::exit(0);//stop
         }

      }
   }

//--------------------------------------------------------------------------------
// Construct least-squares matrix for node-centered schemes.
//
//        o     o
//         \   / 
//          \ /
//     o-----*-----o
//          /|
//         / |
//        /  o        *: node in interest
//       o            o: neighbors (edge-connected nghbrs)
//

// Check the number of neighbor nodes (must have at least 2 neighbors)
   cout << " --- Node neighbor data:" << endl;

   ave_nghbr = node[0].nnghbrs;
   min_nghbr = node[0].nnghbrs;
   max_nghbr = node[0].nnghbrs;
      imin = 0;
      imax = 0;
   if (node[0].nnghbrs==2) {
      cout << "--- 2 neighbors for the node = " << 0 << endl;
   }

  //do i = 2, nnodes
   for (size_t i = 1; i < nnodes; i++) {
      ave_nghbr = ave_nghbr + node[i].nnghbrs;
      if (node[i].nnghbrs < min_nghbr) imin = i;
      if (node[i].nnghbrs > max_nghbr) imax = i;
      min_nghbr = std::min(min_nghbr, node[i].nnghbrs);
      max_nghbr = std::max(max_nghbr, node[i].nnghbrs);
      if (node[i].nnghbrs==2) {
         cout <<  "--- 2 neighbors for the node = " << i << endl;
      }
   }

  cout << "      nnodes    = " << nnodes    << endl;
  cout << "      ave_nghbr = " << ave_nghbr << endl;
  cout << "      ave_nghbr = " << ave_nghbr/nnodes << endl;
  cout << "      min_nghbr = " << min_nghbr << " at node " << imin << endl;
  cout << "      max_nghbr = " << max_nghbr << " at node " << imax << endl;
  cout << "" << endl;

//--------------------------------------------------------------------------------
// Cell centered scheme data
//--------------------------------------------------------------------------------

cout << "Generating CC scheme data......" << endl;

   //do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) { 
      //allocate(elm(i).edge( elm(i).nnghbrs ) )
      elm[i].edge = new Array2D<int>( elm[i].nnghbrs , 1 );
   }

   //edges3 : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {


      e1 = edge[i].e1;
      e2 = edge[i].e2;

      // Left element
      if (e1 > -1) {
         //do k = 1, elm[e1].nnghbrs;
         for (size_t k = 0; k < elm[e1].nnghbrs; k++) {
            if ( (*elm[e1].nghbr)(k)==e2) (*elm[e1].edge)(k) = i;
         }
      }

      // Right element
      if (e2 > -1) {
         //do k = 1, elm[e2].nnghbrs;
         for (size_t k = 0; k < elm[e2].nnghbrs; k++) {
            if ( (*elm[e2].nghbr)(k)==e1)  (*elm[e2].edge)(k) = i;
         }
      }

   }//end do edges3

// Face-data for cell-centered (edge-based) scheme.
//
// Loop over elements 4
// Construct face data:
// face is an edge across elements pointing
// element e1 to element e2 (e2 > e1):
//
//       e2
//        \    
//         \ face: e1 -> e2 
//          \
//  n1 o--------------o n2 <-- face
//            \
//             \          n1, n2: end nodes of the face
//              \         e1: element 1
//              e1        e2: element 2  (e2 > e1)
//
// Note: Face data is dual to the edge data.
//       It can be trivially constructed from the edge data, but
//       here the face data is constructed by using the element
//       neighbor data just for an educational purpose.

   nfaces = 0;
   //elements4
   for (size_t i = 0; i < nelms; i++) {
      for (size_t k = 0; k < elm[i].nnghbrs; k++) {
         jelm = (*elm[i].nghbr)(k);
         if (jelm > i) {
            nfaces = nfaces + 1;
         }
      }
   }//end do elements4

   //   allocate(face(nfaces))
   face = new face_type[nfaces];

   nfaces = 0;

   //   elements5
   for ( size_t i = 0; i < nelms; i++) {
      //do k = 1, elm(i).nnghbrs
      for (size_t k = 0; k < elm[i].nnghbrs; k++) {
         jelm = (*elm[i].nghbr)(k);

         if (jelm > i) {

            nfaces = nfaces + 1;

            face[nfaces-1].e1 = i;
            face[nfaces-1].e2 = jelm;

            iedge = (*elm[i].edge)(k);
            v1 = edge[iedge].n1;
            v2 = edge[iedge].n2;

            if (edge[iedge].e1 == jelm) {
               face[nfaces-1].n1 = v1;
               face[nfaces-1].n2 = v2;
            }
            else {
               face[nfaces-1].n1 = v2;
               face[nfaces-1].n2 = v1;
            }
         }
         else if (jelm == -1) {
            // if (elm[jelm].bmark != -1){
            //    cout << "ERROR: this is supposed to be a boundary \n";
            //    cout << "jelm = " << jelm << endl;
            //    cout << "elm[jelm].bmark = " << elm[jelm].bmark << "\n";
            //    cout << "-------------------------------------------"<< endl;
            //    std::exit(0); 
            // }
      //    Skip boundary faces.
         }

      }//end for
   }// elements5

// Loop over faces
// Construct directed area vector.

   //faces
   for (size_t i = 0; i < nfaces; i++) {

      n1 = face[i].n1;
      n2 = face[i].n2;
      e1 = face[i].e1;
      e2 = face[i].e2;

      // Face vector
      face[i].dav(0) = -( node[n2].y - node[n1].y );
      face[i].dav(1) =    node[n2].x - node[n1].x;
      face[i].da     = std::sqrt( face[i].dav(0)*face[i].dav(0) +
                            face[i].dav(1)*face[i].dav(1) );
      face[i].dav    = face[i].dav / face[i].da;
      if (face[i].da < 1.e-10) {
         cout << "ERROR: collapsed face" << endl;
         std::exit(0);
      }

   } //end do faces

// Construct vertex-neighbor data for cell-centered scheme.
//
// For each element, i, collect all elements sharing the nodes
// of the element, i, including face-neighors.
//
//      ___________
//     |     |     |
//     |  o  |  o  |
//     |_____|_____|
//    /\    / \    \
//   / o\ o/ i \  o \
//  /____\/_____\____\
//  \    /      /\    \
//   \o /  o   / o\ o  \
//    \/______/____\____\
//
//          i: Element of interest
//          o: Vertex neighbors (k = 1,2,...,9)

   cout << " --- Vertex-neighbor (vertex of neighbor element) data:" << endl;

   //do i = 1, nelms
   // for (size_t i = 0; i < nelms; i++) {
   //    elm[i].nvnghbrs = 1; //TLM set 0 here to speed up by avoiding loop below.
   //    //call my_alloc_int_ptr(elm(i).vnghbr, 1) //TLM TODO: redo with reallocation
   //    //Array2D<int> save4 = (*elm[i].vnghbr);
      
   //    //elm[i].vnghbr = new Array2D<int>( 1,1);

   //    // for (size_t ra; ra < save4.storage_size; ra++) {
   //    //    (*elm[i].vnghbr).array[ra] = save4.array[ra];
   //    // }
   //    //(*elm[i].vnghbr) = -2;  //intentional bad data for checking against
   // }

   ave_nghbr = 0;
   min_nghbr = 10000;
   max_nghbr =-10000;
         imin = 0;
         imax = 0;

   // Initialization
   //elements6 
   for (size_t i =0; i< nelms; i++) {
      elm[i].nvnghbrs = 0;
   }//elements6

//--------------------------------------------------------------------------------
// Collect vertex-neighbors
//
//  
   //elements7 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {


      // (1)Add face-neighbors
      //do k = 1, elm(i).nnghbrs
      for (size_t k = 0; k < elm[i].nnghbrs; k++) {
         if ( (*elm[i].nghbr)(k) > -1 ) {
            elm[i].nvnghbrs = elm[i].nvnghbrs + 1;
            if (elm[i].nvnghbrs == 1) {
               elm[i].vnghbr = new Array2D<int>( 1, 1);
            }
            else{
               Array2D<int> save5 = (*elm[i].vnghbr); //saves garbage on first step
               elm[i].vnghbr = new Array2D<int>(elm[i].nvnghbrs,1);
               for (size_t ra = 0; ra < save5.storage_size; ra++) {
                  (*elm[i].vnghbr).array[ra] = save5.array[ra]; //writes in garbage on 1st step
               }
               if (save5.storage_size == (*elm[i].vnghbr).storage_size) {
                  cout << "ERROR: at save 5 -- array size const ";
                  std::exit(0);
               }
            }
            (*elm[i].vnghbr)(elm[i].nvnghbrs-1) = (*elm[i].nghbr)(k); //eliminates garbage on 1st step
            //if (i<10*maxprint) cout << "setting " << (*elm[i].nghbr)(k)<< endl;
            //if (i<10*maxprint) cout << " to elm[i] vnghbr = " << elm[i].nvnghbrs-1 << endl;
         }
      }



      // (2)Add vertex-neighbors
      //do k = 1, elm[i].nvtx
      for (size_t k = 0; k < elm[i].nvtx; k++) {
         v1 = (*elm[i].vtx)(k);

         //velms : doj = 1, node[v1).nelms
         for (size_t j = 0; j < node[v1].nelms; j++) {

            //if (i<10*maxprint) cout << " good start HIGH after break" << endl;
            e1 = (*node[v1].elm)(j);
            if (e1 == i)  {
               //if (i<10*maxprint) cout << " e1 == i = " << i << " " << e1 <<  endl;
               continue; //velms;
            }

      //    Check if the element is already added.
            found = false;
            //do ii = 1, elm[i].nvnghbrs
            for (size_t ii = 0; ii < elm[i].nvnghbrs; ii++) {
               //if (i<10*maxprint) cout << " checking  elm " << i << "  vnghbr(  " << ii << "  ) = " <<(*elm[i].vnghbr)(ii) << endl;
               if ( e1 == (*elm[i].vnghbr)(ii) ) {
                  found = true;
                  // if (i<10*maxprint) cout << "Found element match e1 = " << e1 << " " << (*elm[i].vnghbr)(ii) << endl;
                  // if (i<10*maxprint) cout << " break" << endl;
                  break;
               }
               //if (i<10*maxprint) cout << " skip after break" << endl;
            }
            //if (i<10*maxprint) cout << " good start LOW after break" << endl;

   //       Add the element, e1, if not added yet.
            if (not found) {
               
               // if (i<10*maxprint) cout << "NO element match e1 = " << e1 << endl;
               
               elm[i].nvnghbrs = elm[i].nvnghbrs + 1;
               if (elm[i].nvnghbrs == 1) {
                  elm[i].vnghbr = new Array2D<int>( 1, 1);
               }
               else{
                  Array2D<int> save6 = (*elm[i].vnghbr);
                  elm[i].vnghbr = new Array2D<int>(elm[i].nvnghbrs, 1);
                  for (size_t ra = 0; ra < save6.storage_size; ra++) {
                     (*elm[i].vnghbr).array[ra] = save6.array[ra];
                  }
                  if (save6.storage_size == (*elm[i].vnghbr).storage_size) {
                     cout << "ERROR: at save 6 -- array size const ";
                     std::exit(0);
                  }
               }
               (*elm[i].vnghbr)(elm[i].nvnghbrs-1) = e1;
            }
         }//velms loop

      }//end elm[i].nvtx loop
      // cout << "--- elm[798].vnghbr" << endl;
      // print( (*elm[798].vnghbr) );
      // cout << "--- elm[804].vnghbr" << endl;
      // print( (*elm[804].vnghbr) );

      ave_nghbr = ave_nghbr + elm[i].nvnghbrs;
      if (elm[i].nvnghbrs < min_nghbr) imin = i;
      if (elm[i].nvnghbrs > max_nghbr) imax = i;
      min_nghbr = std::min(min_nghbr, elm[i].nvnghbrs);
      max_nghbr = std::max(max_nghbr, elm[i].nvnghbrs);
      if (elm[i].nvnghbrs < 3) {
         cout << "--- Not enough neighbors: elm = " << i << 
                  "elm[i].nvnghbrs= " << elm[i].nvnghbrs << endl;
         //std::exit(0);
      }

   }// elements7 loop
   cout << "      ave_nghbr(sum) = " << ave_nghbr << " nelms = " << nelms << endl;
   cout << "      ave_nghbr = " << ave_nghbr/nelms << endl;
   cout << "      min_nghbr = " << min_nghbr << " elm = " << imin << endl;
   cout << "      max_nghbr = " << max_nghbr << " elm = " << imax << endl;
   cout << " "  << endl;


   //do i = 1, nelms
   for ( int i = 0; i < nelms; ++i ) {
      elm[i].bmark = -1;
   }

   //bc_loop : do i = 1, nbound
   for ( int i = 0; i < nbound; ++i ) {
      //cout << bound[i].bc_type << endl;
      if ( trim( bound[i].bc_type ) == "dirichlet") {
         cout << "Found dirichlet condition " << endl;
         //do j = 1, bound[i].nbfaces
         for (size_t i = 0; i < bound[i].nbfaces; i++) {
            elm[ (*bound[i].belm)(j) ].bmark = 1;
         }//end do
      }

      // if (  trim( bound[i].bc_type ) == "freestream") {
      //    cout << "Found freestream condition " << endl;
      // }
   }// end do bc_loop

//--------------------------------------------------------------------------------



// //--------------------------------------------------------------------------------
// // Collect vertex-neighbors (alt attempt)
// //
// //  
//    //elements7 : do i = 1, nelms
//    for (size_t i = 0; i < nelms; i++) {
//    // (1)Add face-neighbors
//       //do k = 1, elm(i).nnghbrs
//       for (size_t k = 0; k < elm[i].nnghbrs; k++) {
//          if ( (*elm[i].nghbr)(k) > 0 ) {
//             elm[i].nvnghbrs = elm[i].nvnghbrs + 1;

//             //call my_alloc_int_ptr(elm[i].vnghbr, elm[i].nvnghbrs)
//             //elm[i].vnghbr = new Array2D<int>(elm[i].nvnghbrs,1);
//             //(*elm[i].vnghbr)(elm[i].nvnghbrs-1) = (*elm[i].nghbr)(k);
//          }
//       }
//    } // end elements 7 loop here

//    for (size_t i = 0; i < nelms; i++) {
//       //TLM loops to avoid realloc:
//       if ( elm[i].nvnghbrs > 0 ) {
//             elm[i].vnghbr = new Array2D<int>(elm[i].nvnghbrs,1);
//          }
//    } // end elements 7 loop here

//    //elements7
//    for (size_t i = 0; i < nelms; i++) {
//       for (size_t k = 0; k < elm[i].nnghbrs; k++) {
//          if ( (*elm[i].nghbr)(k) > 0 ) {
//             (*elm[i].vnghbr)(elm[i].tracked_index) = (*elm[i].nghbr)(k);
//             elm[i].tracked_index += 1;
//          }
//       }




//    // (2)Add vertex-neighbors
//       //do k = 1, elm[i].nvtx
//       for (size_t k = 0; k < elm[i].nvtx; k++) {
//          v1 = (*elm[i].vtx)(k);

//          //velms : do j = 1, node[v1).nelms
//          for (size_t j =0; j < node[v1].nelms; j++) {
//             e1 = (*node[v1].elm)(j);
//             if (e1 == i) continue; //velms;

//       //    Check if the element is already added.
//             found = false;
//             //do ii = 1, elm[i].nvnghbrs
//             for (size_t ii = 0; ii < elm[i].tracked_index; ii++) {
//                if ( e1 == (*elm[i].vnghbr)(ii) ) {
//                   found = true;
//                   break;
//                }
//             }

//    //       Add the element, e1, if not added yet.
//             if (not found) {
//                elm[i].tracked_index = elm[i].tracked_index + 1;
//                //call my_alloc_int_ptr(elm[i].vnghbr, elm[i].nvnghbrs)
//                elm[i].vnghbr = new Array2D<int>(elm[i].tracked_index, 1);
//                (*elm[i].vnghbr)(elm[i].tracked_index-1) = e1;
//             }
//          }//velms loop

//          }//end elm[i].nvtx loop

//       ave_nghbr = ave_nghbr + elm[i].tracked_index;
//       if (elm[i].tracked_index < min_nghbr) imin = i;
//       if (elm[i].tracked_index > max_nghbr) imax = i;
//       min_nghbr = std::min(min_nghbr, elm[i].tracked_index);
//       max_nghbr = std::max(max_nghbr, elm[i].tracked_index);
//       if (elm[i].tracked_index < 3) {
//          cout << "--- Not enough neighbors: elm = " << i << 
//                   "elm[i].nvnghbrs= " << elm[i].tracked_index << endl;
//       }

//    }// elements7 loop

//    cout << "      ave_nghbr = " << ave_nghbr/nelms << endl;
//    cout << "      min_nghbr = " << min_nghbr << " elm = " << imin << endl;
//    cout << "      max_nghbr = " << max_nghbr << " elm = " << imax << endl;
//    cout << " "  << endl;


//    //do i = 1, nelms
//    for ( int i = 0; i < nelms; ++i ) {
//       elm[i].bmark = 0;
//    }

//    //bc_loop : do i = 1, nbound
//    for ( int i = 0; i < nbound; ++i ) {
//       if ( bound[i].bc_type == "dirichlet") {
//          //do j = 1, bound[i].nbfaces
//          for (size_t i = 0; i < bound[i].nbfaces; i++) {
//             elm[ (*bound[i].belm)(j) ].bmark = 1;
//          }//end do
//       }
//    }// end do bc_loop

//--------------------------------------------------------------------------------




} //  end function construct_grid_data

//********************************************************************************






//********************************************************************************
//* Check the grid data.
//*
//* 1. Directed area must sum up to zero around every node.
//* 2. Directed area must sum up to zero over the entire grid.
//* 3. Global sum of the boundary normal vectors must vanish.
//* 4. Global sum of the boundary face normal vectors must vanish.
//* 5. Check element volumes which must be positive.
//* 6. Check dual volumes which must be positive.
//* 7. Global sum of the dual volumes must be equal to the sum of element volumes.
//*
//* Add more tests you can think of.
//*
//********************************************************************************
void EulerSolver2D::MainData2D::check_grid_data() {
//  use edu2d_my_main_data  , only : nnodes, node,  nelms,   elm, nedges, edge, &
//                             nbound, bound
//  use edu2d_constants     , only : p2, zero, half

   // //Local variables
   int i, j, n1, n2, ierr, k;
   Array2D<real>*  sum_dav_i;
   Array2D<real>  sum_dav(2,1), sum_bn(2,1);
   Array2D<real>  sum_bfn(2,1);
   real                   sum_volc, sum_vol;
   real                   mag_dav, mag_bn;
   real                   vol_min, vol_max, vol_ave;

   cout << "Checking grid data...." << endl;

   mag_dav = zero;
   mag_bn  = zero;

   sum_dav_i = new Array2D<real>(nnodes,2);

//--------------------------------------------------------------------------------
// Directed area sum check
//--------------------------------------------------------------------------------

// Compute the sum of the directed area for each node.

   (*sum_dav_i) = zero;
   //print( (*sum_dav_i) );
   for (size_t i = 0; i < nedges; i++) {
      n1 = edge[i].n1;
      n2 = edge[i].n2;
      (*sum_dav_i)(n1,0) = (*sum_dav_i)(n1,0) + edge[i].dav(0)*edge[i].da;
      (*sum_dav_i)(n1,1) = (*sum_dav_i)(n1,1) + edge[i].dav(1)*edge[i].da;
      (*sum_dav_i)(n2,0) = (*sum_dav_i)(n2,0) - edge[i].dav(0)*edge[i].da;
      (*sum_dav_i)(n2,1) = (*sum_dav_i)(n2,1) - edge[i].dav(1)*edge[i].da;
      mag_dav = mag_dav + edge[i].da;
      // cout << (*sum_dav_i)(n1,0) << endl;
      // cout << (*sum_dav_i)(n2,0) << endl;
      // cout << (*sum_dav_i)(n1,1) << endl;
      // cout << (*sum_dav_i)(n2,1) << endl;

      // cout << edge[i].dav(0) << endl;
      // cout << edge[i].dav(1) << endl;
      // cout << " edge da = " << edge[i].da << endl;
   }
   //ok:
   cout << "sum edge[i].da =" << mag_dav << endl;
   mag_dav = mag_dav/real(nedges);
   cout << "sum dav i, mag = " << mag_dav << endl;

// Add contribution from boundary edges.
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {

      n1 = (*bound[i].bnode)(j);
      n2 = (*bound[i].bnode)(j+1);

      (*sum_dav_i)(n1,0) = (*sum_dav_i)(n1,0) + half*(*bound[i].bfnx)(j)*(*bound[i].bfn)(j);
      (*sum_dav_i)(n1,1) = (*sum_dav_i)(n1,1) + half*(*bound[i].bfny)(j)*(*bound[i].bfn)(j);

      (*sum_dav_i)(n2,0) = (*sum_dav_i)(n2,0) + half*(*bound[i].bfnx)(j)*(*bound[i].bfn)(j);
      (*sum_dav_i)(n2,1) = (*sum_dav_i)(n2,1) + half*(*bound[i].bfny)(j)*(*bound[i].bfn)(j);

      }
   }

// Compute also the sum of the boundary normal vector (at nodes).

   sum_bn = 0;
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {
         k = (*bound[i].bnode)(j);
         if (j > 0 and k==(*bound[i].bnode)(0)) continue; //Skip if the last node is equal to the first node).
         sum_bn(0)      = sum_bn(0)      + (*bound[i].bnx)(j)*(*bound[i].bn)(j);
         sum_bn(1)      = sum_bn(1)      + (*bound[i].bny)(j)*(*bound[i].bn)(j);
         mag_bn = mag_bn + abs((*bound[i].bn)(j));
      }
      mag_bn = mag_bn/real(bound[i].nbnodes);//
   }

// Global sum of boundary normal vectors must vanish.

 if (sum_bn(0) > 1.0e-12*mag_bn and sum_bn(1) > 1.0e-12*mag_bn) {

   cout << "--- Global sum of the boundary normal vector:" << endl;
   cout << "    sum of bn_x = " << sum_bn(0)  << endl;
   cout << "    sum of bn_y = " << sum_bn(1)  << endl;
   cout << "Error: boundary normal vectors do not sum to zero..." << endl;
   std::exit(0);//stop program
}

// Sum of the directed area vectors must vanish at every node.

   for (size_t i = 0; i < nnodes; i++) {
         if ( std::abs( (*sum_dav_i)(i,0) ) > 1.0e-4 * mag_dav or 
               std::abs( (*sum_dav_i)(i,1) ) > 1.0e-4 * mag_dav)
         {
            cout << "Error: Sum of the directed area vectors large sum_dav...\n";
            cout  << " --- node=" << i << "\n"
                  << " (x,y)=" << node[i].x << " , " << node[i].y << "\n"
                  << " sum_dav=" << (*sum_dav_i)(i,0) << " , " << (*sum_dav_i)(i,1) << endl;
            cout << "-------------------------------------------------------" << endl;
            //std::exit(0);//stop program
         }
   }
   //print( (*sum_dav_i) );

   cout << "--- Max sum of directed area vector around a node:" << endl;
   // cout << "  max(sum_dav_i_x) = " <<  maxval((*sum_dav_i)(:,0)) << endl;
   // cout << "  max(sum_dav_i_y) = " <<  maxval((*sum_dav_i)(:,1)) << endl;
   cout << "  max(sum_dav_i_x) = " <<  MaxColVal( (*sum_dav_i),0) << endl;
   cout << "  max(sum_dav_i_y) = " <<  MaxColVal( (*sum_dav_i),1) << endl;




   if (MaxColVal( abs( (*sum_dav_i) ) , 0 ) > 1.0e-12 * mag_dav or
      MaxColVal( abs( (*sum_dav_i) ) , 1 ) > 1.0e-12 * mag_dav)   
   {
      cout << "--- Max sum of directed area vector around a node:" << endl;
      cout << "  max(sum_dav_i_x) = " <<  MaxColVal( (*sum_dav_i), 0) << endl;
      cout << "  max(sum_dav_i_y) = " <<  MaxColVal( (*sum_dav_i), 1) << endl;
      cout << "Error: directed area vectors do not sum to zero..." << endl;
      std::exit(0);//stop
   }

// Of course, the global sum of the directed area vector sum must vanish.
   sum_dav = zero;
   cout << "sum_dav 0 = " << sum_dav(0) << endl;
   cout << "sum_dav 1 = " << sum_dav(1) << endl;
   for (size_t i = 0; i < nnodes; i++) {
      sum_dav(0) = sum_dav(0) + (*sum_dav_i)(i,0) ;
      sum_dav(1) = sum_dav(1) + (*sum_dav_i)(i,1) ;
      // cout << "sum_dav i,0 =" << (*sum_dav_i)(i,0)  << endl;
      // cout << "sum_dav i,1 =" << (*sum_dav_i)(i,1) << endl;
      //cout << "sum_dav i =" << (*sum_dav_i)(i,0) << sum_dav(0) << endl;
      //cout << "sum_dav i =" << (*sum_dav_i)(i,1) << sum_dav(1) << endl;
      //cout << " add them = " << sum_dav(0) + (*sum_dav_i)(i,0) << endl;
   }

   cout << "--- Global sum of the directed area vector:" << endl;
   cout << "    sum of dav_x = " <<  sum_dav(0) << endl;
   cout << "    sum of dav_y = " <<  sum_dav(1) << endl;

  if (sum_dav(1) > 1.0e-12 *mag_dav and sum_dav(2) > 1.0e-12 *mag_dav)   {
   cout << "Error: directed area vectors do not sum globally to zero..." << endl;
   cout << "--- Global sum of the directed area vector:" << endl;
   cout << "    sum of dav_x = " <<  sum_dav(0) << endl;
   cout << "    sum of dav_y = " <<  sum_dav(1) << endl;
   std::exit(0);//stop
  }



//--------------------------------------------------------------------------------
// Global sum check for boundary face vector
//--------------------------------------------------------------------------------
   sum_bfn = 0;
   // do i = 1, nbound
   //    do j = 1, bound[i].nbfaces
   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbfaces; j++) {
         sum_bfn(0) =  sum_bfn(0) + (*bound[i].bfnx)(j)*(*bound[i].bfn)(j);
         sum_bfn(1) =  sum_bfn(1) + (*bound[i].bfny)(j)*(*bound[i].bfn)(j);
      }
   }

   cout << "--- Global sum of the boundary face vector:" << endl;
   cout << "    sum of bfn_x = " <<  sum_bfn(0) << endl;
   cout << "    sum of bfn_y = " <<  sum_bfn(1) << endl;

  if (sum_bfn(1) > 1.0e-12 *mag_bn and sum_bfn(2) > 1.0e-12 *mag_bn)   {
   cout << "Error: boundary face normals do not sum globally to zero..." << endl;
   cout << "--- Global sum of the boundary face normal vector:" << endl;
   cout << "    sum of bfn_x = " <<  sum_bfn(0) << endl;
   cout << "    sum of bfn_y = " <<  sum_bfn(1) << endl;
   std::exit(0);//stop
  }

//--------------------------------------------------------------------------------
// Volume check
//--------------------------------------------------------------------------------
// (1)Check the element volume: make sure there are no zero or negative volumes

   vol_min =  1.0e+15;
   vol_max = -1.0;
   vol_ave =  zero;

      ierr = 0;
   sum_volc = zero;
   for ( int i = 0; i < nelms; ++i ) {

      vol_min = std::min(vol_min, elm[i].vol);
      vol_max = std::max(vol_max, elm[i].vol);
      vol_ave = vol_ave + elm[i].vol;

   sum_volc = sum_volc + elm[i].vol;

   if (elm[i].vol < zero)   {
     cout << "Negative volc=" << elm[i].vol <<   "elm= " << i <<  " stop..." << endl;
     ierr = ierr + 1;
   }

   if ( std::abs(elm[i].vol) < 1.0e-14 )   {
     cout << "Vanishing volc= "  << elm[i].vol << " " 
                                 << std::abs(elm[i].vol) 
                                 <<   " elm= " << i 
                                 << " stop..." << endl;
     cout << " val = " << abs(elm[i].vol) << " tol = " << 1.0e-14 << endl;
     ierr = ierr + 1;
   }

  }

   vol_ave = vol_ave / real(nelms);

   cout << " " << endl;
   cout << "    minimum element volume = " <<  vol_min << endl;
   cout << "    maximum element volume = " <<  vol_max << endl;
   cout << "    average element volume = " <<  vol_ave << endl;
   cout << " " << endl;

//--------------------------------------------------------------------------------
// (2)Check the dual volume (volume around a node)

   vol_min =  1.0e+15;
   vol_max = -1.0;
   vol_ave =  zero;

      ierr = 0;
   sum_vol = zero;
   for (size_t i = 0; i < nnodes; i++) {

      vol_min = std::min(vol_min,node[i].vol);
      vol_max = std::max(vol_max,node[i].vol);
      vol_ave = vol_ave + node[i].vol;

      sum_vol = sum_vol + node[i].vol;

      if (node[i].vol < zero)   {
         cout << "Negative vol=" << node[i].vol <<   node << " << i <<   stop..." << endl;
         ierr = ierr + 1;
      }

      if (std::abs(node[i].vol) < 1.0e-14 )   {
         cout << "Vanishing vol=" << node[i].vol <<   node << " << i <<   stop..." << endl;
         ierr = ierr + 1;
      }
   }

   vol_ave = vol_ave / real(nnodes);

   cout << " " << endl;
   cout << "    minimum dual volume = " <<  vol_min << endl;
   cout << "    maximum dual volume = " <<  vol_max << endl;
   cout << "    average dual volume = " <<  vol_ave << endl;
   cout << " " << endl;


   if (ierr > 0) std::exit(0);//stop

   if (abs(sum_vol-sum_volc) > 1.0e-08 *sum_vol)   {
      cout << "--- Global sum of volume: must be the same" << endl;
      cout << "    sum of volc = " <<  sum_volc << endl;
      cout << "    sum of vol  = " <<  sum_vol << endl;
      cout << " sum_vol-sum_volc  = " <<  sum_vol-sum_volc << endl;
      cout << "Error: sum of dual volumes and cell volumes do not match..." << endl;
      std::exit(0);//stop
   }

   check_skewness_nc();
   compute_ar();

   cout << " " << endl;
   cout << "Grid data look good\n" << endl;

   cout << " " << endl;
   cout << "    nfaces = " <<  nfaces << endl;
   cout << "    nedges = " <<  nedges << endl;
   cout << "    nbound = " <<  nbound << endl;
   cout << "    nelms = " <<  nelms << endl;
   cout << " " << endl;
} //end  check_grid_data




//*******************************************************************************
//* Skewness computation for edges.
//*******************************************************************************
void EulerSolver2D::MainData2D::check_skewness_nc() {

   int i;
   real e_dot_n, e_dot_n_min, e_dot_n_max, alpha;

      e_dot_n = zero;
   e_dot_n_min = 100000.0;
   e_dot_n_max =-100000.0;

   //edges : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      alpha = edge[i].ev(0) * edge[i].dav(0) + edge[i].ev(1) * edge[i].dav(1);
      e_dot_n     = e_dot_n + std::abs(alpha);
      e_dot_n_min = std::min(e_dot_n_min, std::abs(alpha) );
      e_dot_n_max = std::max(e_dot_n_max, std::abs(alpha) );

   } //end do edges

   e_dot_n = e_dot_n / real(nedges);

   cout << " " <<  "\n";
   cout << " ------ Skewness check (NC control volume) ----------\n";
   cout << "   L1(e_dot_n) = " << e_dot_n << "\n";
   cout << "  Min(e_dot_n) = " << e_dot_n_min << "\n";
   cout << "  Max(e_dot_n) = " << e_dot_n_max << "\n";
   cout << " ----------------------------------------------------" << endl;

 }//end subroutine check_skewness_nc


// //*******************************************************************************
// //* Control volume aspect ratio
// //*******************************************************************************
void EulerSolver2D::MainData2D::compute_ar() {
//  use edu2d_my_main_data  , only : node, nnodes, elm, nelms
//  use edu2d_constants     , only : p2, zero, one, half, two


   int i, n1, n2;
   real ar, ar_min, ar_max, nnodes_eff;

   int k;
   real side_max, side_mid, side_min, height;
   Array2D<real> side(4,1);

   // Initialization

   //node1 : do i = 1, nnodes
   for (size_t i = 0; i < nnodes; i++) {
      node[i].ar     = zero;
   } //node1

// Compute element aspect-ratio: longest_side^2 / vol

   //elm1: do i = 1, nelms
   for (size_t i =0; i < nelms; i++) {

      side_max = -one;
   
      //do k = 1, elm[i].nvtx
      for (size_t k = 0; k < elm[i].nvtx; k ++) {

         n1 = (*elm[i].vtx)(k);
         if (k == elm[i].nvtx-1) {
            n2 = (*elm[i].vtx)(0);
         }
         else {
            n2 = (*elm[i].vtx)(k+1);
         }

         side(k) = std::sqrt( (node[n2].x-node[n1].x) * (node[n2].x-node[n1].x) \
                           + (node[n2].y-node[n1].y) * (node[n2].y-node[n1].y) );
         side_max =  std::max(side_max, side(k));

      }//end do

      if (elm[i].nvtx == 3) {

   // AR for triangle:  Ratio of a half of a square with side_max to volume
         elm[i].ar = (half*side_max*side_max) / elm[i].vol;

         if (side(0) >= side(1) and side(0) >= side(2)) {

            side_max = side(0);
            if (side(1) >= side(2)) {
               side_mid = side(1); 
               side_min = side(2);
            }
            else {
               side_mid = side(2); 
               side_min = side(1);
            }
         }
         else if (side(1) >= side(0) and side(1) >= side(2)) {

            side_max = side(1);
            if (side(0) >= side(2)) {
               side_mid = side(0); 
               side_min = side(2);
            }
            else {
               side_mid = side(2); 
               side_min = side(0);
            }
         }
         else {

            side_max = side(2);
            if (side(0) >= side(1)) {
               side_mid = side(0);
               side_min = side(1);
            } 
            else {
               side_mid = side(1);
               side_min = side(0);
            }

         }

         height = two*elm[i].vol/side_mid;
         elm[i].ar = side_mid/height;

      }

      else {

   // AR for quad: Ratio of a square with side_max to volume
      elm[i].ar = side_max*side_max / elm[i].vol;

      } //endif

   }//end do elm1

   // Compute the aspect ratio:
   //node2 : do i = 1, nnodes
   for (size_t i = 0; i < nnodes; i++) {

      node[i].ar = zero;
      //do k = 1, node[i].nelms
      for (size_t k = 0; k < node[i].nelms; k++) {
         node[i].ar = node[i].ar + elm[ (*node[i].elm)(k) ].ar;
      } //end do

      node[i].ar = node[i].ar / real(node[i].nelms);

   } //end do node2;

// // Compute the min/max and L1 of AR

   nnodes_eff= zero;
         ar = zero;
      ar_min = 100000.0;
      ar_max =-100000.0;

   //node3: do i = 1, nnodes
   for (size_t i = 0; i < nnodes; i++) {
      if (node[i].bmark != -1) { continue; }//cycle node3 
      ar     = ar + std::abs(node[i].ar);
      ar_min = std::min(ar_min, std::abs(node[i].ar) );
      ar_max = std::max(ar_max, std::abs(node[i].ar) );
      nnodes_eff = nnodes_eff + one;
   }//end do node3

   ar = ar / nnodes_eff;

   cout << " \n";
   cout << " ------ Aspect ratio check (NC control volume) ----------\n";
   cout << " Interior nodes only \n";
   cout << "   L1(AR) = " << ar << "\n";
   cout << "  Min(AR) = " << ar_min << "\n";
   cout << "  Max(AR) = " << ar_max << endl;

   nnodes_eff= zero;
         ar = zero;
      ar_min = 100000.0;
      ar_max =-100000.0;

   //node4: do i = 1, nnodes
   for (size_t i = 0; i < nnodes; i++) {
      if (node[i].bmark == -1) {continue;} //cycle node4
      ar     = ar + std::abs(node[i].ar);
      ar_min = std::min(ar_min, std::abs(node[i].ar) );
      ar_max = std::max(ar_max, std::abs(node[i].ar) );
      nnodes_eff = nnodes_eff + one;
   }//end do node4

   ar = ar / nnodes_eff;

   cout << " " << "\n";
   cout << " Boundary nodes only" << "\n";
   cout << "   L1(AR) = " << ar << "\n";
   cout << "  min(AR) = " << ar_min << "\n";
   cout << "  Max(AR) = " << ar_max << "\n";
   cout << " --------------------------------------------------------" << endl;

} // compute_ar



//
//*******************************
// plot grids
void EulerSolver2D::MainData2D::write_tecplot_file(const std::string& datafile) {

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
      outfile  << node[i].x << '\t' 
               << node[i].y << "\n"; 
   }
   //--------------------------------------------------------------------------------
   for ( int i = 0; i < nelms; ++i ) {
      //Triangles
      if (elm[i].nvtx==3) {
         outfile  << (*elm[i].vtx)(0,0) << '\t' 
                  << (*elm[i].vtx)(1,0) << '\t' 
                  << (*elm[i].vtx)(2,0) << '\t' 
                  << (*elm[i].vtx)(2,0) <<  "\n"; //The last one is a dummy.
      }

      //Quadrilaterals
      else if (elm[i].nvtx==4) {
         outfile  << (*elm[i].vtx)(0,0) << '\t' 
                  << (*elm[i].vtx)(1,0) << '\t' 
                  << (*elm[i].vtx)(2,0) << '\t' 
                  << (*elm[i].vtx)(3,0) <<  "\n"; //The last one is a dummy.

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
void EulerSolver2D::MainData2D::write_grid_file(const std::string& datafile) {//(char* datafile)
    int i,j,os;
//--------------------------------------------------------------------------------
    ofstream outfile;
    outfile.open (datafile);

//--------------------------------------------------------------------------------
// Grid size: # of nodes, # of triangles, # of quadrilaterals
    outfile <<  nnodes << '\t' << ntria << '\t' << nquad << "\n";

//--------------------------------------------------------------------------------
// Node data
   for (int i=0; i<nnodes; ++i) {
   outfile  << node[i].x << '\t' 
            << node[i].y << "\n";
   }

//--------------------------------------------------------------------------------
   for ( int i = 0; i < nelms; ++i ) {
      //Triangles
      if (elm[i].nvtx==3) {
         outfile  << (*elm[i].vtx)(0,0) << '\t' 
                  << (*elm[i].vtx)(1,0) << '\t' 
                  //<< (*elm[i].vtx)(2,0) << '\t' 
                  << (*elm[i].vtx)(2,0) <<  "\n"; //The last one is a dummy.
      }

      //Quadrilaterals
      else if (elm[i].nvtx==4) {
         outfile  << (*elm[i].vtx)(0,0) << '\t' 
                  << (*elm[i].vtx)(1,0) << '\t' 
                  << (*elm[i].vtx)(2,0) << '\t' 
                  << (*elm[i].vtx)(3,0) <<  "\n"; //The last one is a dummy.

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
   outfile <<  nbound << "\n"; //5



   //  outfile <<  (ny-1)/2+1  << "\n"; //Inflow
   //  outfile <<  (ny-1)/2+1  << "\n"; //Left Wall
   //  outfile <<   nx << "\n";         //Bottom Outflow
   //  outfile <<   ny  << "\n";        //Right  Outflow
   //  outfile <<   nx  << "\n";        //Top Wall

   //  outfile << "\n";


//    std::cout << " Boundary conditions:" << std::endl;
   for (size_t i = 0; i < nbound; i++) {
      //std::cout << " boundary" << i << "  bc_type = " << bound[i].bc_type << std::endl;
      outfile << bound[i].nbnodes << "\n";
   }

    outfile << "\n";

   for (size_t i = 0; i < nbound; i++) {
      for (size_t j = 0; j < bound[i].nbnodes; j++) {
         outfile << (*bound[i].bnode)(j,0) << "\n";
      }
   }
// // Inflow boundary
//     //do j = ny, (ny-1)/2+1, -1
//     for (int j=ny-1; j<(ny-1)/2+1; ++j) {
//         i = 1;
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

// // // Left wall boundary
//     //do j = (ny-1)/2+1, 1, -1
//     for (int j=(ny-1)/2+1; j>=1; --j) {
//         i = 1;
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

// // // Bottom outflow boundary
// //     do i = 1, nx
//     for (int i=nx-1; i<nx; ++i) {
//         j = 1;
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

// // // Right outflow boundary
// //     do j = 1, ny
//     for (int j=1; j<ny; ++j) {
//         i = nx;
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

// // // Top wall boundary
// //     do i = nx, 1, -1
//     for (int i=0; i<nx; ++i) {
//         j = ny;
//         outfile <<  i + (j-1)*nx  << "\n";  
//     }

//--------------------------------------------------------------------------------
    outfile.close();
}
//********************************************************************************
