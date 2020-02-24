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
         // Fix indices for 0 indexed code!
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
         // Fix indices for 0 indexed code!
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
         //  were starting from 1, but our code starts form 0!
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
   vL = 0;
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
   cout << "nelms = " << nnodes << "\n" << endl;

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

   for ( int i = 0; i < nelms; ++i ) {
      v1 = (*elm[i].vtx)(0,0);
      v2 = (*elm[i].vtx)(1,0);
      v3 = (*elm[i].vtx)(2,0);
      node[v1].nelms = node[v1].nelms + 1;
      node[v2].nelms = node[v2].nelms + 1;
      node[v3].nelms = node[v3].nelms + 1;

   }
   for ( int i = 0; i < nelms; ++i ) {
      v1 = (*elm[i].vtx)(0,0);
      v2 = (*elm[i].vtx)(1,0);
      v3 = (*elm[i].vtx)(2,0);
      node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
      node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
      node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
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
      * DESIGN CHOICE HERE -- use std<vector?> or similar?
      * I don't think we need a true dynamic (i.e. resizable) array
      */

      // node[v1].nelms = node[v1].nelms + 1;
      //attempt copy constructor //Array2D<int> save1( (*node[v1].elm) );
      // node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
      // (*node[v1].elm)(node[v1].nelms-1, 0) = i;

      (*node[v1].elm)((*node[v1].elm).tracked_index, 0) = i;
      (*node[v1].elm).tracked_index +=1;

      // node[v1].nelms = node[v1].nelms + 1;
      // //Array2D<int> save1 = Array2D<int>(node[v1].nelms, 1);
      // // invoke the copy constructor:
      // Array2D<int> save1( (*node[v1].elm) );
      // //save1.print();
      // //delete [] node[v1].elm;
      // node[v1].elm = new Array2D<int>(node[v1].nelms, 1);
      // for (size_t ra = 0; ra < save1.storage_size; ra++){
      //    (*node[v1].elm)(ra) = save1(ra);
      // }
      // (*node[v1].elm)(node[v1].nelms-1, 0) = i;

      //(*node[v1].elm)(node[v1].nelms-1, 0) = i;
      // if ( 80190 < i < 80210) cout << "index to node: " << v1 << " " <<  
      //                   (*node[v1].elm)[node[v1].nelms-1][0] << endl;

      if (i == 160400) {
         cout << "----------------------------------" << endl;
         cout << (*node[v1].elm)(node[v1].nelms-1) << " " << i << endl;
         cout << "v1 = " << v1 << ", node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         ielm = (*node[v1].elm)(node[v1].nelms-1, 0);
         cout << "ielm = " << ielm << endl;
         cout << "node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         cout << "----------------------------------" << endl;
      }
      if (i == 319201) {
         cout << "----------------------------------" << endl;
         cout << (*node[v1].elm)(node[v1].nelms-1) << " " << i << endl;
         cout << "v1 = " << v1 << ", node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         ielm = (*node[v1].elm)(node[v1].nelms-1, 0);
         cout << "ielm = " << ielm << endl;
         cout << "node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         cout << "----------------------------------" << endl;
      }

      if (i == 80200) {
         cout << "----------------------------------" << endl;
         cout << (*node[v1].elm)(node[v1].nelms-1) << " " << i << endl;
         cout << "v1 = " << v1 << ", node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         ielm = (*node[v1].elm)(node[v1].nelms-1, 0);
         cout << "ielm = " << ielm << endl;
         cout << "node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         cout << "----------------------------------" << endl;
      }


      if (i == 159201) {
         cout << "----------------------------------" << endl;
         cout << (*node[v1].elm)(node[v1].nelms-1) << " " << i << endl;
         cout << "v1 = " << v1 << ", node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         ielm = (*node[v1].elm)(node[v1].nelms-1, 0);
         cout << "ielm = " << ielm << endl; 
         cout << "node[v1].nelms-1 = " << node[v1].nelms-1 << endl;
         cout << "----------------------------------" << endl;
      }

      // node[v2].nelms = node[v2].nelms + 1;
      // // node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
      // //
      // // invoke the copy constructor:
      // Array2D<int> save2( (*node[v2].elm) );
      // // call resize copy?
      // // node[v1].elm = new Array2D<int>(*node[v1].elm, node[v1].nelms, 0);
      
      // delete [] node[v2].elm;
      // node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
      // for (size_t ra = 0; ra < save2.storage_size; ra++){
      //    (*node[v2].elm)(ra) = save2(ra);
      // }
      // (*node[v2].elm)(node[v2].nelms-1 , 0) = i;


      // node[v3].nelms = node[v3].nelms + 1;
      // // // node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
      // // Array2D<int> save3 = Array2D<int>(node[v3].nelms, 1);
      // // invoke the copy constructor:
      // Array2D<int> save3( (*node[v3].elm) );
      // delete [] node[v3].elm;
      // node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
      // for (size_t ra = 0; ra < save3.storage_size; ra++){
      //    (*node[v3].elm)(ra) = save3(ra);
      // }
      // (*node[v3].elm)(node[v3].nelms-1 , 0) = i;

      // simplest, but ineffective:
      // node[v2].nelms = node[v2].nelms + 1;
      // node[v2].elm = new Array2D<int>(node[v2].nelms, 1);
      // (*node[v2].elm)[node[v2].nelms-1][0] = i;

      // node[v3].nelms = node[v3].nelms + 1;
      // node[v3].elm = new Array2D<int>(node[v3].nelms, 1);
      // (*node[v3].elm)[node[v3].nelms-1][0] = i;


      (*node[v2].elm)((*node[v2].elm).tracked_index, 0) = i;
      (*node[v2].elm).tracked_index +=1;

      (*node[v3].elm)((*node[v3].elm).tracked_index, 0) = i;
      (*node[v3].elm).tracked_index +=1;


         // for (size_t jw = 0; j < node[v1].nelms; jw++) {
         //    jelm = (*node[v1].elm)(jw);
            
         //    //tlm issue seen here:
         //    cout << "jelm = " << jelm << endl;
         // }

      // node[v2].elm = new Array2D<int>(node[v2].nelms+1, 1);
      // (*node[v2].elm)[node[v2].nelms][0] = i;
      // node[v2].nelms = node[v2].nelms + 1;

      // node[v3].elm = new Array2D<int>(node[v3].nelms+1, 1);
      // (*node[v3].elm)[node[v3].nelms][0] = i;
      // node[v3].nelms = node[v3].nelms + 1;

      //FIX here!
      //(*node[v4].elm)[ node[v4].nelms ][0] = i;

      // Compute the cell center and cell volume.
      //tri_or_quad : if (elm(i).nvtx==3) then
      if (elm[i].nvtx==3) {

      // Triangle centroid and volume
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
         if (tri_area(x1,x2,xc,y1,y2,yc)<=zero) {
            cout << " Centroid outside the quad element 12c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            //stop
         }

         if (tri_area(x2,x3,xc,y2,y3,yc)<=zero) {
            cout << " Centroid outside the quad element 23c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            //stop
         }

         if (tri_area(x3,x4,xc,y3,y4,yc)<=zero) {
            cout << " Centroid outside the quad element 34c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1 << endl;
            cout << "  (x2,y2)=" << x2 << y2 << endl;
            cout << "  (x3,y3)=" << x3 << y3 << endl;
            cout << "  (x4,y4)=" << x4 << y4 << endl;
            cout << "  (xc,yc)=" << xc << yc << endl;
            //stop
         }

         if (tri_area(x4,x1,xc,y4,y1,yc)<=zero) {
            cout << " Centroid outside the quad element 41c: i=" << i << endl;
            cout << "  (x1,y1)=" << x1 << y1  << endl;
            cout << "  (x2,y2)=" << x2 << y2  << endl;
            cout << "  (x3,y3)=" << x3 << y3  << endl;
            cout << "  (x4,y4)=" << x4 << y4  << endl;
            cout << "  (xc,yc)=" << xc << yc  << endl;
            //stop
         }

      //  Distribution of element number to the 4th node of the quadrilateral
         // node[v4].nelms = node[v4].nelms + 1;
         // node[v4].elm = new Array2D<int>(node[v4].nelms, 1);
         // (*node[v4].elm)[ node[v4].nelms-1 ][0] = i;

         // easier way
         (*node[v4].elm)((*node[v4].elm).tracked_index, 0) = i;
         (*node[v4].elm).tracked_index +=1;

         // node[v4].elm = new Array2D<int>(node[v4].nelms+1, 1);
         // (*node[v4].elm)( node[v4].nelms , 0) = i;
         // node[v4].nelms = node[v4].nelms + 1;

      }//    endif tri_or_quad

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

   }

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
      //elms_around_vR : do j = 1, node(vR).nelms;
         for (size_t j = 0; j < node[vR].nelms; j++) {
            jelm = (*node[vR].elm)(j);
            
            //tlm issue seen here:
            //cout << "jelm = " << jelm << endl;
            
            // cout << vR << " " << j << "   " << (*node[vR].elm)(j) << endl;
            //if (i < 20) cout << "vR , jelm = "<< vR << "   " << jelm << endl;

            //edge_matching : do ii = 1, elm(jelm).nvtx;
            for (size_t ii = 0; ii < elm[jelm].nvtx; ii++) {
               
               v1 = (*elm[jelm].vtx)(ii);
               //cout << ii << endl;
               if (ii  > 0) { 
                  v2 = (*elm[jelm].vtx)(ii-1); 
               }
               else if (ii == 0) { 
                  v2 = (*elm[jelm].vtx)(elm[jelm].nvtx-1); 
               } //TLM fix: array bounds overrun fixed here
               
               if (v1==vR and v2==vL) {
                  found = true;
                  //if (k < 2) cout << " v = " << vR << "  " << v1 << "  " << vL << "  " << v2 << endl;
               
                  //cout << "found v1==VR, v2==VL " << v1 << " " << vR << "   " << v2 << " " <<  vL << endl;
                  im = ii+1;
                  if (im > (elm[jelm].nvtx-1)) { 
                     im = im - (elm[jelm].nvtx-0); 
                  }
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
         //(*elm[   i].nghbr)(in) = 0;
         (*elm[   i].nghbr)(in) = 0;
      }

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
//       count. Zero element number indicates that it is outside
//       the domain (boundary face).

   //   elements0 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {

      v1 = (*elm[i].vtx)(0);
      v2 = (*elm[i].vtx)(1);
      v3 = (*elm[i].vtx)(2);

//    tri_quad0 : if (elm[i].nvtx==3) then
      if (elm[i].nvtx==3) {

         if ( (*elm[i].nghbr)(2) > i  or (*elm[i].nghbr)(2)==0 ) {
            nedges = nedges + 1;
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0)==0 ) {
            nedges = nedges + 1;
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1)==0 ) {
            nedges = nedges + 1;
         }
      }
      else if (elm[i].nvtx==4) {

      v4 = (*elm[i].vtx)(3);

      if ( (*elm[i].nghbr)(2) > i or (*elm[i].nghbr)(2) ==0 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(3) > i or (*elm[i].nghbr)(3) ==0 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0) ==0 ) {
         nedges = nedges + 1;
      }

      if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1) ==0 ) {
       nedges = nedges + 1;
      }

      }//    endif tri_quad0

   }//   end do elements0

// Allocate the edge array.
   edge = new edge_type[nedges];
   for (size_t i = 0; i < nedges; i++) {
      edge[i].e1 = 0;
      edge[i].e2 = 0;
   }
   nedges = -1; //TLM fence post fix!

// Construct the edge data:
//  two end nodes (n1, n2), and left and right elements (e1, e2)

   //elements3 : do i = 1, nelms
   for (size_t i = 0; i < nelms; i++) {

      v1 = (*elm[i].vtx)(0);
      v2 = (*elm[i].vtx)(1);
      v3 = (*elm[i].vtx)(2);

   
   // Triangular element
      //tri_quad2 : 
      if (elm[i].nvtx==3) {

         if ( (*elm[i].nghbr)(2) > i  or (*elm[i].nghbr)(2)==0 ) {

            nedges = nedges + 1;
            edge[nedges].n1 = v1;
            edge[nedges].n2 = v2;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(2);
            // assert(v1 != v2 && "v1 should not be equal to v2 -1");
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0)==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v2;
            edge[nedges].n2 = v3;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(0);
            // assert(v1 != v2 && "v1 should not be equal to v2 -2" );
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1)==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v3;
            edge[nedges].n2 = v1;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(1);
            // assert(v1 != v2 && "v1 should not be equal to v2 -3");
         }
      }
   //  Quadrilateral element
      else if (elm[i].nvtx==4) {

         v4 = (*elm[i].vtx)(3);

         if ( (*elm[i].nghbr)(2) > i or (*elm[i].nghbr)(2) ==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v1;
            edge[nedges].n2 = v2;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(2);
            // assert(v1 != v2 && "v1 should not be equal to v2 -q1");
         }

         if ( (*elm[i].nghbr)(3) > i or (*elm[i].nghbr)(3) ==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v2;
            edge[nedges].n2 = v3;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(3);
            // assert(v1 != v2 && "v1 should not be equal to v2 -q2");
         }

         if ( (*elm[i].nghbr)(0) > i or (*elm[i].nghbr)(0) ==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v3;
            edge[nedges].n2 = v4;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(0);
            // assert(v1 != v2 && "v1 should not be equal to v2 -q3");
         }

         if ( (*elm[i].nghbr)(1) > i or (*elm[i].nghbr)(1) ==0 ) {
            nedges = nedges + 1;
            edge[nedges].n1 = v4;
            edge[nedges].n2 = v1;
            edge[nedges].e1 = i;
            edge[nedges].e2 = (*elm[i].nghbr)(1);
            // assert(v1 != v2 && "v1 should not be equal to v2 -q4");
         }

      }//    endif tri_quad2

   }//   end do elements3

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
//   edges : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      n1 = edge[i].n1;
      n2 = edge[i].n2;
      e1 = edge[i].e1;
      e2 = edge[i].e2;
      // edge centroids:
      xm = half*( node[n1].x + node[n2].x );
      ym = half*( node[n1].y + node[n2].y );

      
      edge[i].dav = zero;

/*
Q, should below be

(*edge[i].dav)(0) 


A: no, those were statically allocated, not pointers to
*/
      // Contribution from the left element
      if (e1 > 0) {
         xc = elm[e1].x;
         yc = elm[e1].y;
         edge[i].dav(0) = -(ym-yc);
         edge[i].dav(1) =   xm-xc;
      }

      // Contribution from the right element
      if (e2 > 0) {
         xc = elm[e2].x;
         yc = elm[e2].y;
         edge[i].dav(0) = edge[i].dav(0) -(yc-ym);
         edge[i].dav(1) = edge[i].dav(1) + xc-xm;
      }

      if (e1 < 0 and e2 < 0) {
         cout << "ERROR: e1 and e2 are both negative... " << endl;
      }

      // Magnitude and unit vector
      edge[i].da  = sqrt( edge[i].dav(0) * edge[i].dav(0) + \
                                 edge[i].dav(1) * edge[i].dav(1) );
      edge[i].dav = edge[i].dav / edge[i].da;

      // Edge vector

      edge[i].ev(0) = node[n2].x - node[n1].x;
      edge[i].ev(1) = node[n2].y - node[n1].y;
      edge[i].e     = sqrt( edge[i].ev(0) * edge[i].ev(0) + \
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

   for (size_t i = 0; i < nnodes; i++) {
      node[i].nnghbrs = 0;
   }

// Loop over edges and distribute the node numbers:

   //edges4 : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {

      n1 = edge[i].n1;
      n2 = edge[i].n2;
      

      // (1) Add n1 to the neighbor list of n2
      node[n1].nnghbrs = node[n1].nnghbrs + 1;
      node[n1].nghbr = new Array2D<int>(node[n1].nnghbrs, 1);
      (*node[n1].nghbr)( node[n1].nnghbrs - 1 ) = n2; // more fence post trickery


      // (2) Add n2 to the neighbor list of n1
      node[n2].nnghbrs = node[n2].nnghbrs + 1;
      node[n2].nghbr = new Array2D<int>(node[n2].nnghbrs, 1);
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
         x1 = node[ (*bound[i].bnode)[j  ][0] ].x;
         y1 = node[ (*bound[i].bnode)[j  ][0] ].y;

         x2 = node[ (*bound[i].bnode)[j+1][0] ].x;
         y2 = node[ (*bound[i].bnode)[j+1][0] ].y;

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

         (*bound[i].bn)(j)  = sqrt( (*bound[i].bnx)(j) * (*bound[i].bnx)(j) + \
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
            v2 = (*bound[i].bnode)[j][0];
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

     dsL = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
     dsR = sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) );
      dx = dsR*x1/(dsL*(-dsL-dsR))-x2/dsR+x2/dsL+dsL*x3/((dsR+dsL)*dsR);
      dy = dsR*y1/(dsL*(-dsL-dsR))-y2/dsR+y2/dsL+dsL*y3/((dsR+dsL)*dsR);

     ds  = sqrt( dx*dx + dy*dy );
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
//        /  o        Note: k-th neighbor is given by "node(j).nghbr(k)"
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
//    n2 = node[n1].nghbr(edge[i].kth_nghbr_of_1)
//    n1 = node[n2].nghbr(edge[i].kth_nghbr_of_2)
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
      node[i].bmark   = 0;
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

         (*bound[i].bfn)(j,0)  =  sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
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

         (*bound[i].bfn)(j)  =  sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
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
      //for (size_t j = 0; j < 1; j++) {

         //   bface is defined by the nodes v1 and v2.
         v1 = (*bound[i].bnode)(j) ;
         v2 = (*bound[i].bnode)(j+1);

         found = false;

      //   Find the element having the bface from the elements
      //   around the node v1.

         //do k = 1, node[v1).nelms
         for (size_t k = 0; k < node[v1].nelms; k ++) {
            //k = node[v1].nelms-1;

            ielm = (*node[v1].elm)(k, 0); //[k][0];


            //cout << "v1, k, ielm  " << v1 << "    " << k << "    " << ielm << endl;
            //do ii = 1, elm[ielm].nvtx;
            for (size_t ii = 0; ii < elm[ielm].nvtx; ii++) {


               in = ii;
               im = ii+1;
               if (im > elm[ielm].nvtx-1 ) { im = 0;}//im - (elm[ielm].nvtx-0); }//return to 0? (cannot use im = 0; }//)
              
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
                  //if (j < 2) cout << "found! " << endl;;//" " << in << " " << im << endl;
                  //if (j < 2) cout << "    " << ielm << "    " << im << "    " << in << endl;
                  // if (j < 2) cout << " v = " << vt1 << "  " << v1 << "  " << vt2 << "  " << v2 << endl;
                  //cout << "break 1" << endl;
                  break; //continue; //exit
               }
               // cout << "break 2" << endl;
               //if (found) {break;} //exit  //extra break needed to account for exit behavior!
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
            //std::exit(0);//stop
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

   ave_nghbr = node[1].nnghbrs;
   min_nghbr = node[1].nnghbrs;
   max_nghbr = node[1].nnghbrs;
      imin = 1;
      imax = 1;
   if (node[1].nnghbrs==2) {
      cout << "--- 2 neighbors for the node = " << 1 << endl;
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
      elm[i].edge = new Array2D<int>[ elm[i].nnghbrs ];
   }

   //edges3 : do i = 1, nedges
   for (size_t i = 0; i < nedges; i++) {


      e1 = edge[i].e1;
      e2 = edge[i].e2;

      // Left element
      if (e1 > 0) {
         //do k = 1, elm[e1].nnghbrs;
         for (size_t k = 0; k < elm[e1].nnghbrs; k++) {
            if ( (*elm[e1].nghbr)(k)==e2) (*elm[e1].edge)(k) = i;
         }
      }

      // Right element
      if (e2 > 0) {
         //do k = 1, elm[e2].nnghbrs;
         for (size_t k = 0; k < elm[e2].nnghbrs; k++) {
            if ( (*elm[e2].nghbr)(k)==e1)  (*elm[e2].edge)(k) = i;
         }
      }

   }//end do edges3

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
//     if (jelm > i) {
//      nfaces = nfaces + 1
//     }
//    end do
//   end do elements4

//   allocate(face(nfaces))

//   nfaces = 0

//   elements5 : do i = 1, nelms
//    do k = 1, elm(i).nnghbrs
//    jelm = elm(i).nghbr(k)

//     if (jelm > i) {

//      nfaces = nfaces + 1

//      face(nfaces).e1 = i
//      face(nfaces).e2 = jelm

//      iedge = elm(i).edge(k)
//      v1 = edge(iedge).n1
//      v2 = edge(iedge).n2

//      if (edge(iedge).e1 == jelm) {
//       face(nfaces).n1 = v1
//       face(nfaces).n2 = v2
//      else
//       face(nfaces).n1 = v2
//       face(nfaces).n2 = v1
//      }
   
//     elseif (jelm == 0) {
// //    Skip boundary faces.
//     }

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
//   face(i).dav(1) = -( node[n2].y - node[n1].y )
//   face(i).dav(2) =    node[n2].x - node[n1].x
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
//     if ( elm(i).nghbr(k) > 0 ) {
//      elm(i).nvnghbrs = elm(i).nvnghbrs + 1
//      call my_alloc_int_ptr(elm(i).vnghbr, elm(i).nvnghbrs)
//      elm(i).vnghbr(elm(i).nvnghbrs) = elm(i).nghbr(k)
//     }
//    end do

// // (2)Add vertex-neighbors
//    do k = 1, elm(i).nvtx
//     v1 = elm(i).vtx(k)

//     velms : do j = 1, node[v1).nelms
//      e1 = node[v1).elm(j)
//      if (e1 == i) cycle velms

// //    Check if the element is already added.
//        found = .false.
//      do ii = 1, elm(i).nvnghbrs
//       if ( e1 == elm(i).vnghbr(ii) ) {
//        found = .true.
//        exit
//       }
//      end do

// //    Add the element, e1, if not added yet.
//      if (.not.found) {
//       elm(i).nvnghbrs = elm(i).nvnghbrs + 1
//       call my_alloc_int_ptr(elm(i).vnghbr, elm(i).nvnghbrs)
//       elm(i).vnghbr(elm(i).nvnghbrs) = e1
//      }
//     end do velms

//    end do

//    ave_nghbr = ave_nghbr + elm(i).nvnghbrs
//    if (elm(i).nvnghbrs < min_nghbr) imin = i
//    if (elm(i).nvnghbrs > max_nghbr) imax = i
//    min_nghbr = min(min_nghbr, elm(i).nvnghbrs)
//    max_nghbr = max(max_nghbr, elm(i).nvnghbrs)
//    if (elm(i).nvnghbrs < 3) {
//     write(*,*) "--- Not enough neighbors: elm = ", i, &
//                "elm(i).nvnghbrs=",elm(i).nvnghbrs
//    }

//   end do elements7

//   write(*,*) "      ave_nghbr = ", ave_nghbr/nelms
//   write(*,*) "      min_nghbr = ", min_nghbr, " elm = ", imin
//   write(*,*) "      max_nghbr = ", max_nghbr, " elm = ", imax
//   write(*,*)


//   do i = 1, nelms
//    elm(i).bmark = 0
//   end do

//   bc_loop : do i = 1, nbound
//    if (trim(bound[i].bc_type) == "dirichlet") {
//     do j = 1, bound[i].nbfaces
//      elm( bound[i].belm(j) ).bmark = 1
//     end do
//    }
//   end do bc_loop

// //--------------------------------------------------------------------------------

   return;
};
//  end subroutine construct_grid_data

//********************************************************************************

