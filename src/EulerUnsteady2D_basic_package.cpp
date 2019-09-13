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
//*  1. module edu2d_constants      - Some numerical values, e.g., zero, one, pi, etc.
//*  2. module edu2d_grid_data_type - Grid data types: node, edge, face, element, etc.
//*  3. module edu2d_my_main_data   - Parameters and arrays mainly used in a solver.
//*  4. module edu2d_my_allocation  - Subroutines for dynamic allocation
//*  5. module edu2d_grid_data      - Subroutines for reading/constructing/checking grid data
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
//* 1. module edu2d_constants
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

using std::cout;
using std::endl;

edu2d_my_main_data::MainData2D::MainData2D(){}
edu2d_my_main_data::MainData2D::~MainData2D(){}



//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//********************************************************************************
//* 5. module edu2d_grid_data
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
//namespace edu2d_grid_data{



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
//*    elm(1:nelms)%nvtx   =  Number of vertices of each element
//*    elm(1:nelms)%vtx(:) = Pointer to vertices of each element
//*
//* 3. Node data: nodes are stored in a 1D array
//*    node(1:nnodes)%x     = x-coordinate of the nodes
//*    node(1:nnodes)%y     = y-coordinate of the nodes
//*
//* 4. Boundary Data:
//*    nbound                   = Number of boundary segments
//*    bound(1:nbound)%nbnodes  = Number of nodes in each segment
//*    bound(1:nbound)%bnode(:) = List of node numbers for each segment
//*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
//*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
//*
//********************************************************************************
void edu2d_my_main_data::MainData2D::read_grid(std::string datafile_grid_in, 
                                               std::string datafile_bcmap_in)
{

    //use edu2d_my_main_data, only : nnodes, node, ntria, nquad, nelms, elm, nbound, bound

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
    // while (std::getline(infile, line))
    // {
    //     std::istringstream iss(line);
    //     int a, b;
    //     if (!(iss >> a >> b)) { break; } //error
        
    //     //process data here
    // }

    // READ: Get the size of the grid.
    //read(1,*) nnodes, ntria, nquad
    //nelms = ntria + nquad
    std::getline(infile, line);
    std::istringstream iss(line);
    iss >> nnodes >> ntria >> nquad;
    nelms = ntria + nquad;
    cout << "found data" << endl;
    cout << "nnodes = " << nnodes << endl;
    cout << "ntria = " << ntria << endl;
    cout << "nquad = " << nquad << endl;
    cout << "nelms = " << nelms << endl;
    

    // //  Allocate node and element arrays.
    // allocate(node(nnodes))
    // allocate(elm(  nelms))
    //node = new node_type[nnodes]; //??
    edu2d_grid_data_type::node_type* node = new edu2d_grid_data_type::node_type[nnodes]; //??
    edu2d_grid_data_type::elm_type*  elm = new edu2d_grid_data_type::elm_type[nelms];


    // // READ: Read the nodal coordinates
    // do i = 1, nnodes
    // read(1,*) node(i)%x, node(i)%y
    // end do

    // while (std::getline(infile, line))
    // {
    //     std::istringstream iss(line);
    //     int a, b;
    //     if (!(iss >> a >> b)) { break; } //error
    for (size_t i = 1; i < nnodes; i++) {
        std::getline(infile, line);
        std::istringstream iss(line);
        iss >> node[i-1].x >> node[i-1].y ;
    }
        

    // // Read element-connectivity information

    // // Triangles: assumed that the vertices are ordered counterclockwise
    // //
    // //         v3
    // //         /\
    // //        /  \
    // //       /    \
    // //      /      \
    // //     /        \
    // //    /__________\
    // //   v1           v2

    // // READ: read connectivity info for triangles
    if (ntria > 0) {
        for (size_t i = 1; i < ntria; i++) {
            std::getline(infile, line);
            std::istringstream in(line);
            elm[i].nvtx = 3;
            //allocate(elm(i)%vtx(3))
            elm[i].vtx = new Array2D<int>(3,1) ;
            //edu2d_grid_data_type::elm_type[i].vtx* = new Array2D<int>(3,1) ;
            //edu2d_grid_data_type::elm_type*  elm = new edu2d_grid_data_type::elm_type[nelms];

            std::string type;
            in >> type;                  //and read the first whitespace-separated token


            float x, y, z;
            in >> x >> y >> z;       //now read the whitespace-separated floats
            //read(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3)
            //iss >> elm[i].vtx[0] >> elm[i].vtx[1] >> elm[i].vtx[2];
            elm[i].vtx[0] = x;
            elm[i].vtx[1] = y;
            elm[i].vtx[2] = z;
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
    // if ( nquad > 0 ) then
    // do i = 1, nquad
    //     elm(ntria+i)%nvtx = 4
    //     allocate( elm(ntria+i)%vtx(4))
    //     read(1,*) elm(ntria+i)%vtx(1), elm(ntria+i)%vtx(2), &
    //             elm(ntria+i)%vtx(3), elm(ntria+i)%vtx(4)

    if (nquad > 0) {
        for (size_t i = 1; i < nquad; i++) {
            std::getline(infile, line);
            std::istringstream in(line);
            elm[i].nvtx = 4;
            //allocate(elm(i)%vtx(3))
            elm[i].vtx = new Array2D<int>(4,1) ;
            //edu2d_grid_data_type::elm_type[i].vtx* = new Array2D<int>(3,1) ;
            //edu2d_grid_data_type::elm_type*  elm = new edu2d_grid_data_type::elm_type[nelms];

            std::string type;
            in >> type;                  //and read the first whitespace-separated token


            float x1,x2,x3,x4;
            in >> x1 >> x2 >> x3 >> x4;       //now read the whitespace-separated floats
            //read(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3)
            //iss >> elm[i].vtx[0] >> elm[i].vtx[1] >> elm[i].vtx[2];
            elm[i].vtx[0] = x1;
            elm[i].vtx[1] = x2;
            elm[i].vtx[2] = x3;
            elm[i].vtx[2] = x4;
        }
    }
    // end do
    // endif



    // //  Write out the grid data.

    // write(*,*)
    // write(*,*) " Total numbers:"
    // write(*,*) "      nodes = ", nnodes
    // write(*,*) "  triangles = ", ntria
    // write(*,*) "      quads = ", nquad
    // write(*,*)
    cout << " Total numbers:" << endl;
    cout << "       nodes = " << nnodes << endl;
    cout << "   triangles = " << ntria << endl;
    cout << "       nquad = " << nquad << endl;
    cout << "       nelms = " << nelms << endl;
    

    // // Read the boundary grid data

    // // READ: Number of boundary condition types
    // read(1,*) nbound
    // allocate(bound(nbound))

    // // READ: Number of Boundary nodes (including the starting one at the end if
    // // it is closed such as an airfoil.)
    // do i = 1, nbound
    // read(1,*) bound(i)%nbnodes
    // allocate(bound(i)%bnode(bound(i)%nbnodes))
    // end do

    // // READ: Read boundary nodes
    // do i = 1, nbound
    // do j = 1, bound(i)%nbnodes
    // read(1,*) bound(i)%bnode(j)
    // end do
    // end do

    // //  Print the boundary grid data.
    // write(*,*) " Boundary nodes:"
    // write(*,*) "    segments = ", nbound
    //     do i = 1, nbound
    //     write(*,'(a9,i3,2(a11,i5))') " boundary", i, "  bnodes = ", bound(i)%nbnodes, &
    //                                                 "  bfaces = ", bound(i)%nbnodes-1
    //     end do
    // write(*,*)

    // close(1)

    // // End of Read grid file>: datafile_grid_in
    // //--------------------------------------------------------------------------------

    // //--------------------------------------------------------------------------------
    // // 2. Read the boundary condition data file

    // write(*,*)
    // write(*,*) "Reading the boundary condition file....", datafile_bcmap_in
    // write(*,*)

    // // Open the input file.
    // open(unit=2, file=datafile_bcmap_in, status="unknown", iostat=os)

    //     read(2,*) 

    // // READ: Read the boundary condition type
    // do i = 1, nbound
    //     read(2,*) dummy_int, bound(i)%bc_type
    // end do

    // //  Print the data
    //     write(*,*) " Boundary conditions:"
    // do i = 1, nbound
    //     write(*,'(a9,i3,a12,a35)') " boundary", i, "  bc_type = ", trim(bound(i)%bc_type)
    // end do

    //     i = dummy_int //Never mind. Just to avoid a compilation warning.

    //     write(*,*)

    // close(2)
    infile.close();

    // End of Read the boundary condition data file
    //--------------------------------------------------------------------------------
    return;

 } // end function read_grid
