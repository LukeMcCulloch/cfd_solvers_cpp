


//======================================
// my simple array class template (type)
#include "tests_array.hpp"
#include "array_template.hpp"
#include "arrayops.hpp"



// //======================================
// // 2D Euler sovler
// 2D Eiuler approximate Riemann sovler
#include "EulerUnsteady2D.h"


//======================================
#include "EulerUnsteady2D_basic_package.h"

//======================================
// string trimfunctions
#include "StringOps.h"

// EulerSolver2D::MainData2D::MainData2D() {


//    M_inf             = 0.0;                  // Freestream Mach number to be set in the function
//                                              //    -> "initial_solution_shock_diffraction"
//                                              //    (Specify M_inf here for other problems.)
//    gamma             = 1.4;                  // Ratio of specific heats
//    CFL               = 0.95;                 // CFL number
//    t_final           = 0.18;                 // Final time to stop the calculation.
//    time_step_max     = 5000;                 // Max time steps (just a big enough number)
//    inviscid_flux     = "rhll";               // = Rotated-RHLL      , "roe"  = Roe flux
//    limiter_type      = "vanalbada";          // = Van Albada limiter, "none" = No limiter
//    nq        = 4;                    // The number of equtaions/variables in the target equtaion.
//    gradient_type     = "linear";             // or "quadratic2 for a quadratic LSQ.
//    gradient_weight   = "none";               // or "inverse_distance"
//    gradient_weight_p =  EulerSolver2D::one;  // or any other real value
// }




// EulerSolver2D::MainData2D::~MainData2D() {
//     delete[] node;
//     }

EulerSolver2D::Solver::Solver(){}


EulerSolver2D::Solver::~Solver(){
   printf("destruct Solver");
}


//********************************************************************************
//********************************************************************************
//********************************************************************************
//* Euler solver: Node-Centered Finite-Volume Method (Edge-Based)
//*
//* - Node-centered finite-volume method for unstructured grids(quad/tri/mixed)
//* - Roe flux with an entropy fix and Rotated-RHLL flux
//* - Reconstruction by unweighted least-squares method (2x2 system for gradients)
//* - Van Albada slope limiter to the primitive variable gradients
//* - 2-Stage Runge-Kutta time-stepping
//*
//********************************************************************************
void EulerSolver2D::Solver::euler_solver_main(EulerSolver2D::MainData2D& E2Ddata ){

   // //Local variables
   Array2D<real> res_norm(4,3); //Residual norms(L1,L2,Linf)
   Array2D<real>*  u0;       //Saved solution
   real dt, time;    //Time step and actual time
   int i_time_step; //Number of time steps
   int i;

   // // Allocate the temporary solution array needed for the Runge-Kutta method.
   // allocate(u0(nnodes,4));

   // These parameters are set in main. Here just print them on display.
   cout << " \n";
   cout << "Calling the Euler solver...";
   cout << " \n";
   cout << "                  M_inf = " <<  E2Ddata.M_inf << " \n";
   cout << "                    CFL = " <<  E2Ddata.CFL << " \n";
   cout << "             final time = " <<  E2Ddata.t_final << " \n";
   cout << "          time_step_max = " <<  E2Ddata.time_step_max << " \n";
   cout << "          inviscid_flux = " <<  trim(E2Ddata.inviscid_flux) << " \n";
   cout << "           limiter_type = " <<  trim(E2Ddata.limiter_type) << " \n";
   cout << " \n";

   //--------------------------------------------------------------------------------
   // First, make sure that normal mass flux is zero at all solid boundary nodes.
   // NOTE: Necessary because initial solution may generate the normal component.
   //--------------------------------------------------------------------------------

}
//********************************************************************************
// End of program
//********************************************************************************




void EulerSolver2D::Solver::compute_lsq_coeff_nc(EulerSolver2D::MainData2D& E2Ddata) {

      int i, in, ell, ii, k;

      cout << " \n";
      cout << " " "Constructing LSQ coefficients...\n";

      // 1. Coefficients for the linear LSQ gradients

      cout << "---(1) Constructing Linear LSQ coefficients...\n";

      // nnodes
      for (size_t i = 0; i < E2Ddata.nnodes; i++) {

         //my_alloc_p2_ptr(node[i].lsq2x2_cx,node[i].nnghbrs)
         E2Ddata.node[i].lsq2x2_cx = new Array2D<real>( E2Ddata.node[i].nnghbrs, 1 );
         //my_alloc_p2_ptr(node[i].lsq2x2_cy,node[i].nnghbrs)
         E2Ddata.node[i].lsq2x2_cy = new Array2D<real>( E2Ddata.node[i].nnghbrs, 1 );
         lsq01_2x2_coeff_nc(E2Ddata, i);

      }


// 2. Coefficients for the quadratic LSQ gradients (two-step method)

   cout << "---(2) Constructing Quadratic LSQ coefficients...\n";

   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {

      ii = 0;
      //loop node[i].nnghbrs;
      for (size_t k = 0; k < E2Ddata.node[i].nnghbrs; k++) {
         in = (*E2Ddata.node[i].nghbr)(k);
         //loop E2Ddata.node[in].nnghbrs;
         for (size_t ell = 0; ell < E2Ddata.node[in].nnghbrs; ell++) {
            ii = ii + 1;
         }//nghbr_nghbr
      }//nghbr

   //   call my_alloc_p2_ptr(node(i)%lsq5x5_cx, ii)
      E2Ddata.node[i].lsq5x5_cx = new Array2D<real>(ii,1);
   //   call my_alloc_p2_ptr(node(i)%lsq5x5_cy, ii)
      E2Ddata.node[i].lsq5x5_cy = new Array2D<real>(ii,1);
   //   call my_alloc_p2_ptr(node(i)%dx,node(i)%nnghbrs)
      //E2Ddata.node[i].dx = new Array2D<real>(E2Ddata.node[i].nnghbrs+1,1);
      E2Ddata.node[i].dx = new Array2D<real>(E2Ddata.node[i].nnghbrs,1);
   //   call my_alloc_p2_ptr(node(i)%dy,node(i)%nnghbrs)
      //E2Ddata.node[i].dy = new Array2D<real>(E2Ddata.node[i].nnghbrs+1,1);
      E2Ddata.node[i].dy = new Array2D<real>(E2Ddata.node[i].nnghbrs,1);
   //   call my_alloc_p2_matrix_ptr(node(i)%dw, nq,node(i)%nnghbrs)
      //E2Ddata.node[i].dw = new Array2D<real>(E2Ddata.nq, E2Ddata.node[i].nnghbrs+1);
      E2Ddata.node[i].dw = new Array2D<real>(E2Ddata.nq, E2Ddata.node[i].nnghbrs);

   }//end do

   lsq02_5x5_coeff2_nc(E2Ddata);



} //  end compute_lsq_coeff_nc

// //********************************************************************************
// //* This subroutine verifies the implementation of LSQ gradients.
// //*
// //* 1. Check if the linear LSQ gradients are exact for linear functions.
// //* 2. Check if the quadratic LSQ gradients are exact for quadratic functions.
// //*
// //* Note: Here, we use only the first component of u=(u1,u2,u3), i.e., ivar=1.
// //*
// //********************************************************************************
void EulerSolver2D::Solver::check_lsq_coeff_nc(EulerSolver2D::MainData2D& E2Ddata) {
   //  use edu2d_constants   , only : p2, one, two
   //  use edu2d_my_main_data, only : nnodes, node

   int       i, ix, iy, ivar;
   std::string grad_type_temp;
   real error_max_wx, error_max_wy, x, y;
   real x_max_wx, y_max_wx, x_max_wy, y_max_wy, wx, wxe, wy, wye;
   real a0, a1, a2, a3, a4, a5;

   ix = 0;
   iy = 1;

// We only use w(0) for this test.
  ivar = 1;

//---------------------------------------------------------------------
// 1. Check linear LSQ gradients
//---------------------------------------------------------------------
  cout << " \n";
  cout << "---------------------------------------------------------\n";
  cout << "---------------------------------------------------------\n";
  cout << "- Checking Linear LSQ gradients..." << endl;

//  (1). Store a linear function in w(ivar) = x + 2*y.
//       So the exact gradient is grad(w(ivar)) = (0,1).

   cout << "- Storing a linear function values... \n";
   // nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      x = E2Ddata.node[i].x;
      y = E2Ddata.node[i].y;
      (*E2Ddata.node[i].w)(ivar) = one*x + two*y;
   }

//  (2). Compute the gradient by linear LSQ

   cout << "- Computing linear LSQ gradients..\n";
   grad_type_temp = 'linear';
   compute_gradient_nc(E2Ddata, ivar, grad_type_temp);

//  (3). Compute the relative errors (L_infinity)

   cout << "- Computing the relative errors (L_infinity)..\n";
   error_max_wx = -one;
   error_max_wy = -one;
   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      error_max_wx = max( std::abs( (*E2Ddata.node[i].gradw)(ivar,ix) - one )/one, error_max_wx );
      error_max_wy = max( std::abs( (*E2Ddata.node[i].gradw)(ivar,iy) - two )/two, error_max_wy );
   }

   cout << " Max relative error in wx =  " << error_max_wx;
   cout << " Max relative error in wy =  " << error_max_wy;
   cout << "---------------------------------------------------------\n";
   cout << "---------------------------------------------------------\n";


//---------------------------------------------------------------------
// 2. Check quadratic LSQ gradients
//---------------------------------------------------------------------
   cout << "\n";
   cout << "---------------------------------------------------------\n";
   cout << "---------------------------------------------------------\n";
   cout << "- Checking Quadratic LSQ gradients...\n";

//  (1). Store a quadratic function in w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
//       So the exact gradient is grad(w(ivar)) = (a1+2*a3*x+a4*y, a2+2*a5*y+a4*x)

   a0 =    21.122;
   a1 =     1.000;
   a2 = -   1.970;
   a3 =   280.400;
   a4 = -2129.710;
   a5 =   170.999;

   cout << "- Storing a quadratic function values...\n";
   // loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      x = E2Ddata.node[i].x;
      y = E2Ddata.node[i].y;
      (*E2Ddata.node[i].w)(ivar) = a0 + a1*x + a2*y + a3*x*x + a4*x*y + a5*y*y;
   }

//  (2). Compute the gradient by linear LSQ

   cout << "- Computing quadratic LSQ gradients..\n";
   grad_type_temp = 'quadratic2';
   compute_gradient_nc(E2Ddata,ivar,grad_type_temp);

//  (3). Compute the relative errors (L_infinity)

   cout << "- Computing the relative errors (L_infinity)..\n";
   error_max_wx = -one;
   error_max_wy = -one;
   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      x = E2Ddata.node[i].x;
      y = E2Ddata.node[i].y;

      if ( std::abs( (*E2Ddata.node[i].gradw)(ivar,ix) - (a1+2.0*a3*x+a4*y) )/(a1+2.0*a3*x+a4*y) >  error_max_wx )  {
         wx  = (*E2Ddata.node[i].gradw)(ivar,ix);
         wxe = a1+2.0*a3*x+a4*y;
         error_max_wx = std::abs( wx - wxe )/wxe;
         x_max_wx = x;
         y_max_wx = y;
      }

      if ( std::abs( (*E2Ddata.node[i].gradw)(ivar,iy) - (a2+2.0*a5*y+a4*x) )/(a2+2.0*a5*y+a4*x) >  error_max_wy )  {
         wy  = (*E2Ddata.node[i].gradw)(ivar,iy);
         wye = a2+2.0*a5*y+a4*x;
         error_max_wy = std::abs( wy - wye )/wye;
         x_max_wy = x;
         y_max_wy = y;
      }

   }//end do

  cout << " Max relative error in wx = " <<  error_max_wx <<  " at (x,y) = " <<  x_max_wx <<  y_max_wx << "\n";
  cout << "   At this location, LSQ ux = " <<  wx <<  ": Exact ux = " <<  wxe << "\n";
  cout << " Max relative error in wy = " <<  error_max_wy <<  " at (x,y) = " <<  x_max_wy <<  y_max_wy << "\n";
  cout << "   At this location, LSQ uy = " <<  wy <<  ": Exact uy = " <<  wye << "\n";
  cout << "---------------------------------------------------------\n";
  cout << "---------------------------------------------------------\n";
  cout << " " << endl;


} //  end check_lsq_coeff_nc
//********************************************************************************
//*
//********************************************************************************




//********************************************************************************
//* This subroutine computes gradients at nodes for the variable u(ivar),
//* where ivar = 1,2,3, ..., or nq.
//*
//* ------------------------------------------------------------------------------
//*  Input: node[:).u(ivar)
//*
//* Output: node[i].gradu(ivar,1:2) = ( du(ivar)/dx, du(ivar)/dy )
//* ------------------------------------------------------------------------------
//********************************************************************************
void EulerSolver2D::Solver::compute_gradient_nc(EulerSolver2D::MainData2D& E2Ddata,
                                                   int ivar, std::string grad_type) {

   //use edu2d_my_main_data, only : node, nnodes

   //integer, intent(in) :: ivar
   //std::string grad_type

   int i, k, in;

   cout << "computing gradient nc type " << grad_type << endl;

  if (trim(grad_type) == "none") return;

   //-------------------------------------------------
   //  Two-step quadratic LSQ 5x5 system
   //  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at
   //        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.

   if (trim(grad_type) == "quadratic2") {

//  Perform Step 1 as below (before actually compute the gradient).

      //do i = 1, nnodes
      for (size_t i = 0; i < E2Ddata.nnodes; i++) {


         //nghbr0 : do k = 1, node[i].nnghbrs;
         for (size_t k = 0; k < E2Ddata.node[i].nnghbrs; k++) {
                       in      = (*E2Ddata.node[i].nghbr)(k);
            (*E2Ddata.node[i].dx)(k)      = E2Ddata.node[in].x       - E2Ddata.node[i].x;
            (*E2Ddata.node[i].dy)(k)      = E2Ddata.node[in].y       - E2Ddata.node[i].y;
            (*E2Ddata.node[i].dw)(ivar,k) = (*E2Ddata.node[in].w)(ivar) - (*E2Ddata.node[i].w)(ivar);
         } //end loop nghbr0 nnghbrs
      }//end loop nnodes

   }
//-------------------------------------------------

//------------------------------------------------------------
//------------------------------------------------------------
//-- Compute LSQ Gradients at all nodes.
//------------------------------------------------------------
//------------------------------------------------------------

   //nodes : loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {

      //-------------------------------------------------
      // Linear LSQ 2x2 system
      if (trim(grad_type) == "linear") {

         lsq_gradients_nc(E2Ddata, i, ivar);
      }
      //-------------------------------------------------
      // Two-step quadratic LSQ 5x5 system
      //  Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at
      //        http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
      else if (trim(grad_type) == "quadratic2") {

         lsq_gradients2_nc(E2Ddata, i, ivar);
      }
      //-------------------------------------------------
      else {

         cout << " Invalid input value -> " << trim(grad_type);
         std::exit(0); //stop

      }
      //-------------------------------------------------

   } //end loop nodes

} // end compute_gradient_nc





//********************************************************************************
//* Compute the gradient, (wx,wy), for the variable u by Quadratic LSQ.
//*
//* ------------------------------------------------------------------------------
//*  Input:            inode = Node number at which the gradient is computed.
//*                     ivar =   Variable for which the gradient is computed.
//*          node[:).w(ivar) = Solution at nearby nodes.
//*
//* Output:  node[inode].gradu = gradient of the requested variable
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
void EulerSolver2D::Solver::lsq_gradients2_nc(
   EulerSolver2D::MainData2D& E2Ddata, int inode, int ivar) {

   //Local variables
   int  in;
   int  ix, iy, ii, ell, k;
   real da, ax, ay;

   ix = 0;
   iy = 1;
   ax = zero;
   ay = zero;

//   Loop over neighbors

   ii = 0;

   //nghbr : loop node[inode].nnghbrs
   for (size_t k = 0; k < E2Ddata.node[inode].nnghbrs; k++) {
            //in = node[inode].nghbr(k)
            in = (*E2Ddata.node[inode].nghbr)(k);

      //nghbr_nghbr : do ell = 1, node[in].nnghbrs
      for (size_t ell = 0; ell < E2Ddata.node[in].nnghbrs; ell++ ) {

         da = (*E2Ddata.node[in].w)(ivar) - (*E2Ddata.node[inode].w)(ivar) +\
                                          (*E2Ddata.node[in].dw)(ivar,ell);

         if ( (*E2Ddata.node[in].nghbr)(ell) == inode ) {
            da = (*E2Ddata.node[in].w)(ivar) - (*E2Ddata.node[inode].w)(ivar);
         }

         ii = ii + 1;

         ax = ax + (*E2Ddata.node[inode].lsq5x5_cx)(ii) * da;
         ay = ay + (*E2Ddata.node[inode].lsq5x5_cy)(ii) * da;

      } // end loop nghbr_nghbr

   } // end loop nghbr

   (*E2Ddata.node[inode].gradw)(ivar,ix) = ax;  //<-- dw(ivar)/dx;
   (*E2Ddata.node[inode].gradw)(ivar,iy) = ay;  //<-- dw(ivar)/dy;

} //end lsq_gradients2_nc
//--------------------------------------------------------------------------------






//********************************************************************************
//* --- LSQ Coefficients for 2x2 Linear Least-Squares Gradient Reconstruction ---
//*
//* ------------------------------------------------------------------------------
//*  Input:  inode = node number at which the gradient is computed.
//*
//* Output:  node[inode].lsq2x2_cx(:)
//*          node[inode].lsq2x2_cx(:)
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
void EulerSolver2D::Solver::lsq01_2x2_coeff_nc(EulerSolver2D::MainData2D& E2Ddata , int inode) {


//  use edu2d_my_main_data           , only : node
//  use edu2d_constants              , only : p2, zero


//  int, intent(in) :: inode

//Local variables
//  real a(2,2), dx, dy, det, w2, w2dvar
//  int  :: k, inghbr, ix=1,iy=2
//  real(p2), dimension(2,2) :: local_lsq_inverse
   Array2D<real> a(2,2), local_lsq_inverse(2,2);
   real dx, dy, det, w2, w2dvar;
   int k, inghbr, ix,iy;
   ix = 0;
   iy = 1;

   a = zero;

//  Loop over the neighbor nodes:  node[inode].nnghbrs
   for (size_t k = 0; k < E2Ddata.node[inode].nnghbrs; k++) {
      inghbr = (*E2Ddata.node[inode].nghbr)(k);

      if (inghbr == inode) {
         cout << "ERROR: nodes must differ//" << endl;
         std::exit(0);
      }
      dx = E2Ddata.node[inghbr].x - E2Ddata.node[inode].x;
      dy = E2Ddata.node[inghbr].y - E2Ddata.node[inode].y;

      w2 = lsq_weight(E2Ddata, dx, dy);
      w2 = w2 * w2;

      a(0,0) = a(0,0) + w2 * dx*dx;
      a(0,1) = a(0,1) + w2 * dx*dy;

      a(1,0) = a(1,0) + w2 * dx*dy;
      a(1,1) = a(1,1) + w2 * dy*dy;


   }//end do

   det = a(0,0)*a(1,1) - a(0,1)*a(1,0);
   if (std::abs(det) < 1.0e-14) {
      cout << " Singular: LSQ det = " <<  det <<  " i= " <<  inode;
      std::exit(0);
    }

   // invert andE2Ddatainverse(0,1) = -a(1,0)/det;
   local_lsq_inverse(1,0) = -a(0,1)/det;
   local_lsq_inverse(1,1) =  a(0,0)/det;

   //  Now compute the coefficients for neighbors.

   //nghbr : loop node[inode].nnghbrs
   for (size_t k = 0; k < E2Ddata.node[inode].nnghbrs; k++) {
      inghbr = (*E2Ddata.node[inode].nghbr)(k);

      dx = E2Ddata.node[inghbr].x - E2Ddata.node[inode].x;
      dy = E2Ddata.node[inghbr].y - E2Ddata.node[inode].y;

      w2dvar = lsq_weight(E2Ddata, dx, dy);
      w2dvar = w2dvar * w2dvar;

      (*E2Ddata.node[inode].lsq2x2_cx)(k)  = local_lsq_inverse(ix,0)*w2dvar*dx \
                                 + local_lsq_inverse(ix,1)*w2dvar*dy;

      (*E2Ddata.node[inode].lsq2x2_cy)(k)  = local_lsq_inverse(iy,0)*w2dvar*dx \
                                 + local_lsq_inverse(iy,1)*w2dvar*dy;

   } //end nghbr loop

}//lsq01_2x2_coeff_nc
//********************************************************************************
//*
//********************************************************************************



//********************************************************************************
//* --- LSQ Coefficients for 5x5 Quadratic Least-Squares Gradient Reconstruction ---
//*
//* Note: See Nishikawa, JCP2014v273pp287-309 for details, which is available at
//*
//* http://www.hiroakinishikawa.com/My_papers/nishikawa_jcp2014v273pp287-309_preprint.pdf.
//*
//* ------------------------------------------------------------------------------
//*  Input:
//*
//* Output:  node[:).lsq5x5_cx(ii)
//*          node[:).lsq5x5_cx(ii)
//*
//* Note: This subroutine computes the LSQ coefficeints at all nodes.
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
void EulerSolver2D::Solver::lsq02_5x5_coeff2_nc(
   EulerSolver2D::MainData2D& E2Ddata ) {

//  use edu2d_my_main_data           , only : node, nnodes
//  use edu2d_constants              , only : p2, zero, half

//  implicit none

// //Local variables
//  real a(5,5), ainv(5,5), dx, dy
//  real dummy1(5), dummy2(5)
//  real w2
//  int  :: istat, ix=1, iy=2
   Array2D<real> a(5,5), ainv(5,5), dummy1(5,1), dummy2(5,1);
   real dx, dy, det, w2;
//  int i, k, ell, in, ii
   int istat;
   int ell, in, ii;
   int ix = 0;
   int iy = 1;

   cout << "     lsq02_5x5_coeff2_nc " << endl;
   // Step 1

   a = zero;

   //   node1 : loop nnodes
   //    nghbr : loop node[i].nnghbrs
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      for (size_t k = 0; k < E2Ddata.node[i].nnghbrs; k++) {
         in      = (*E2Ddata.node[i].nghbr)(k);
         // if (in >= E2Ddata.node[i].nnghbrs) {
         //    cout << "ERROR: in >= nnghbrs " << in << " " << E2Ddata.node[i].nnghbrs << endl;
         //    std::exit(0);
         // }
         //cout << " i = " << i << " k = " << k << endl;
         (*E2Ddata.node[i].dx)(k)   = E2Ddata.node[in].x - E2Ddata.node[i].x;
         (*E2Ddata.node[i].dy)(k)   = E2Ddata.node[in].y - E2Ddata.node[i].y;


      }//    end loop nghbr
   }//   end loop node1

   // Step 2

   //   node2 : loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {

     a = zero;

      //  Get dx, dy, and dw

      //nghbr2 loop k ,E2Ddata.node[i].nnghbrs
      for (size_t k = 0; k < E2Ddata.node[i].nnghbrs; k++) {

         in = (*E2Ddata.node[i].nghbr)(k);

         //nghbr_nghbr : do ell = 1, E2Ddata.node[in].nnghbrs
         for (size_t ell = 0; ell < E2Ddata.node[in].nnghbrs; ell++) {

            dx = E2Ddata.node[in].x - E2Ddata.node[i].x + (*E2Ddata.node[in].dx)(ell);
            dy = E2Ddata.node[in].y - E2Ddata.node[i].y + (*E2Ddata.node[in].dy)(ell);

            if ( (*E2Ddata.node[in].nghbr)(ell) == i ) {

               dx = E2Ddata.node[in].x - E2Ddata.node[i].x;
               dy = E2Ddata.node[in].y - E2Ddata.node[i].y;

               if ( std::abs(dx) + std::abs(dy) < 1.0e-13 ) {
                  cout << " Zero distance found at lsq02_5x5_coeff2_nc...\n";
                  cout << "    dx = " << dx << " \n";
                  cout << "    dy = " << dy << " \n";
                  cout << "- Centered node = " << i << "\n";
                  cout << "          (x,y) = " << E2Ddata.node[in].x << E2Ddata.node[in].y << "\n";
                  cout << "- Neighbor node = " << in << "\n";
                  cout << "          (x,y) = " << E2Ddata.node[in].x << E2Ddata.node[in].y << "\n";
                  std::exit(0);
               }

            }


            w2 = lsq_weight(E2Ddata, dx, dy);
            w2 = w2*w2;

            a(0,0) = a(0,0) + w2 * dx         *dx;
            a(0,1) = a(0,1) + w2 * dx         *dy;
            a(0,2) = a(0,2) + w2 * dx         *dx*dx * half;
            a(0,3) = a(0,3) + w2 * dx         *dx*dy;
            a(0,4) = a(0,4) + w2 * dx         *dy*dy * half;

            //     a(1,0) = a(1,0) + w2 * dy         *dx;
            a(1,1) = a(1,1) + w2 * dy         *dy;
            a(1,2) = a(1,2) + w2 * dy         *dx*dx * half;
            a(1,3) = a(1,3) + w2 * dy         *dx*dy;
            a(1,4) = a(1,4) + w2 * dy         *dy*dy * half;

            //     a(2,0) = a(2,0) + w2 * half*dx*dx *dx;
            //     a(2,1) = a(2,1) + w2 * half*dx*dx *dy;
            a(2,2) = a(2,2) + w2 * half*dx*dx *dx*dx * half;
            a(2,3) = a(2,3) + w2 * half*dx*dx *dx*dy;
            a(2,4) = a(2,4) + w2 * half*dx*dx *dy*dy * half;

            //     a(3,0) = a(3,0) + w2 *      dx*dy *dx;
            //     a(3,1) = a(3,1) + w2 *      dx*dy *dy;
            //     a(3,2) = a(3,2) + w2 *      dx*dy *dx*dx * half;
            a(3,3) = a(3,3) + w2 *      dx*dy *dx*dy;
            a(3,4) = a(3,4) + w2 *      dx*dy *dy*dy * half;

            //     a(4,0) = a(4,0) + w2 * half*dy*dy *dx;
            //     a(4,1) = a(4,1) + w2 * half*dy*dy *dy;
            //     a(4,2) = a(4,2) + w2 * half*dy*dy *dx*dx * half;
            //     a(4,3) = a(4,3) + w2 * half*dy*dy *dx*dy;
            a(4,4) = a(4,4) + w2 * half*dy*dy *dy*dy * half;

         } //end do nghbr_nghbr

      } //end do nghbr2

      //   Fill symmetric part

      a(1,0) = a(0,1);
      a(2,0) = a(0,2);  a(2,1) = a(1,2);
      a(3,0) = a(0,3);  a(3,1) = a(1,3); a(3,2) = a(2,3);
      a(4,0) = a(0,4);  a(4,1) = a(1,4); a(4,2) = a(2,4); a(4,3) = a(3,4);

      //   Invert the matrix

      dummy1 = zero;
      dummy2 = zero;
      //gewp_solve(a,dummy1,dummy2,ainv,istat, 5);
      //ainv = 0.;// 
      GSinv(a,dummy1,dummy2);

      // if (ainv.istat>=0) {
      //    cout << "Problem in solving the linear system//: Quadratic_LSJ_Matrix \n";
      //    std::exit(0);
      // }

      //  Now compute the coefficients for neighbors.

      ii = -1;

      //nghbr3 : do k = 1, E2Ddata.node[i].nnghbrs
      for (size_t k = 0; k < E2Ddata.node[i].nnghbrs; k++) {
         in = (*E2Ddata.node[i].nghbr)(k);

         //nghbr_nghbr2 : do ell = 1, E2Ddata.node[in].nnghbrs
         for (size_t ell = 0; ell < E2Ddata.node[in].nnghbrs; ell++) {

            dx = E2Ddata.node[in].x - E2Ddata.node[i].x + (*E2Ddata.node[in].dx)(ell);
            dy = E2Ddata.node[in].y - E2Ddata.node[i].y + (*E2Ddata.node[in].dy)(ell);

            if ( (*E2Ddata.node[in].nghbr)(ell) == i )  {
               dx = E2Ddata.node[in].x - E2Ddata.node[i].x;
               dy = E2Ddata.node[in].y - E2Ddata.node[i].y;
            }

            ii = ii + 1;

            w2 = lsq_weight(E2Ddata, dx, dy);
            w2 = w2 * w2;

 //  Multiply the inverse LSQ matrix to get the coefficients: cx(:) and cy(:):

            (*E2Ddata.node[i].lsq5x5_cx)(ii)  =   ainv(ix,0)*w2*dx  \
                                                + ainv(ix,1)*w2*dy             \
                                                + ainv(ix,2)*w2*dx*dx * half   \
                                                + ainv(ix,3)*w2*dx*dy          \
                                                + ainv(ix,4)*w2*dy*dy * half;

            (*E2Ddata.node[i].lsq5x5_cy)(ii)  =   ainv(iy,0)*w2*dx  \
                                                + ainv(iy,1)*w2*dy             \
                                                + ainv(iy,2)*w2*dx*dx * half   \
                                                + ainv(iy,3)*w2*dx*dy          \
                                                + ainv(iy,4)*w2*dy*dy * half;
         }//end do nghbr_nghbr2

      }//end do nghbr3

   }//   end do node2

} //end  lsq02_5x5_coeff2_nc
//********************************************************************************
//*
//********************************************************************************


//****************************************************************************
//* Compute the LSQ weight
//*
//* Note: The weight computed here is the square of the actual LSQ weight.
//*****************************************************************************
real EulerSolver2D::Solver::lsq_weight(EulerSolver2D::MainData2D& E2Ddata, real dx, real dy) {

   //  use edu2d_constants   , only : p2, one
   //  use edu2d_my_main_data, only : gradient_weight, gradient_weight_p

   //Output
   real lsq_weight;

   //Local
   real distance;

   if (trim(E2Ddata.gradient_weight) == "none") {

      lsq_weight = one;
   }
   else if (trim(E2Ddata.gradient_weight) == "inverse_distance") {

      distance = sqrt(dx*dx + dy*dy);

      real val = pow(distance, E2Ddata.gradient_weight_p);
      if (val < 1.e-6) {
         cout << "ERROR: distance is 0" << endl;
         std::exit(0);
      }

      lsq_weight = one / val;

   }
   return lsq_weight;
} //end lsq_weight






//********************************************************************************
//* Compute the gradient, (wx,wy), for the variable u by Linear LSQ.
//*
//* ------------------------------------------------------------------------------
//*  Input:            inode = Node number at which the gradient is computed.
//*                     ivar =   Variable for which the gradient is computed.
//*          node[:).w(ivar) = Solution at nearby nodes.
//*
//* Output:  node[inode].gradu = gradient of the requested variable
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
void EulerSolver2D::Solver::lsq_gradients_nc(EulerSolver2D::MainData2D& E2Ddata, int inode, int ivar) {

   //  use edu2d_my_main_data           , only : node
   //  use edu2d_constants              , only : p2, zero

   //  implicit none

   //  integer, intent(in) :: inode, ivar

   //Local variables
   int in, inghbr;
   int ix, iy;
   real da, ax, ay;

   ix = 0;
   iy = 1;
   ax = zero;
   ay = zero;

//   Loop over neighbors

   // loop in = 1, E2Ddata.node[inode].nnghbrs
   for (size_t in = 0; in < E2Ddata.node[inode].nnghbrs; in ++) {
      inghbr = (*E2Ddata.node[inode].nghbr)(in);

      da = (*E2Ddata.node[inghbr].w)(ivar) - (*E2Ddata.node[inode].w)(ivar);

      ax = ax + (*E2Ddata.node[inode].lsq2x2_cx)(in)*da;
      ay = ay + (*E2Ddata.node[inode].lsq2x2_cy)(in)*da;

   }

   (*E2Ddata.node[inode].gradw)(ivar,ix) = ax;  //<-- du(ivar)/dx
   (*E2Ddata.node[inode].gradw)(ivar,iy) = ay;  //<-- du(ivar)/dy

}// end lsq_gradients_nc
//--------------------------------------------------------------------------------



//****************************************************************************
//* ---------------------------- GAUSS Seidel -------------------------------
//*
//*  This computes the inverse of an (Nm)x(Nm) matrix "ai".
//*
//*****************************************************************************
Array2D<real> EulerSolver2D::Solver::GSinv(const Array2D<real>& a,
                                    Array2D<real>& m, const Array2D<real>& b) {
   return GaussSeidelInv(a,m,b);
}