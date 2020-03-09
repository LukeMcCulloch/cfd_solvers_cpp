
// //======================================
// // 2D Euler sovler
// 2D Eiuler approximate Riemann sovler
#include "../include/EulerUnsteady2D.h"


//======================================
#include "../include/EulerUnsteady2D_basic_package.h"

//======================================
// string trimfunctions
//#include "../include/StringOps.h" 

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
   // cout << "          inviscid_flux = " <<  trim(E2Ddata.inviscid_flux) << " \n";
   // cout << "           limiter_type = " <<  trim(E2Ddata.limiter_type) << " \n";
   cout << " \n";

   //--------------------------------------------------------------------------------
   // First, make sure that normal mass flux is zero at all solid boundary nodes.
   // NOTE: Necessary because initial solution may generate the normal component.
   //--------------------------------------------------------------------------------

}
//********************************************************************************
// End of program
//********************************************************************************




void EulerSolver2D::Solver::compute_lsq_coeff_nc() {

} //  end compute_lsq_coeff_nc

// !********************************************************************************
// !* This subroutine verifies the implementation of LSQ gradients.
// !*
// !* 1. Check if the linear LSQ gradients are exact for linear functions.
// !* 2. Check if the quadratic LSQ gradients are exact for quadratic functions.
// !*
// !* Note: Here, we use only the first component of u=(u1,u2,u3), i.e., ivar=1.
// !*
// !********************************************************************************
//void EulerSolver2D::Solver::check_lsq_coeff_nc() {
   

//  use edu2d_constants   , only : p2, one, two
//  use edu2d_my_main_data, only : nnodes, node

//  integer       :: i, ix, iy, ivar
//  character(80) :: grad_type_temp
//  real(p2)      :: error_max_wx, error_max_wy, x, y
//  real(p2)      :: x_max_wx, y_max_wx, x_max_wy, y_max_wy, wx, wxe, wy, wye
//  real(p2)      :: a0, a1, a2, a3, a4, a5

//   ix = 1
//   iy = 2

// ! We only use w(1) for this test.
//   ivar = 1
 
// !---------------------------------------------------------------------
// ! 1. Check linear LSQ gradients
// !---------------------------------------------------------------------
//   write(*,*)
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "- Checking Linear LSQ gradients..."

// !  (1). Store a linear function in w(ivar) = x + 2*y.
// !       So the exact gradient is grad(w(ivar)) = (1,2).

//    write(*,*) "- Storing a linear function values..."
//    do i = 1, nnodes
//     x = node(i)%x
//     y = node(i)%y
//     node(i)%w(ivar) = one*x + two*y
//    end do

// !  (2). Compute the gradient by linear LSQ

//    write(*,*) "- Computing linear LSQ gradients.."
//    grad_type_temp = 'linear'
//    call compute_gradient_nc(ivar,grad_type_temp)

// !  (3). Compute the relative errors (L_infinity)

//    write(*,*) "- Computing the relative errors (L_infinity).."
//    error_max_wx = -one
//    error_max_wy = -one
//    do i = 1, nnodes
//     error_max_wx = max( abs( node(i)%gradw(ivar,ix) - one )/one, error_max_wx )
//     error_max_wy = max( abs( node(i)%gradw(ivar,iy) - two )/two, error_max_wy )
//    end do

//   write(*,*) " Max relative error in wx = ", error_max_wx
//   write(*,*) " Max relative error in wy = ", error_max_wy
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "---------------------------------------------------------"


// !---------------------------------------------------------------------
// ! 2. Check quadratic LSQ gradients
// !---------------------------------------------------------------------
//   write(*,*)
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "- Checking Quadratic LSQ gradients..."

// !  (1). Store a quadratic function in w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
// !       So the exact gradient is grad(w(ivar)) = (a1+2*a3*x+a4*y, a2+2*a5*y+a4*x)

//    a0 =    21.122_p2
//    a1 =     1.000_p2
//    a2 = -   1.970_p2
//    a3 =   280.400_p2
//    a4 = -2129.710_p2
//    a5 =   170.999_p2

//    write(*,*) "- Storing a quadratic function values..."
//    do i = 1, nnodes
//     x = node(i)%x
//     y = node(i)%y
//     node(i)%w(ivar) = a0 + a1*x + a2*y + a3*x**2 + a4*x*y + a5*y**2
//    end do

// !  (2). Compute the gradient by linear LSQ

//    write(*,*) "- Computing quadratic LSQ gradients.."
//    grad_type_temp = 'quadratic2'
//    call compute_gradient_nc(ivar,grad_type_temp)

// !  (3). Compute the relative errors (L_infinity)

//    write(*,*) "- Computing the relative errors (L_infinity).."
//    error_max_wx = -one
//    error_max_wy = -one
//    do i = 1, nnodes
//     x = node(i)%x
//     y = node(i)%y

//     if ( abs( node(i)%gradw(ivar,ix) - (a1+2.0_p2*a3*x+a4*y) )/(a1+2.0_p2*a3*x+a4*y) >  error_max_wx ) then
//       wx  = node(i)%gradw(ivar,ix)
//       wxe = a1+2.0_p2*a3*x+a4*y
//       error_max_wx = abs( wx - wxe )/wxe
//       x_max_wx = x
//       y_max_wx = y
//     endif

//     if ( abs( node(i)%gradw(ivar,iy) - (a2+2.0_p2*a5*y+a4*x) )/(a2+2.0_p2*a5*y+a4*x) >  error_max_wy ) then
//       wy  = node(i)%gradw(ivar,iy)
//       wye = a2+2.0_p2*a5*y+a4*x
//       error_max_wy = abs( wy - wye )/wye
//       x_max_wy = x
//       y_max_wy = y
//     endif

//    end do

//   write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wx = ", error_max_wx, " at (x,y) = ", x_max_wx, y_max_wx
//   write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ ux = ", wx, ": Exact ux = ", wxe
//   write(*,'(a,es20.3,a,2es12.5)') " Max relative error in wy = ", error_max_wy, " at (x,y) = ", x_max_wy, y_max_wy
//   write(*,'(a,es20.10,a,es20.10)')  "   At this location, LSQ uy = ", wy, ": Exact uy = ", wye
//   write(*,*) "---------------------------------------------------------"
//   write(*,*) "---------------------------------------------------------"
//   write(*,*)


//}//  end check_lsq_coeff_nc
