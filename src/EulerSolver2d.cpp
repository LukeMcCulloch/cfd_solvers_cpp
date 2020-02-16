
// //======================================
// // 2D Euler sovler
// #include "../include/EulerUnsteady2D_basic_package.h"
// #include "../include/EulerUnsteady2D.h"



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