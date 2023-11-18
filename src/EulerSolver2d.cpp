
//======================================
// eigen for fast matrix ops
#include "GetEigen.h" //since you overload real in your code, you must include Eigen first

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

// math
#include <math.h>       /* copysign */

//using Eigen::Dynamic;
using Eigen::MatrixXd;

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
   Array2D<real> u0(E2Ddata.nnodes,E2Ddata.nq);  //Saved solution
   real dt, time;    //Time step and actual time
   //int i_time_step; //Number of time steps
   int i;

   std::string  timestep_tec      = "shock_tecplot.dat"; 

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

   eliminate_normal_mass_flux(E2Ddata);

   //--------------------------------------------------------------------------------
   // Time-stepping toward the final time
   //--------------------------------------------------------------------------------
   time = zero;


   for ( int i_time_step = 0; i_time_step < 10; ++i_time_step ) {
      
      //------------------------------------------------------
      // Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
      //  1. u^*     = u^n - (dt/vol)*Res(u^n)
      //  2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*Res(u^*)]
      //------------------------------------------------------

      
      //-----------------------------
      //- 1st Stage of Runge-Kutta:
      //-----------------------------
      compute_residual_ncfv(E2Ddata);

      
      //   Compute Res(u^n)
      residual_norm(E2Ddata, res_norm);
      
      if (i_time_step==0) {
         timestep_tec = "shock_time_s_"+ std::to_string(i_time_step-1) + ".dat";
         E2Ddata.write_tecplot_file(timestep_tec);
         std::cout << "Density    X-momentum  Y-momentum   Energy" << std::endl;
         std::cout << "t= " << time << " steps= " << i_time_step << " L1(res)= " << res_norm(0,0) << " " << res_norm(1,0) << " " << res_norm(2,0) << " " << res_norm(3,0) << std::endl;
      } else { 
         std::cout << "t= " << time << " steps= " << i_time_step << " L1(res)= " << res_norm(0,0) << " " << res_norm(1,0) << " " << res_norm(2,0) << " " << res_norm(3,0) << std::endl;
      }



      
      //   Stop if the final time is reached.
      // if ( abs(time-t_final) < 1.0e-12_p2 ) exit time_step  //tlm todo exit loop

      // }

      for (size_t i=0; i<E2Ddata.nnodes; ++i) {
         u0(i,0) = (*E2Ddata.node[i].u)(0);
         u0(i,1) = (*E2Ddata.node[i].u)(1);
         u0(i,2) = (*E2Ddata.node[i].u)(2);
         u0(i,3) = (*E2Ddata.node[i].u)(3);
      }
      
         timestep_tec = "bc_update_"+ std::to_string(i_time_step) + ".dat";
         E2Ddata.write_tecplot_file(timestep_tec);

      //   Compute the time step (local and global)
      compute_time_step(E2Ddata, dt);
      //   Adjust dt so as to finish exactly at the final time
      if (time + dt > E2Ddata.t_final) dt = E2Ddata.t_final - time;

      
      // Update the solution
      // 1st Stage => u^* = u^n - dt/dx*Res(u^n)
      update_solution(E2Ddata, one,dt,E2Ddata.CFL);

      //-----------------------------
      //- 2nd Stage of Runge-Kutta:
      //-----------------------------

      //   Compute Res(u^*)
      compute_residual_ncfv(E2Ddata);

      //   Compute 1/2*(u^n + u^*)
      for (size_t i=0; i<E2Ddata.nnodes; ++i) {
         (*E2Ddata.node[i].u)(0) = half*( (*E2Ddata.node[i].u)(0) + u0(i,0) );
         (*E2Ddata.node[i].u)(1) = half*( (*E2Ddata.node[i].u)(1) + u0(i,1) );
         (*E2Ddata.node[i].u)(2) = half*( (*E2Ddata.node[i].u)(2) + u0(i,2) );
         (*E2Ddata.node[i].u)(3) = half*( (*E2Ddata.node[i].u)(3) + u0(i,3) );
      }//end do

      //   2nd Stage => u^{n+1} = 1/2*(u^n + u^*) - 1/2*dt/dx*Res(u^*)
      update_solution(E2Ddata, half,dt,E2Ddata.CFL);

      time = time + dt;

      //if ( i_time_step % 5 .eq. 0.0) {
         //write(picCount,*) i_time_step
         timestep_tec = "shock_time_s_"+ std::to_string(i_time_step) + ".dat";
         E2Ddata.write_tecplot_file(timestep_tec);
      //}


   }//end loop time_step
   //--------------------------------------------------------------------------------
   // End of Time-stepping to the final time.
   //--------------------------------------------------------------------------------

//********************************************************************************
// End of program
//********************************************************************************
}



//********************************************************************************
//* Prepararion for Tangency condition (slip wall):
//*
//* Eliminate normal mass flux component at all solid-boundary nodes at the
//* beginning. The normal component will never be changed in the solver: the
//* residuals will be constrained to have zero normal component.
//*
//********************************************************************************
void EulerSolver2D::Solver::eliminate_normal_mass_flux( EulerSolver2D::MainData2D& E2Ddata ) {

   int i, j, inode;
   real normal_mass_flux;
   Array2D<real> n12(2,1);

   //  bc_loop : loop nbound
   for (size_t i = 0; i < E2Ddata.nbound; i++) {

      // only_slip_wall : if (trim(E2Ddata.bound[i].bc_type) == "slip_wall") then
      if ( trim(E2Ddata.bound[i].bc_type) == "slip_wall" ) {

         cout << " Eliminating the normal momentum on slip wall boundary " << i << " \n";

         // bnodes_slip_wall : loop E2Ddata.bound[i].nbnodes
         for (size_t j = 0; j < E2Ddata.bound[i].nbnodes; j++ ) {

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
            //
            // NOTE: This is a corner point between the inflow boundary and
            //       the lower-left wall. Enforce zero y-momentum, which is
            //       not ensured by the standard BCs.
            //       This special treatment is necessary because the domain
            //       is rectangular (the left boundary is a straight ine) and
            //       the midpoint node on the left boundary is actually a corner.
            //
            //       Our computational domain:
            //
            //                 ---------------
            //          Inflow |             |
            //                 |             |  o: Corner node
            //          .......o             |
            //            Wall |             |  This node is a corner.
            //                 |             |
            //                 ---------------
            //
            //       This is to simulate the actual domain shown below:
            //      
            //         -----------------------
            // Inflow  |                     |
            //         |                     |  o: Corner node
            //         --------o             |
            //            Wall |             |
            //                 |             |
            //                 ---------------
            //      In effect, we're simulating this flow by a simplified
            //      rectangular domain (easier to generate the grid).
            //      So, an appropriate slip BC at the corner node needs to be applied,
            //      which is "zero y-momentum", and that's all.
            //
            if (i==1 && j==0) {
               inode                       = (*E2Ddata.bound[i].bnode)(j);
               (*E2Ddata.node[inode].u)(1) = zero;                                  // Make sure zero y-momentum.
               (*E2Ddata.node[inode].w)    = u2w( (*E2Ddata.node[inode].u) , E2Ddata );// Update primitive variables
               
               continue; // cycle bnodes_slip_wall // That's all we neeed. Go to the next.

            }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            inode = (*E2Ddata.bound[i].bnode)(j);
            n12(0) = (*E2Ddata.bound[i].bnx)(j);
            n12(1) = (*E2Ddata.bound[i].bny)(j);

            normal_mass_flux = (*E2Ddata.node[inode].u)(0)*n12(0) + (*E2Ddata.node[inode].u)(1)*n12(0); //tlm fixed

            (*E2Ddata.node[inode].u)(1) = (*E2Ddata.node[inode].u)(1) - normal_mass_flux * n12(0);
            (*E2Ddata.node[inode].u)(2) = (*E2Ddata.node[inode].u)(2) - normal_mass_flux * n12(1);

            (*E2Ddata.node[inode].w)    = u2w( (*E2Ddata.node[inode].u) , E2Ddata );

         }//end loop bnodes_slip_wall

         std::cout << " Finished eliminating the normal momentum on slip wall boundary " << i << " \n";
         std::cout << " \n";

      }//end if only_slip_wall

   }//end loop bc_loop

 } // end eliminate_normal_mass_flux

//--------------------------------------------------------------------------------




//********************************************************************************
//* Initial solution for the shock diffraction problem:
//*
//* NOTE: So, this is NOT a general purpose subroutine.
//*       For other problems, specify M_inf in the main program, and
//*       modify this subroutine to set up an appropriate initial solution.
//
//  Shock Diffraction Problem:
//
//                             Wall
//                     --------------------
// Post-shock (inflow) |                  |
// (rho,u,v,p)_inf     |->Shock (M_shock) |            o: Corner node
//    M_inf            |                  |
//              .......o  Pre-shock       |Outflow
//                Wall |  (rho0,u0,v0,p0) |
//                     |                  |
//                     |                  |
//                     --------------------
//                           Outflow
//
//********************************************************************************
void EulerSolver2D::Solver::initial_solution_shock_diffraction(
                                             EulerSolver2D::MainData2D& E2Ddata ) {

   //Local variablesinitial_solution_shock_diffraction
   int i;
   real M_shock, u_shock, rho0, u0, v0, p0;
   real gamma = E2Ddata.gamma;

   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {

      // Pre-shock state: uniform state; no disturbance has reahced yet.

      rho0  = one;
      u0    = zero;
      v0    = zero;
      p0    = one/gamma;

   // Incoming shock speed

      M_shock = 5.09;
      u_shock = M_shock * sqrt(gamma*p0/rho0);

      // Post-shock state: These values will be used in the inflow boundary condition.
       E2Ddata.rho_inf = rho0 * (gamma + one)*M_shock*M_shock/( (gamma - one)*M_shock*M_shock + two );
         E2Ddata.p_inf =   p0 * (   two*gamma*M_shock*M_shock - (gamma - one) )/(gamma + one);
         E2Ddata.u_inf = (one - rho0/E2Ddata.rho_inf)*u_shock;
         E2Ddata.M_inf = E2Ddata.u_inf / sqrt(gamma*E2Ddata.p_inf/E2Ddata.rho_inf);
         E2Ddata.v_inf = zero;

      // Set the initial solution: set the pre-shock state inside the domain.

      (*E2Ddata.node[i].w)(0) = rho0;
      (*E2Ddata.node[i].w)(1) = u0;
      (*E2Ddata.node[i].w)(2) = v0;
      (*E2Ddata.node[i].w)(3) = p0;
      (*E2Ddata.node[i].u) = w2u( (*E2Ddata.node[i].w), E2Ddata);
      //
   }
 }  // end function initial_solution_shock_diffraction






//********************************************************************************
//* This subroutine computes the residual for a node-centered finite-volume method
//*
//* ------------------------------------------------------------------------------
//*  Input: the current solution
//*
//* Output: node(:)%res = the residual computed by the current solution.
//* ------------------------------------------------------------------------------
//*
//* Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
//*       the integral of (dF/dx + dG/dy), and at the end negative sign is added
//*       so that we have dU/dt = Res at every node.
//* 
//********************************************************************************
void EulerSolver2D::Solver::compute_residual_ncfv( EulerSolver2D::MainData2D& E2Ddata ) {

   //Local variables
   Array2D<real> num_flux   = Array2D<real>(E2Ddata.nq,1);        //Numerical flux
   Array2D<real> wL = Array2D<real>(E2Ddata.nq,1);          //Left and right face values
   Array2D<real> wR = Array2D<real>(E2Ddata.nq,1);          //Left and right face values
   Array2D<real> dwL = Array2D<real>(E2Ddata.nq,1);        //Slope at left and right nodes
   Array2D<real> dwR = Array2D<real>(E2Ddata.nq,1);        //Slope at left and right nodes
   Array2D<real> dwm  = Array2D<real>(E2Ddata.nq,1);  //Re-defined slopes to be limited
   Array2D<real> dwp  = Array2D<real>(E2Ddata.nq,1);  //Re-defined slopes to be limited
   Array2D<real> dwij   = Array2D<real>(E2Ddata.nq,1);  //Re-defined slopes to be limited
   Array2D<real> e12   = Array2D<real>(2,1);                      //Unit edge vector
   Array2D<real> n12   = Array2D<real>(2,1);                      //Unit face normal vector
   real mag_e12;                                                  //Magnitude of the edge vector
   real mag_n12;                                                  //Magnitude of the face-normal vector
   real wsn;                                                      //Scaled maximum wave speed
   real norm_momentum;                                            //Normal component of the momentum
   Array2D<real> bfluxL   = Array2D<real>(E2Ddata.nq,1);  //Boundary flux at left/right nodes
   Array2D<real> bfluxR   = Array2D<real>(E2Ddata.nq,1);  //Boundary flux at left/right nodes
   int node1, node2;                                              //Left and right nodes of each edge
   int boundary_elm;                                              //Element adjacent to boundary face
   int n1, n2;                                                    //Left and right nodes of boundary face
   int ix=0, iy=1;

   //MatrixXd me(5,5),ime(5,5);
   
   //-------------------------------------------------------------------------
   // Gradient computations for second-order accuracy

   //  Initialization of the gradient of the primitive variables.

   for (size_t i=0; i<E2Ddata.nnodes; ++i) {
      for (size_t j=0; j<E2Ddata.nq; ++j) {
         (*E2Ddata.node[i].gradw)(j) = zero;
      }
   }

   //  Perform LSQ gradient computations in the premitive variables w=[rho,u,v,p]:
   compute_gradient_nc(E2Ddata, 0, E2Ddata.gradient_type);     // Density gradient: grad(rho)
   compute_gradient_nc(E2Ddata, 1, E2Ddata.gradient_type);     // Velocity gradient: grad(u)
   compute_gradient_nc(E2Ddata, 2, E2Ddata.gradient_type);     // Velocity gradient: grad(v)
   compute_gradient_nc(E2Ddata, 3, E2Ddata.gradient_type);     // Pressure gradient: grad(p)

   //-------------------------------------------------------------------------
   // Residual computation: interior fluxes

   for (size_t i=0; i<E2Ddata.nnodes; ++i) {
      for (size_t j=0; j<E2Ddata.nq; ++j) {
         (*E2Ddata.node[i].res)(j) = zero;
         //std::cout << "E2Ddata.node[i].res(j)= " << (*E2Ddata.node[i].res)(j) << "\n";
         //std::cout << "E2Ddata.node[i].gradw(j)= " << (*E2Ddata.node[i].gradw)(j) << "\n";
      }
      E2Ddata.node[i].wsn = zero;
   }
   // Flux computation across internal edges (to be accumulated in res(:))
   //
   //   node2              1. Extrapolate the solutions to the edge-midpoint
   //       o                 from the nodes, n1 and n2.
   //        \   face      2. Compute the numerical flux
   //         \ -------c2  3. Add it to the residual for n1, and subtract it from
   //        / \              the residual for n2.
   //   face/   \ edge
   //      /     o         Directed area is the sum of the left and the right faces.
   //    c1    node1       Left/right face is defined by the edge-midpoint and
   //                      the centroid of the left/right element.
   //                      Directed area is positive in n1 -> n2
   // 
   // (c1, c2: element centroids)
   // 
   //--------------------------------------------------------------------------------
   for (int i=0; i<E2Ddata.nedges; ++i) {
      //--------------------------------------------------------------------------------

      // Left and right nodes of the i-th edge

      node1   = E2Ddata.edge[i].n1;  // Left node of the edge
      node2   = E2Ddata.edge[i].n2;  // Right node of the edge
      n12     = E2Ddata.edge[i].dav; // This is the directed area vector (unit vector) 
      
      mag_n12 = E2Ddata.edge[i].da;  // Magnitude of the directed area vector
      e12     = E2Ddata.edge[i].ev;  // This is the vector along the edge (uniti vector)
      mag_e12 = E2Ddata.edge[i].e;   // Magnitude of the edge vector (Length of the edge)

      // Solution gradient projected along the edge
      //
      //  NOTE: The gradient is multiplied by the distance.
      //        So, it is equivalent to the solution difference.

      for  (int j=0; j< E2Ddata.nq; ++j){ 
         dwL = (*E2Ddata.node[node1].gradw)(j,ix)*e12(ix) + (*E2Ddata.node[node1].gradw)(j,iy) * e12(iy) *half*mag_e12;
         dwR = (*E2Ddata.node[node2].gradw)(j,ix)*e12(ix) + (*E2Ddata.node[node2].gradw)(j,iy) * e12(iy) *half*mag_e12;
      }
      //  It is now limiter time!
      
      // if (node1 >= E2Ddata.nnodes-1) {
      //    std::cout << "ix = " << ix << "\n";
      //    std::cout << "iy = " << iy << "\n";
      //    std::cout << "node1 = " << node1 << "\n";
      //    std::cout << "node2 = " << node2 << "\n";
      //    std::cout << "w = \n";
      //    (*E2Ddata.node[node2].w).print();

      //    // std::cout << "dwL = \n";
      //    // dwL.print();
      //    // std::cout << "dwR = \n";
      //    // dwR.print();
         
      //    // std::cout << "wL = \n";
      //    // wL.print();
      //    // std::cout << "wR = \n";
      //    // wR.print();

      // }

      //  (1) No limiter (good for smooth solutions)

      //limiter : if (trim(limiter_type) == "none") then
      if (E2Ddata.limiter_type == "none") {

         //  Simple linear extrapolation
         wL = (*E2Ddata.node[node1].w) + dwL;
         wR = (*E2Ddata.node[node2].w) - dwR;


      //  (2) UMUSCL-type limiters: simple 1D limiting.
      } else if (E2Ddata.limiter_type == "vanalbada") {

         //       In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} => limited_slope = limiter(dwm,dwp)
         //
         //       We can do the same in 2D as follows.
         //
         //       In 2D:    dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
         //              => limited_slope = limiter(dwm,dwp)
         //
         //      NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
         //            so that the average (dwm+dwp)/2 will be the central-difference just like in 1D.

         //     Edge derivative
         //dwij = half*(node(node2)%w - node(node1)%w)
         dwij = half*( (*E2Ddata.node[node2].w) - (*E2Ddata.node[node1].w) );

         // //     Left face value (wL) with the Van Albada limiter
         dwm  = two*dwL-dwij;
         dwp  = dwij;
         //wL  = node(node1)%w + va_slope_limiter(dwm,dwp,mag_e12)
         wL  = (*E2Ddata.node[node1].w) + va_slope_limiter(E2Ddata, dwm,dwp,mag_e12);

         //     Right face value (wR) with the Van Albada limiter
         dwm  = -(two*dwR-dwij);
         dwp  = -dwij;
         //wR  = node(node2)%w + va_slope_limiter(dwm,dwp,mag_e12)
         wR = (*E2Ddata.node[node2].w) + va_slope_limiter(E2Ddata, dwm,dwp,mag_e12);

      //  (3) No other limiters available.
      } else {

      //std::cout << " Invalid input for limiter_type = ", E2Ddata.limiter_type << std::endl;
      std::cout << " Choose none or vanalbada, and try again."<< "\n";
      std::cout <<  " ... Stop."<< "\n";
      //stop

      } //endif limiter





      
   //  Compute the numerical flux for given wL and wR.

   //  (1) Roe flux (carbuncle is expected for strong shocks)
   if     (E2Ddata.inviscid_flux == "roe") {
      
      roe(E2Ddata,wL,wR,n12, num_flux,wsn);
      // if (node1 >= E2Ddata.nnodes-1) {
      //    std::cout << "num_flux = \n";
      //    num_flux.print();
      //    std::cout << "wsn = " << wsn << "\n";
      // }

   //  (2) Rotated-RHLL flux (no carbuncle is expected for strong shocks)
   } else if (E2Ddata.inviscid_flux =="rhll") {
      
      //TLM todo: implement the rotated rhll flux:
      rotated_rhll(E2Ddata,wL,wR,n12, num_flux,wsn);

   } else {
      
      std::cout <<  " Invalid input for inviscid_flux = " << E2Ddata.inviscid_flux << std::endl;
      std::cout <<  " Choose roe or rhll, and try again." << std::endl;
      std::cout <<  " ... Stop. \n";
      //stop

   }

   //  Add the flux multiplied by the magnitude of the directed area vector to node1,
   //  and accumulate the max wave speed quantity for use in the time step calculation.

   (*E2Ddata.node[node1].res) = (*E2Ddata.node[node1].res)  +  num_flux * mag_n12;
   E2Ddata.node[node1].wsn = E2Ddata.node[node1].wsn  +       wsn * mag_n12;

   // Subtract the flux multiplied by the magnitude of the directed area vector from node2,
   // and accumulate the max wave speed quantity for use in the time step calculation.
   //
   // NOTE: Subtract because the outward face normal is -n12 for the node2.

   (*E2Ddata.node[node2].res) = (*E2Ddata.node[node2].res)  -  num_flux * mag_n12;
   E2Ddata.node[node2].wsn = E2Ddata.node[node2].wsn  +       wsn * mag_n12;



   // std::cout << "E2Ddata.edge[i].dav(1) = " << E2Ddata.edge[i].dav(0) << "\n";
   // std::cout << "E2Ddata.edge[i].dav(2) = " << E2Ddata.edge[i].dav(1) << "\n";
   // std::cout << "n12 = " << n12(0) << "\n";
   // std::cout << "n12 = " << n12(1) << "\n";
   // std::cout << "ix = " << ix << "\n";
   // std::cout << "iy = " << iy << "\n";

   // std::cout << "node1 = " << node1 << "\n";
   // std::cout << "node2 = " << node2 << "\n";
   // std::cout << "mag_n12 = " << mag_n12 << "\n";
   // std::cout << "num_flux = \n";
   // num_flux.print();
   // std::cout << "wave speed = " << wsn << "\n";
   // std::cout << "res[node1] =  \n";
   // (*E2Ddata.node[node1].res).print();
   // std::cout << "res[node2] =  \n";
   // (*E2Ddata.node[node2].res).print();
   // std::cout << "w[node1] = " << E2Ddata.node[node1].wsn << "\n";
   // std::cout << "w[node1] = " << E2Ddata.node[node2].wsn << "\n";
   // std::cout << "--------------------------- \n\n";

   // std::cout << "dwL = \n";
   // dwL.print();
   // std::cout << "dwR = \n";
   // dwR.print();
   
   // std::cout << "wL = \n";
   // wL.print();
   // std::cout << "wR = \n";
   // wR.print();

  
   //--------------------------------------------------------------------------------
   } //end loop edges
   //--------------------------------------------------------------------------------

   //std::exit(0);




   //-------------------------------------------------------------------------
   //Close with the boundary flux using the element-based formula that is
   //exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
   //that ensure the linear exactness for 2D/3D elements).
   //
   //     |  Interior Domain          |
   //     |        .........          |
   //     |        .       .          |
   //     |        .       .          |
   //     o--o--o-----o---------o--o--o  <- Boundary segment
   //                 n1   |   n2
   //                      v
   //                    n12 (unit face normal vector)
   //
   //NOTE: We visit each boundary face, defined by the nodes n1 and n2,
   //      and compute the flux across the boundary face: left half for node1,
   //      and the right half for node2. In the above figure, the dots indicate
   //      the control volume around the node n1. Clearly, the flux across the
   //      left half of the face contributes to the node n1. Similarly for n2.
   //
   //
   //--------------------------------------------------------------------------------
   for (size_t i=0; i<E2Ddata.nbound; ++i) { //bc_loop
      //--------------------------------------------------------------------------------

      //------------------------------------------------
      //  BC: Upwind flux via freestream values
      //
      //      NOTE: If the final solution at the boundary node is far from
      //            the freestream values, then the domain is probably is not large enough.
      if (E2Ddata.bound[i].bc_type == "freestream") {
         for (size_t j=0; j<E2Ddata.bound[i].nbfaces; ++j) { //freestream condition

            n1 = (*E2Ddata.bound[i].bnode)(j  );  //Left node
            n2 = (*E2Ddata.bound[i].bnode)(j+1);  //Right node
            n12(0) = (*E2Ddata.bound[i].bfnx)(j);     //x-component of the unit face normal vector
            n12(1) = (*E2Ddata.bound[i].bfny)(j);     //y-component of the unit face normal vector
            mag_e12 = (*E2Ddata.bound[i].bfn)(j)*half; //Half length of the boundary face, j.
            
            

            //  1. Left node
            wL = (*E2Ddata.node[n1].w);
            wR(0) = E2Ddata.rho_inf;
            wR(1) = E2Ddata.u_inf;
            wR(2) = E2Ddata.v_inf;
            wR(3) = E2Ddata.p_inf;
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxL = num_flux;
            E2Ddata.node[n1].wsn = E2Ddata.node[n1].wsn + wsn*mag_e12;

            //   2. Right node
            wL = (*E2Ddata.node[n2].w);
            wR(0) = E2Ddata.rho_inf;
            wR(1) = E2Ddata.u_inf;
            wR(2) = E2Ddata.v_inf;
            wR(3) = E2Ddata.p_inf;
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxR = num_flux;
            E2Ddata.node[n2].wsn = E2Ddata.node[n2].wsn + wsn*mag_e12;

            //   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
            boundary_elm = (*E2Ddata.bound[i].belm)(j);

            if  (E2Ddata.elm[boundary_elm].nvtx == 3) { //Triangle

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + (5.0*bfluxL + bfluxR)/6.0*mag_e12;
               (*E2Ddata.node[n2].res )= (*E2Ddata.node[n2].res) + (5.0*bfluxR + bfluxL)/6.0*mag_e12;

            } else if (E2Ddata.elm[boundary_elm].nvtx == 4) { //Quad

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + bfluxL*mag_e12;
               (*E2Ddata.node[n2].res) = (*E2Ddata.node[n2].res) + bfluxR*mag_e12;

            } else {

               std::cout <<  " Error: Element is neither tria nor quad. Stop. "; 
               //stop

            }

         }  //end do bnodes_numerical_flux_via_freestream

      } //end freestream boundary type
         

      //------------------------------------------------
      //  BC: Solid body and Supersonic outflow
      //
      //      NOTE: Basically, simply compute the physical flux, which
      //            can be done by calling the Roe flux with wR = wL.
      //            It is equivalent to the interior-extrapolation condition.
      //
      //      NOTE: Tangency condition for solid body will be applied later.


      if ( (E2Ddata.bound[i].bc_type == "slip_wall")  || ((E2Ddata.bound[i].bc_type == "outflow_supersonic") ) ) {
         
         for (size_t j=0; j<E2Ddata.bound[i].nbfaces; ++j) { //slip wall(part of it anyway) and outflow supersonic
            
            //bnodes_slip_wall
            n1 = (*E2Ddata.bound[i].bnode)(j);   // Left node
            n2 = (*E2Ddata.bound[i].bnode)(j+1); // Right node
            n12(0) = (*E2Ddata.bound[i].bfnx)(j);
            n12(1) = (*E2Ddata.bound[i].bfny)(j);
            mag_e12 = (*E2Ddata.bound[i].bfn)(j)*half;

            //   1. Left node
            wL = (*E2Ddata.node[n1].w);
            wR = wL;
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxL = num_flux;
            E2Ddata.node[n1].wsn = E2Ddata.node[n1].wsn + wsn*mag_e12;

            //   2. Right node
            wL = (*E2Ddata.node[n2].w);
            wR = wL;
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxR = num_flux;
            E2Ddata.node[n2].wsn = E2Ddata.node[n2].wsn + wsn*mag_e12;

            //   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
            boundary_elm = (*E2Ddata.bound[i].belm)(j);

            if     (E2Ddata.elm[boundary_elm].nvtx == 3) { //Triangle

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + (5.0*bfluxL + bfluxR)/6.0*mag_e12;
               (*E2Ddata.node[n2].res) = (*E2Ddata.node[n2].res) + (5.0*bfluxR + bfluxL)/6.0*mag_e12;

            } else if (E2Ddata.elm[boundary_elm].nvtx == 4) { //Quad

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + bfluxL*mag_e12;
               (*E2Ddata.node[n2].res) = (*E2Ddata.node[n2].res) + bfluxR*mag_e12;

            } else {

               std::cout << " Error: Element is neither tria nor quad. Stop. \n"; //stop

            }

         } //end do bnodes_slip_wall
      


      //------------------------------------------------
      //  BC: Subsonic Outflow - Fixed Back Pressure
      //
      //      NOTE: Fix the pressure as freestream pressure
      //            on the right side of the face (outside the domain).
      //            Assumption is that the outflow boundary is far from the body.

      } else if (E2Ddata.bound[i].bc_type == "outflow_back_pressure") {

         //bnodes_outflow : do j = 1, E2Ddata.bound[i].nbfaces
         for (size_t j=0; j<E2Ddata.bound[i].nbfaces; ++j) { 

            n1 = (*E2Ddata.bound[i].bnode)(j);   // Left node
            n2 = (*E2Ddata.bound[i].bnode)(j+1); // Right node
            n12(0) = (*E2Ddata.bound[i].bfnx)(j);
            n12(2) = (*E2Ddata.bound[i].bfny)(j);
            mag_e12 = (*E2Ddata.bound[i].bfn)(j)*half;

            //   1. Left node
            wL = (*E2Ddata.node[n1].w);
            wR    = wL;
            wR(3) = E2Ddata.p_inf; //Fix the pressure
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxL = num_flux;
            E2Ddata.node[n1].wsn = E2Ddata.node[n1].wsn + wsn*mag_e12;

            //   2. Right node
            wL = (*E2Ddata.node[n2].w);
            wR = wL;
            wR(3) = E2Ddata.p_inf; //Fix the pressure
            roe(E2Ddata,wL,wR,n12, num_flux,wsn);
            bfluxR = num_flux;
            E2Ddata.node[n2].wsn = E2Ddata.node[n2].wsn + wsn*mag_e12;

            //   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
            
            //   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
            boundary_elm = (*E2Ddata.bound[i].belm)(j);

            if (E2Ddata.elm[boundary_elm].nvtx == 3) { //Triangle

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + (5.0*bfluxL + bfluxR)/6.0*mag_e12;
               (*E2Ddata.node[n2].res) = (*E2Ddata.node[n2].res) + (5.0*bfluxR + bfluxL)/6.0*mag_e12;

            } else if (E2Ddata.elm[boundary_elm].nvtx == 4) { //Quad

               (*E2Ddata.node[n1].res) = (*E2Ddata.node[n1].res) + bfluxL*mag_e12;
               (*E2Ddata.node[n2].res) = (*E2Ddata.node[n2].res) + bfluxR*mag_e12;

            } else {

               std::cout << " Error: Element is neither tria nor quad. Stop. \n"; //stop

            }

         }//end loop bnodes_outflow

      //------------------------------------------------
      }//endif bc

   //--------------------------------------------------------------------------------
   }//end do bc_loop
   //--------------------------------------------------------------------------------


   //------------------------------------------------
   //  BC: Solid body - Slip condition (Tangency condition).
   //
   //  NOTE: It is good to enforce this condition after everything else has been done.

   //bc_loop2 : do i = 1, nbound
   for (size_t i=0; i<E2Ddata.nbound; ++i) { //bc_loop

         //only_slip_wall : if (trim(E2Ddata.bound[i].bc_type) == "slip_wall") then
      if (E2Ddata.bound[i].bc_type == "slip_wall") {

         //bnodes_slip_wall2 : do j = 1, E2Ddata.bound[i].nbnodes
         for (size_t j=0; j<E2Ddata.bound[i].nbnodes; ++j) { 

            n1 = (*E2Ddata.bound[i].bnode)(j);   // Left node
            n12(0) = (*E2Ddata.bound[i].bfnx)(j);
            n12(2) = (*E2Ddata.bound[i].bfny)(j);

            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            // THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
            // Same as in the subroutine "eliminate_normal_mass_flux" above.

            if (i==1 && j==0) {
               std::cout << "shock y mmtm 0, n1 = " << n1 << "\n";
               (*E2Ddata.node[n1].res)(2) = zero; // Make sure no updates to y-momentum.
               //cycle bnodes_slip_wall2 // That's all we neeed. Go to the next. TLM TODO: c++ cycle
            }//endif
            //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            //    Subtract the normal component of the mass flow for tangency, so that
            //    normal mass flow through the boundary will not be created at nodes
            //    after the solution update.
               
            norm_momentum = (*E2Ddata.node[n1].res)(1)*n12(0) + (*E2Ddata.node[n1].res)(2)*n12(1);
            (*E2Ddata.node[n1].res)(1) = (*E2Ddata.node[n1].res)(1) - norm_momentum*n12(0);
            (*E2Ddata.node[n1].res)(2) = (*E2Ddata.node[n1].res)(2) - norm_momentum*n12(1);

         } //end loop bnodes_slip_wall2

      } //endif only_slip_wall

   }//end loop bc_loop2

   // Switch the residual sign.

   for (size_t i=0; i<E2Ddata.nnodes; ++i) {

      (*E2Ddata.node[i].res) = -(*E2Ddata.node[i].res);

   }//end do nodes3


}



//********************************************************************************
//* This subroutine computes the explicit time-step: the minimum dt over nodes.
//*
//* ------------------------------------------------------------------------------
//*  Input: node(i)%vol = Dual volume
//*         node(i)%wsn = Sum of the max wave speed multiplied by the face length
//*
//* Output:         dt  = global time step
//*         node(:)%dt  =  local time step
//* ------------------------------------------------------------------------------
//*
//* NOTE: Local time step is computed and stored at every node, but not used.
//*       For steady problems, it can be used to accelerate the convergence.
//*
//********************************************************************************
 void EulerSolver2D::Solver::compute_time_step( EulerSolver2D::MainData2D& E2Ddata,real dt) {

   //Local variables
   real dt_min;

   dt_min = 1.0e+05;

   //--------------------------------------------------------------------------------
   for (size_t i=0; i<E2Ddata.nnodes; ++i) {
   //--------------------------------------------------------------------------------

      // Local time step: dt = volume/sum(0.5*max_wave_speed*face_area).

      E2Ddata.node[i].dt = E2Ddata.node[i].vol / E2Ddata.node[i].wsn;
      // std::cout << " E2Ddata.node[i].vol = " << E2Ddata.node[i].vol << "\n";
      // std::cout << " E2Ddata.node[i].wsn = " << E2Ddata.node[i].wsn << "\n";
      // std::cout << " E2Ddata.node[i].dt = " << E2Ddata.node[i].dt << "\n";

      // Keep the minimum dt

      if (i==0) dt_min = E2Ddata.node[i].dt;
      dt_min = std::min( dt_min, E2Ddata.node[i].dt );

   //--------------------------------------------------------------------------------
   }
   //--------------------------------------------------------------------------------

   // Global time-step

      std::cout << " dt_min = " << dt_min << "\n";
   dt = dt_min;

 }
//--------------------------------------------------------------------------------


//********************************************************************************
//* This subroutine updates the solution.
//*
//* ------------------------------------------------------------------------------
//*  Input:       coeff = coefficient for RK time-stepping
//*                  dt = global time step (not used if local time stepping)
//*                 CFL = CFL number
//*          node(:)res = the residual
//*
//* Output:   node(:)u  = updated conservative variables 
//*           node(:)w  = updated primitive variables 
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
void EulerSolver2D::Solver::update_solution( EulerSolver2D::MainData2D& E2Ddata, 
                                             real coeff, real dt, real CFL){
   //

   for (size_t i=0; i<E2Ddata.nnodes; ++i) {

      //   Solution change based on the global time stepping 
      (*E2Ddata.node[i].du) = (CFL*coeff*dt/E2Ddata.node[i].vol) * (*E2Ddata.node[i].res);

      //   Solution update
      (*E2Ddata.node[i].u) = (*E2Ddata.node[i].u) + (*E2Ddata.node[i].du); // make changes                                // Make sure zero y-momentum.
      (*E2Ddata.node[i].w) = u2w( (*E2Ddata.node[i].u) , E2Ddata );// Update primitive variables
          
   }

   
   // Check the density/pressure positivity.
   //
   //  NOTE: Well, let's just print the warning message.
   //        We may stop immediately, though.

}

//********************************************************************************
//* This subroutine computes the residual norms: L1, L2, L_infty
//*
//* ------------------------------------------------------------------------------
//*  Input:  node(:)res = the residuals
//*
//* Output:  res_norm   = residual norms (L1, L2, Linf)
//* ------------------------------------------------------------------------------
//*
//* NOTE: It is not done here, but I advise you to keep the location of the
//*       maximum residual (L_inf).
//*
//********************************************************************************
void EulerSolver2D::Solver::residual_norm( EulerSolver2D::MainData2D& E2Ddata, Array2D<real>& res_norm_data) {
   //Array2D<real> res_norm_data(4,3); // residual norms L1, L2, Linfinity
   Array2D<real> residual(4,1);      // residual

   for (size_t i=0; i<E2Ddata.nq; ++i){
      res_norm_data(i,0) =  zero;
      res_norm_data(i,1) =  zero;
      res_norm_data(i,2) = - one;
   }

   
   //--------------------------------------------------------------------------------
   for (size_t i=0; i<E2Ddata.nnodes; ++i) {
   //--------------------------------------------------------------------------------
      for (size_t j=0; j<E2Ddata.nq; ++j){
         residual(j) = std::abs( (*E2Ddata.node[i].res)(j)/E2Ddata.node[i].vol );      //Divided residual
         res_norm_data(j,0) = res_norm_data(j,0)    + residual(j);               //L1   norm
         res_norm_data(j,1) = res_norm_data(j,1)    + residual(j)*residual(j);   //L2   norm
         res_norm_data(j,2) = std::max(res_norm_data(j,2), residual(j));         //Linf norm
      }
   //--------------------------------------------------------------------------------
   }//end loop nodes
   //--------------------------------------------------------------------------------
   for (size_t j=0; j<E2Ddata.nq; ++j){
      res_norm_data(j,0) = res_norm_data(j,0)/real(E2Ddata.nnodes);
      res_norm_data(j,1) = sqrt(res_norm_data(j,1)/real(E2Ddata.nnodes));
   }

}

// *******************************************************************************
//  -- Roe's Flux Function with entropy fix---
// 
//  P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
//  Schemes, Journal of Computational Physics, 43, pp. 357-372.
// 
//  NOTE: 3D version of this subroutine is available for download at
//        http://cfdbooks.com/cfdcodes.html
// 
//  ------------------------------------------------------------------------------
//   Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
//            primR(1:5) = right state (rhoR, uR, vR, pR)
//                njk(2) = Face normal (L -> R). Must be a unit vector.
// 
//  Output:    flux(1:5) = numerical flux
//                   wsn = half the max wave speed
//                         (to be used for time step calculations)
//  ------------------------------------------------------------------------------
// 
// *******************************************************************************
void EulerSolver2D::Solver::roe(EulerSolver2D::MainData2D& E2Ddata,
                                                      const Array2D<real>& primL,
                                                      const Array2D<real>& primR,
                                                      const Array2D<real>& njk,
                                                      Array2D<real>& flux,
                                                      real& wsn) {

   
   //Local variables
   real nx, ny;                  // Normal vector
   real mx, my;                  // Tangent vector: mx*nx+my*ny = 0
   real uL, uR, vL, vR;          // Velocity components.
   real rhoL, rhoR, pL, pR;      // Primitive variables.
   real unL, unR, umL, umR;      // Normal and tangent velocities
   real aL, aR, HL, HR;          // Speeds of sound.
   real RT,rho,u,v,H,a,un, um;   // Roe-averages
   real drho,dun,dum,dp;
   Array2D<real> LdU(4,1);       // Wave strenghs
   Array2D<real> ws(4,1);        // Wave speeds
   Array2D<real> Rv(4,4);        // right-eigevectors
   Array2D<real> fL(4,1);        // Flux left
   Array2D<real> fR(4,1);        // Flux right
   Array2D<real> diss(4,1);      // dissipation term
   Array2D<real> dws(4,1);       // User-specified width for entropy fix

   
   nx = njk(0);
   ny = njk(1);

   //Tangent vector (Actually, Roe flux can be implemented 
   // without any tangent vector. See "I do like CFD, VOL.1" for details.)
   mx = -ny;
   my =  nx;

   //Primitive and other variables.
   //  Left state
    rhoL = primL(0);
      uL = primL(1);
      vL = primL(2);
     unL = uL*nx+vL*ny;
     umL = uL*mx+vL*my;
      pL = primL(3);
      aL = sqrt(E2Ddata.gamma*pL/rhoL);
      HL = aL*aL/(E2Ddata.gamma-one) + half*(uL*uL+vL*vL);
   //  Right state
    rhoR = primR(0);
      uR = primR(1);
      vR = primR(2);
     unR = uR*nx+vR*ny;
     umR = uR*mx+vR*my;
      pR = primR(3);
      aR = sqrt(E2Ddata.gamma*pR/rhoR);
      HR = aR*aR/(E2Ddata.gamma-one) + half*(uR*uR+vR*vR);

   //First compute the Roe Averages
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL;
     u = (uL+RT*uR)/(one+RT);
     v = (vL+RT*vR)/(one+RT);
     H = (HL+RT* HR)/(one+RT);
     a = sqrt( (E2Ddata.gamma-one)*(H-half*(u*u+v*v)) );
    un = u*nx+v*ny;
    um = u*mx+v*my;

   //Wave Strengths
   drho = rhoR - rhoL;
     dp =   pR - pL;
    dun =  unR - unL;
    dum =  umR - umL;

   LdU(0) = (dp - rho*a*dun )/(two*a*a);
   LdU(1) = rho*dum;
   LdU(2) =  drho - dp/(a*a);
   LdU(3) = (dp + rho*a*dun )/(two*a*a);

   //Wave Speed
   ws(0) = abs(un-a);
   ws(1) = abs(un);
   ws(2) = abs(un);
   ws(3) = abs(un+a);

   //Harten's Entropy Fix JCP(1983), 49, pp357-393:
   // only for the nonlinear fields.
   dws(0) = fifth;
   if ( ws(0) < dws(0) ) ws(0) = half * ( ws(0)*ws(0)/dws(0)+dws(0) );
   dws(3) = fifth;
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) );

   //Right Eigenvectors
   Rv(0,0) = one;
   Rv(1,0) = u - a*nx;
   Rv(2,0) = v - a*ny;
   Rv(3,0) = H - un*a;

   Rv(0,1) = zero;
   Rv(1,1) = mx;
   Rv(2,1) = my;
   Rv(3,1) = um;

   Rv(0,2) = one;
   Rv(1,2) = u;
   Rv(2,2) = v;
   Rv(3,2) = half*(u*u+v*v);

   Rv(0,3) = one;
   Rv(1,3) = u + a*nx;
   Rv(2,3) = v + a*ny;
   Rv(3,3) = H + un*a;


   //Dissipation Term
   diss = zero;
   for (size_t i=0; i<E2Ddata.nq; ++i) {
      for (size_t j=0; j<E2Ddata.nq; ++j) {
         diss(i) += ws(j)*LdU(j)*Rv(i,j);
      }
   }

   //Compute the flux.
   fL(0) = rhoL*unL;
   fL(1) = rhoL*unL * uL + pL*nx;
   fL(2) = rhoL*unL * vL + pL*ny;
   fL(3) = rhoL*unL * HL;

   fR(0) = rhoR*unR;
   fR(1) = rhoR*unR * uR + pR*nx;
   fR(2) = rhoR*unR * vR + pR*ny;
   fR(3) = rhoR*unR * HR;

   flux = half * (fL + fR - diss);
   wsn = half*(abs(un) + a);  //Normal max wave speed times half

   // std::cout << "--------------------------- \n\n";
   // std::cout << "num_flux = \n";
   // flux.print();
   // std::cout << "wsn = " << wsn << "\n";
   // std::cout << "--------------------------- \n\n";
}  


//*****************************************************************************
//* -- Rotated-RHLL Flux Function ---
//*
//* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
//* Resolving, Rotated-Hybrid Riemann Solvers,
//* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
//*
//* Robust Riemann solver for nonlinear instability (carbuncle).
//*
//* NOTE: 3D version of this subroutine is available for download at
//*       http://cfdbooks.com/cfdcodes.html
//*
//* ------------------------------------------------------------------------------
//*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
//*           primR(1:5) = right state (rhoR, uR, vR, pR)
//*               njk(2) = Face normal (L -> R). Must be a unit vector.
//*
//* Output:    flux(1:5) = numerical flux
//*                  wsn = half the max wave speed
//*                        (to be used for time step calculations)
//* ------------------------------------------------------------------------------
//*
//*****************************************************************************
void EulerSolver2D::Solver::rotated_rhll(EulerSolver2D::MainData2D& E2Ddata,
                                                      const Array2D<real>& primL,
                                                      const Array2D<real>& primR,
                                                      const Array2D<real>& njk,
                                                      Array2D<real>& flux,
                                                      real wsn) {
   }

//********************************************************************************
//* -- vanAlbada Slope Limiter Function--
//*
//* 'A comparative study of computational methods in cosmic gas dynamics', 
//* Van Albada, G D, B. Van Leer and W. W. Roberts, Astronomy and Astrophysics,
//* 108, p76, 1982
//*
//* ------------------------------------------------------------------------------
//*  Input:   da, db     : two differences
//*
//* Output:   va_limiter : limited difference
//* ------------------------------------------------------------------------------
//*
//********************************************************************************
Array2D<real> EulerSolver2D::Solver::va_slope_limiter(EulerSolver2D::MainData2D& E2Ddata,
                                                      const Array2D<real>& da, 
                                                      const Array2D<real>& db, 
                                                      const real h) {
   Array2D<real> va_slope_limiter_data(4,1);
   real eps2;
   
   eps2 = pow(3, (0.3 * h) );

   // TLM todo: implement a vector (array) version of copysign and make the code below a 1 liner. 
   for (size_t i=0; i<E2Ddata.nq; ++i){
      va_slope_limiter_data(i) = half*( copysign(1.0,da(i)*db(i)) + one ) * 
           ( (db(i)*db(i) + eps2)*da(i) + (da(i)*da(i) + eps2)*db(i) ) /
           (da(i)*da(i) + db(i)*da(i) + two*eps2);
   }
   return va_slope_limiter_data;
}


//********************************************************************************
//* Compute U from W
//*
//* ------------------------------------------------------------------------------
//*  Input:  w =    primitive variables (rho,     u,     v,     p)
//* Output:  u = conservative variables (rho, rho*u, rho*v, rho*E)
//* ------------------------------------------------------------------------------
//* 
//********************************************************************************
Array2D<real>  EulerSolver2D::Solver::w2u(const Array2D<real>& w,
                                             EulerSolver2D::MainData2D& E2Ddata) {
   Array2D<real> u(4,1);

   u(0) = w(0);
   u(1) = w(0)*w(1);
   u(2) = w(0)*w(2);
   u(3) = w(3)/(E2Ddata.gamma-one)+half*w(1)*(w(2)*w(2)+w(3)*w(3));

   return u;

} // end function w2u
//--------------------------------------------------------------------------------

//********************************************************************************
//* Compute U from W
//*
//* ------------------------------------------------------------------------------
//*  Input:  u = conservative variables (rho, rho*u, rho*v, rho*E)
//* Output:  w =    primitive variables (rho,     u,     v,     p)
//* ------------------------------------------------------------------------------
//* 
//********************************************************************************
Array2D<real>  EulerSolver2D::Solver::u2w(const Array2D<real>& u,
                                             EulerSolver2D::MainData2D& E2Ddata) {

   Array2D<real> w(4,1);

   w(0) = u(0);
   w(1) = u(1)/u(0);
   w(2) = u(2)/u(0);
   w(3) = (E2Ddata.gamma-one)*( u(3) - half*w(0)*(w(1)*w(1)+w(2)*w(2)) );

   return w;

}//end function u2w
//--------------------------------------------------------------------------------



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
   // use edu2d_constants   , only : p2, one, two
   // use edu2d_my_main_data, only : nnodes, node

   int       i, ix, iy, ivar;
   std::string grad_type_temp;
   real error_max_wx, error_max_wy, x, y;
   real x_max_wx, y_max_wx, x_max_wy, y_max_wy, wx, wxe, wy, wye;
   real a0, a1, a2, a3, a4, a5;

   ix = 0;
   iy = 1;

   // We only use w(0) for this test.
   ivar = 0;

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
   grad_type_temp = "linear";
   //cout << " grad_type_temp =  " << grad_type_temp << endl;
   compute_gradient_nc(E2Ddata, ivar, grad_type_temp);
   cout << " ivar = " << ivar << endl;

//  (3). Compute the relative errors (L_infinity)

   cout << "- Computing the relative errors (L_infinity)..\n";
   error_max_wx = -one;
   error_max_wy = -one;


   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      //cout << " (*E2Ddata.node[i].gradw)(ivar,ix) = " << (*E2Ddata.node[i].gradw)(ivar,ix) << endl;
      //cout << " (*E2Ddata.node[i].gradw)(ivar,iy) = " << (*E2Ddata.node[i].gradw)(ivar,iy) << endl;
      error_max_wx = max( std::abs( (*E2Ddata.node[i].gradw)(ivar,ix) - one )/one, error_max_wx );
      error_max_wy = max( std::abs( (*E2Ddata.node[i].gradw)(ivar,iy) - two )/two, error_max_wy );
   }

   cout << " Max relative error in wx =  " << error_max_wx << "\n";
   cout << " Max relative error in wy =  " << error_max_wy << "\n";
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
   grad_type_temp = "quadratic2";
   compute_gradient_nc(E2Ddata,ivar,grad_type_temp);

//  (3). Compute the relative errors (L_infinity)

   cout << "- Computing the relative errors (L_infinity)..\n";
   error_max_wx = -one;
   error_max_wy = -one;
   // find node i's of max grad error
   size_t intmax_x = 0;
   size_t intmax_y = 0;
   //loop nnodes
   for (size_t i = 0; i < E2Ddata.nnodes; i++) {
      x = E2Ddata.node[i].x;
      y = E2Ddata.node[i].y;

      if ( std::abs( (*E2Ddata.node[i].gradw)(ivar,ix) - 
            (a1+2.0*a3*x+a4*y) )/(a1+2.0*a3*x+a4*y) > error_max_wx )  {
         wx  = (*E2Ddata.node[i].gradw)(ivar,ix);
         wxe = a1 + 2.0*a3*x + a4*y;
         error_max_wx = std::abs( wx - wxe )/wxe;
         x_max_wx = x;
         y_max_wx = y;
         intmax_x = i;
      }

      if ( std::abs( (*E2Ddata.node[i].gradw)(ivar,iy) - 
            (a2+2.0*a5*y+a4*x) )/(a2+2.0*a5*y+a4*x) > error_max_wy )  {
         wy  = (*E2Ddata.node[i].gradw)(ivar,iy);
         wye = a2 + 2.0*a5*y + a4*x;
         error_max_wy = std::abs( wy - wye )/wye;
         x_max_wy = x;
         y_max_wy = y;
         intmax_y = i;
      }

   }//end do

  cout << " Max relative error in wx = " <<  error_max_wx <<  " at (x,y) = (" 
                                          <<  x_max_wx << " , " << y_max_wx << ")\n";
  cout << "   At this location, LSQ ux = " <<  wx <<  ": Exact ux = " <<  wxe << "\n";
  cout << " Max relative error in wy = " <<  error_max_wy <<  " at (x,y) = (" 
                                          <<  x_max_wy << " , " << y_max_wy << ")\n";
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
void EulerSolver2D::Solver::compute_gradient_nc(
                              EulerSolver2D::MainData2D& E2Ddata,
                              int ivar, std::string grad_type) {

   //use edu2d_my_main_data, only : node, nnodes

   //integer, intent(in) :: ivar
   //std::string grad_type

   int in;

   //cout << "computing gradient nc type " << grad_type << endl;

   if (trim(grad_type) == "none") {
      cout << "trim(grad_type) == none, return " << endl;
      return;
   }
   //   else {
   //      cout << "trim(grad_type) == " << trim(grad_type) << " full =  " << grad_type << endl;
   //   }

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
   //cout << "Compute LSQ Gradients at all nodes. " << endl;
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

         //cout << "(trim(grad_type) == quadratic2) " << endl;
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
void EulerSolver2D::Solver::lsq_gradients_nc(
   EulerSolver2D::MainData2D& E2Ddata, int inode, int ivar) {

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
      
      // if (inode<E2Ddata.maxit-0) {
      //    cout << " ax = " << ax << endl;
      //    cout << " ay = " << ay << endl;
      // }

      // cout << " (*E2Ddata.node[inode].lsq2x2_cx)(in) = "
      //          << (*E2Ddata.node[inode].lsq2x2_cx)(in) << endl;
      // cout << " (*E2Ddata.node[inode].lsq2x2_cy)(in) = "
      //          << (*E2Ddata.node[inode].lsq2x2_cy)(in) << endl;

   }

   (*E2Ddata.node[inode].gradw)(ivar,ix) = ax;  //<-- du(ivar)/dx
   (*E2Ddata.node[inode].gradw)(ivar,iy) = ay;  //<-- du(ivar)/dy
   
   // if (inode >= E2Ddata.nnodes-1) {
   //    std::cout << "ix = " << ix << "\n";
   //    std::cout << "iy = " << iy << "\n";
   //    std::cout << "inode = " << inode << "\n";
   // }


}// end lsq_gradients_nc
//--------------------------------------------------------------------------------



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

   ii = -1;

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

         // if (inode<E2Ddata.maxit-1) {
         //    cout << "(*E2Ddata.node[inode].lsq5x5_cx)(ii) = " << (*E2Ddata.node[inode].lsq5x5_cx)(ii) << endl;
         //    cout << "(*E2Ddata.node[inode].lsq5x5_cy)(ii) = " << (*E2Ddata.node[inode].lsq5x5_cy)(ii) << endl;
         // }
         // if (inode<E2Ddata.maxit-0) {
         //    cout << " = " << ax << endl;
         //    cout << " = " << ay << endl;
         // }

      } // end loop nghbr_nghbr

   } // end loop nghbr

   (*E2Ddata.node[inode].gradw)(ivar,ix) = ax;  //<-- dw(ivar)/dx;
   (*E2Ddata.node[inode].gradw)(ivar,iy) = ay;  //<-- dw(ivar)/dy;

   // if (inode >= E2Ddata.nnodes-1) {
   //    std::cout << "ix = " << ix << "\n";
   //    std::cout << "iy = " << iy << "\n";
   //    std::cout << "inode = " << inode << "\n";
   // }

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
void EulerSolver2D::Solver::lsq01_2x2_coeff_nc(
   EulerSolver2D::MainData2D& E2Ddata , int inode) {


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
   // local_lsq_inverse(1,0) = -a(0,1)/det;
   // local_lsq_inverse(1,1) =  a(0,0)/det;

   local_lsq_inverse(0,0) =  a(1,1)/det;
   local_lsq_inverse(0,1) = -a(1,0)/det;
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

      (*E2Ddata.node[inode].lsq2x2_cx)(k)  = \
                                   local_lsq_inverse(ix,0)*w2dvar*dx \
                                 + local_lsq_inverse(ix,1)*w2dvar*dy;

      (*E2Ddata.node[inode].lsq2x2_cy)(k)  = \
                                   local_lsq_inverse(iy,0)*w2dvar*dx \
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
   cout << "gradient_weight  = " << trim(E2Ddata.gradient_weight) << endl;
   
   
   // Step 1

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
               if (i < E2Ddata.maxit-1 && ell < E2Ddata.maxit-1) {
                  cout << " --------------E2Ddata.node[in].nghbr)(ell) == i " << endl;
                  cout << " (*E2Ddata.node[in].nghbr)(ell) = " << (*E2Ddata.node[in].nghbr)(ell) << " i = " << i << endl;
                  cout << " dx = " << dx << endl;
                  cout << " dy = " << dy << endl;
                  cout << "  "<< endl;
               }

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

            if (i < E2Ddata.maxit-1 & ell < E2Ddata.maxit-1) {
               cout << " ---------TOP---------- " << endl;
               cout << " E2Ddata.node[in].x = " << E2Ddata.node[in].x << endl;
               cout << " E2Ddata.node[in].y = " << E2Ddata.node[in].y << endl;
               cout << " E2Ddata.node[i].x = " << E2Ddata.node[i].x << endl;
               cout << " E2Ddata.node[i].y = " << E2Ddata.node[i].y << endl;
               cout << " dx = " << dx << endl;
               cout << " dy = " << dy << endl;
               cout << "  "<< endl;
            }

            w2 = lsq_weight(E2Ddata, dx, dy);
            w2 = w2*w2;


            if (i < E2Ddata.maxit-1) {
               cout << " ---------MIDDLE---------- " << endl;
               cout << " dx = " << dx << endl;
               cout << " dy = " << dy << endl;
               cout << " w2 = " << w2 << endl;
               //cout << " half = " << half << endl;
               cout << "  "<< endl;
            }
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

      if (i < E2Ddata.maxit-1) {
         cout << " a = " << endl;
         print(a);
      }

      dummy1 = zero;
      dummy2 = zero;
      //gewp_solve(a,dummy1,dummy2,ainv,istat, 5);
      //ainv = 0.;// 
      //cout << " invert i = " << i << endl;
      Array2D ca = a;
      MatrixXd me(5,5),ime(5,5);
      for (size_t ie = 0; ie < 5; ie++) {
         for (size_t je = 0; je < 5; je++) {
            me(ie,je) = a(ie,je);
         }
      }
      ime = me.inverse();
      //ainv = a.invert(); // destructive solve //GSinv(a,dummy1,dummy2);
      //ainv = a.inverse();
      for (size_t ie = 0; ie < 5; ie++) {
         for (size_t je = 0; je < 5; je++) {
            ainv(ie,je) = ime(ie,je);
         }
      }

      if (i < E2Ddata.maxit-1) {
         std::ostringstream out;
         out << " ---------BOTTOM---------- " << endl;
         out << " dx = " << dx << endl;
         out << " dy = " << dy << endl;
         out << " w2 = " << w2 << endl;
         E2Ddata.write_diagnostic(out, "log/out.dat");
      }

      if (i < E2Ddata.maxit-1) {
         std::ostringstream out;
         out << " a = \n" << endl;
         a.print(out);
         out << "\n ainv = \n" << endl;
         ainv.print(out);
         //ainv.print();
         E2Ddata.write_diagnostic(out, "log/out.dat");
      }
      if (i < E2Ddata.maxit-1) {
         std::ostringstream out;
         cout << " check 1" << endl;
         Array2D cka1 = matmul(ainv, ca);
         cout << "ainv * a = " << endl;
         print (cka1);

         Array2D cka = matmul(ainv, ca);
         cout << " check 2" << endl;
         Array2D cka2 = matmul(ca, ainv);
         cout << "a * ainv = " << endl;
         print (cka2);
      }
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
real EulerSolver2D::Solver::lsq_weight(
   EulerSolver2D::MainData2D& E2Ddata, real dx, real dy) {

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

      // cout << "INVERSE DISTANCE  = " << endl;
      // std::exit(0);


      distance = std::sqrt(dx*dx + dy*dy);

      real val = std::pow(distance, E2Ddata.gradient_weight_p);
      if (val < 1.e-6) {
         cout << "ERROR: distance is 0" << endl;
         std::exit(0);
      }

      lsq_weight = one / val;
   }
   
   //cout << " lsq_weight = " << lsq_weight << endl;

   return lsq_weight;
} //end lsq_weight






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
