#ifndef _PROB_H_
#define _PROB_H_

//#include <AMReX_Print.H>
//#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Random.H>

#include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "Constants.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Get random val for initial x-momentum.
  // Consider changing seed to be non-deterministic
  amrex::Real pert = 1.0 + ((amrex::Random() - 0.5) * prob_parm.init_pert_mag);
  amrex::Real u_pert = prob_parm.u_in * pert;
  amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  //amrex::Real u_profile = (8.0 / (prob_hi[1] * prob_hi[1])) * ProbParm::u_in *
  //  (prob_hi[1] - y) * (y - prob_lo[1]);

  // If flag is enabled, fill entire domain with turbulence using Huang's method
  // if (prob_parm.init_turb_fill) {

  //   // Find cartesian location based on indices and dx (computed above)
  //   amrex::Real loc_factor = 2.0 * constants::PI() / 
  //   (prob_parm.Ls_prefactor * prob_parm.turb_length_scale);
  //   amrex::Real x_hat = (prob_lo[0] + (i + 0.5) * dx[0]) * loc_factor;
  //   amrex::Real y_hat = (prob_lo[1] + (j + 0.5) * dx[1]) * loc_factor;
  //   amrex::Real z_hat = (prob_lo[2] + (k + 0.5) * dx[2]) * loc_factor;

  //   for(int m = 0; m < prob_parm.turb_num_modes; m++) {

  //     for(int n = 0; n < prob_parm.sampling_number; n++) {

  //       amrex::Real arg = prob_parm.k[0][n][m] * x_hat
  //       + prob_parm.k[1][n][m] * y_hat
  //       + prob_parm.k[2][n][m] * z_hat
  //       + prob_parm.omega[n][m] * 0.0; //time = 0.0 becase init

  //       // Compute turbulence
  //       for(int i = 0; i < 3; i++){
  //         turb_u[i] += prob_parm.p[i][n][m] * cos(arg)
  //         + prob_parm.q[i][n][m] * sin(arg);
  //       }
  //     }
  //   }
  // }

  // Set the state
  state(i, j, k, URHO) = prob_parm.rho_in;
  if (prob_parm.init_turb_fill) {
    state(i, j, k, UMX) = turb_u[0] * prob_parm.rho_in;
    state(i, j, k, UMY) = turb_u[1] * prob_parm.rho_in;
    state(i, j, k, UMZ) = turb_u[2] * prob_parm.rho_in;
  } else {
    state(i, j, k, UMX) = u_pert * prob_parm.rho_in;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
  }
  state(i, j, k, UEINT) = prob_parm.rho_in * prob_parm.e_in;
  state(i, j, k, UEDEN) = prob_parm.rho_in * 
  (prob_parm.e_in + 0.5 * (u_pert * u_pert));
  state(i, j, k, UTEMP) = prob_parm.T_in;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = prob_parm.rho_in * prob_parm.massfrac[n];
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParmDevice const& prob_parm)
{
  // Coefficient to optionally ramp up inflow velocity
  amrex::Real t_coef = amrex::min(1.0, time / 0.00005);

  if (1) {

    // Initialize the state variables as zeros
    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

    if (sgn == 1) {
      // Inlet/Inflow
      // amrex::Print() << "starting bcnormal. x = " << x[0] << " " << x[1] << " " << x[2] << " " << '\n';

      // Save the time for debug
      // prob_parm.timechk = time;

      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real z_len = prob_hi[2] - prob_lo[2];
      const amrex::Real rho0 = prob_parm.rho_in;

      // const bool periodic_locationcheck = (x[2] >= 0.03 * z_len && x[2] <= 0.97 * z_len);

      // if(prob_parm.do_turb_inlet && periodic_locationcheck) {
      if(prob_parm.do_turb_inlet) {

        // Compute turbulence inflow as described in Huang, 2010 Generalized turbulent inflow

        // initialize variables
        amrex::Real x_hat[AMREX_SPACEDIM]  = {0.0};
        amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

        // Normalize position values
        amrex::Real loc_factor = 2.0 * constants::PI() / 
        (prob_parm.Ls_prefactor * prob_parm.turb_length_scale);
        for(int i = 0; i < AMREX_SPACEDIM; i++) {
          x_hat[i] = x[i] * loc_factor;
        }

        // Correct for velocity at inflow, assuming inflow is approximately 
        // purely x-dimensional (u >> v, w)
        // x_hat[0] -= (prob_parm.u_in * time) * loc_factor;

        //amrex::Print() << "bcnormal check #1 (before loop) " << '\n';

        for(int m = 0; m < prob_parm.turb_num_modes; m++) {

          for(int n = 0; n < prob_parm.sampling_number; n++) {

            // amrex::Real arg = prob_parm.k[0][n][m] * x_hat[0]
            // + prob_parm.k[1][n][m] * x_hat[1]
            // + prob_parm.k[2][n][m] * x_hat[2]
            // + prob_parm.omega[n][m] * time;

            amrex::Real arg = prob_parm.k[0*(prob_parm.turb_num_modes*prob_parm.sampling_number)
            + m*prob_parm.sampling_number + n] * x_hat[0]
            + prob_parm.k[1*(prob_parm.turb_num_modes*prob_parm.sampling_number) 
              + m*prob_parm.sampling_number + n] * x_hat[1]
            + prob_parm.k[2*(prob_parm.turb_num_modes*prob_parm.sampling_number) 
              + m*prob_parm.sampling_number + n] * x_hat[2]
            + prob_parm.omega[prob_parm.turb_num_modes
              * prob_parm.sampling_number + n] * time;

            // Compute turbulence perturbation
            for(int i = 0; i < 3; i++){
              // turb_u[i] += prob_parm.p[i][n][m] * cos(arg)
              // + prob_parm.q[i][n][m] * sin(arg);

              turb_u[i] += prob_parm.p[i*(prob_parm.turb_num_modes*prob_parm.sampling_number)
                + m*prob_parm.sampling_number + n] * cos(arg)
              + prob_parm.q[i*(prob_parm.turb_num_modes*prob_parm.sampling_number)
                + m*prob_parm.sampling_number + n] * sin(arg);
            }
          }
        }

        u = turb_u[0]+prob_parm.u_in;
        v = turb_u[1];
        w = turb_u[2];

        // Compute adjustment for rho to account for velocity increase
        //amrex::Real temp_ratio = prob_parm.u_in / std::sqrt(
        //  std::pow(u, 2.0) + std::pow(v, 2.0) + std::pow(w, 2.0));
        //rho = rho0 * std::pow(temp_ratio, 2.0);
        rho = rho0;

      } else {

        u = prob_parm.u_in;
        v = 0.0;
        w = 0.0;
        rho = rho0;

      }

      eint =  prob_parm.e_in;
      T = prob_parm.T_in;
      for (int n = 0; n < NUM_SPECIES; n++)
        s_ext[UFS + n] = prob_parm.rho_in * prob_parm.massfrac[n];
        //s_ext[UFS + n] = rho * prob_parm.massfrac[n];
      // amrex::Print() << "bcnormal check #3 (finished) " << '\n';

    } else if (sgn == -1) {
      // outflow
      // Following Blazek p 279, eq. 8.23

      // Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_int = s_int[URHO];
      const amrex::Real u_int = s_int[UMX] / rho_int;
      const amrex::Real v_int = s_int[UMY] / rho_int;
      const amrex::Real w_int = s_int[UMZ] / rho_int;
      const amrex::Real T_int = s_int[UTEMP];
      amrex::Real p_int = 0.0, cs_int = 0.0;

      // get species massfracs from interior points.
      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac_int[n] = s_int[UFS + n] / rho_int;

      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2P(rho_int, T_int, massfrac_int.begin(), p_int);
      eos.RTY2Cs(rho_int, T_int, massfrac_int.begin(), cs_int);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p_in;
      const amrex::Real rhob = s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
      const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
      const amrex::Real vb = v_int;
      const amrex::Real wb = w_int;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_int) / (xb - xd) * (x[0] - xd) + rho_int;
      const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;

      eos.RYP2E(rho, massfrac_int.begin(), p, eint);
      eos.EY2T(eint, massfrac_int.begin(), T);

      u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
      v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
      w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;

      // Outflow-specific
      for (int n = 0; n < NUM_SPECIES; n++) {
       s_ext[UFS + n] = s_int[UFS+n];
     }
      // amrex::Print() << "outflow bcnormal check (finished) " << '\n';
   }

   s_ext[URHO] = rho;
   s_ext[UMX] = rho * u;
   s_ext[UMY] = rho * v;
   s_ext[UMZ] = rho * w;
   s_ext[UEINT] = rho * eint;
   s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
   s_ext[UTEMP] = T;

 }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
