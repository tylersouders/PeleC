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
  amrex::Real pert = 1.0 + ((amrex::Random() - 0.5) * PeleC::h_prob_parm_device->init_pert_mag);
  amrex::Real u_pert = PeleC::h_prob_parm_device->u_in * pert;
  amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  //amrex::Real u_profile = (8.0 / (prob_hi[1] * prob_hi[1])) * ProbParm::u_in *
  //  (prob_hi[1] - y) * (y - prob_lo[1]);

  // If flag is enabled, fill entire domain with turbulence using Huang's method
  if (PeleC::h_prob_parm_device->init_turb_fill) {

    amrex::Print() << "Initializing with Turbulence Fill!!" << std::endl;

    // Find cartesian location based on indices and dx (computed above)
    amrex::Real loc_factor = 2.0 * constants::PI() / 
    (PeleC::h_prob_parm_device->Ls_prefactor * PeleC::h_prob_parm_device->turb_length_scale);
    amrex::Real x_hat = (prob_lo[0] + (i + 0.5) * dx[0]) * loc_factor;
    amrex::Real y_hat = (prob_lo[1] + (j + 0.5) * dx[1]) * loc_factor;
    amrex::Real z_hat = (prob_lo[2] + (k + 0.5) * dx[2]) * loc_factor;

    for(int m = 0; m < PeleC::h_prob_parm_device->turb_num_modes; m++) {

      for(int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {

        amrex::Real arg = PeleC::h_prob_parm_device->k[0][n][m] * x_hat
        + PeleC::h_prob_parm_device->k[1][n][m] * y_hat
        + PeleC::h_prob_parm_device->k[2][n][m] * z_hat
        + PeleC::h_prob_parm_device->omega[n][m] * 0.0; //time = 0.0 becase init
        
        // Compute turbulence
        for(int i = 0; i < 3; i++){
          turb_u[i] += PeleC::h_prob_parm_device->p[i][n][m] * cos(arg)
          + PeleC::h_prob_parm_device->q[i][n][m] * sin(arg);
        }
      }
    }
  }

  // Set the state
  state(i, j, k, URHO) = PeleC::h_prob_parm_device->rho_in;
  if (PeleC::h_prob_parm_device->init_turb_fill) {
    state(i, j, k, UMX) = turb_u[0] * PeleC::h_prob_parm_device->rho_in;
    state(i, j, k, UMY) = turb_u[1] * PeleC::h_prob_parm_device->rho_in;
    state(i, j, k, UMZ) = turb_u[2] * PeleC::h_prob_parm_device->rho_in;
  } else {
    state(i, j, k, UMX) = u_pert * PeleC::h_prob_parm_device->rho_in;
    state(i, j, k, UMY) = 0.0;
    state(i, j, k, UMZ) = 0.0;
  }
  state(i, j, k, UEINT) = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->e_in;
  state(i, j, k, UEDEN) = PeleC::h_prob_parm_device->rho_in * 
  (PeleC::h_prob_parm_device->e_in + 0.5 * (u_pert * u_pert));
  state(i, j, k, UTEMP) = PeleC::h_prob_parm_device->T_in;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->massfrac[n];
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

  // Initialize all state variables before any conditional statement
  amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, 
  eint = 0.0, T = 0.0, P = 0.0;
  amrex::Real umx = 0.0, umy = 0.0, umz = 0.0, ueint = 0.0, 
  ueden = 0.0, utemp = 0.0, urho = 0.0;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

  if (1) {

    if (sgn == 1) {
      // Inlet/Inflow
      // amrex::Print() << "starting bcnormal. x = " << x[0] << " " << x[1] << " " << x[2] << " " << '\n';

      // Save the time for debug
      PeleC::h_prob_parm_device->timechk = time;

      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();

      if(PeleC::h_prob_parm_device->do_turb_inlet) {

        // Compute turbulence inflow as described in Huang, 2010 Generalized turbulent inflow

        // initialize variables
        amrex::Real x_hat[AMREX_SPACEDIM]  = {0.0};
        amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

        // Normalize position values
        amrex::Real loc_factor = 2.0 * constants::PI() / 
        (PeleC::h_prob_parm_device->Ls_prefactor * PeleC::h_prob_parm_device->turb_length_scale);
        for(int i = 0; i < AMREX_SPACEDIM; i++) {
          x_hat[i] = x[i] * loc_factor;
        }

        // Correct for velocity at inflow, assuming inflow is approximately 
        // purely x-dimensional (u >> v, w)
        // x_hat[0] -= (PeleC::h_prob_parm_device->u_in * time) * loc_factor;

        //amrex::Print() << "bcnormal check #1 (before loop) " << '\n';

        for(int m = 0; m < PeleC::h_prob_parm_device->turb_num_modes; m++) {

          for(int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {

            amrex::Real arg = PeleC::h_prob_parm_device->k[0][n][m] * x_hat[0]
            + PeleC::h_prob_parm_device->k[1][n][m] * x_hat[1]
            + PeleC::h_prob_parm_device->k[2][n][m] * x_hat[2]
            + PeleC::h_prob_parm_device->omega[n][m] * time;

            // Compute turbulence perturbation
            for(int i = 0; i < 3; i++){
              turb_u[i] += PeleC::h_prob_parm_device->p[i][n][m] * cos(arg)
              + PeleC::h_prob_parm_device->q[i][n][m] * sin(arg);
            }
          }
        }

        u = turb_u[0]+PeleC::h_prob_parm_device->u_in;
        v = turb_u[1];
        w = turb_u[2];

      } else {

        u = PeleC::h_prob_parm_device->u_in;
        v = 0.0;
        w = 0.0;

      }

      // After computing velocity, need to use Blazek pg. 279 for subsonic inlet
      // Get spatial information for extrapolation
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real x_ext = prob_lo[0] - 0.5 * dx[0];
      const amrex::Real x_b = prob_lo[0];
      const amrex::Real x_int = prob_lo[0] + 0.5 * dx[0];

      // Interior Values (point d)
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

      // We know the goal values
      const amrex::Real T_b = PeleC::h_prob_parm_device->T_in;
      const amrex::Real p_b = PeleC::h_prob_parm_device->p_in;
      const amrex::Real u_b = u; // Values computed above
      const amrex::Real v_b = v;
      const amrex::Real w_b = w;

      // Extrapolate Ghost Values
      const amrex::Real T_ext = (T_b - T_int) / (x_b - x_int) * (x_ext - x_int) + T_int;
      const amrex::Real p_ext = (p_b - p_int) / (x_b - x_int) * (x_ext - x_int) + p_int;
      const amrex::Real u_ext = (u_b - u_int) / (x_b - x_int) * (x_ext - x_int) + u_int;
      const amrex::Real v_ext = (v_b - v_int) / (x_b - x_int) * (x_ext - x_int) + v_int;
      const amrex::Real w_ext = (w_b - w_int) / (x_b - x_int) * (x_ext - x_int) + w_int;

      // Use eos to get eint and rho at ghost cell
      auto eos = pele::physics::PhysicsType::eos();
      eos.PYT2RE(p_ext, massfrac_int.begin(), T_ext, rho, eint);

      // Get interior (point d) values
      for (int n = 0; n < NUM_SPECIES; n++)
        s_ext[UFS + n] = PeleC::h_prob_parm_device->rho_in * 
      PeleC::h_prob_parm_device->massfrac[n];

      // Get variables for return
      urho = rho;
      umx = rho * u_ext;
      umy = rho * v_ext;
      umz = rho * w_ext;
      ueint = rho * eint;
      utemp = T_ext;
      ueden = rho * (eint * 0.5 * (u_ext * u_ext + v_ext * v_ext + w_ext * w_ext));

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
      const amrex::Real pb = PeleC::h_prob_parm_device->p_in;
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

      // Get variables for return
     urho = rho;
     umx = rho * u;
     umy = rho * v;
     umz = rho * w;
     ueint = rho * eint;
     utemp = T;
     ueden = rho * (eint * 0.5 * (u * u + v * v + w * w));

   }

   s_ext[URHO] = urho;
   s_ext[UMX] = umx;
   s_ext[UMY] = umy;
   s_ext[UMZ] = umz;
   s_ext[UEINT] = ueint;
   s_ext[UEDEN] = ueden;
   s_ext[UTEMP] = utemp;

 }

}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
