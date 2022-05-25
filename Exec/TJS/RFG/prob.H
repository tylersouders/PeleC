#ifndef _PROB_H_
#define _PROB_H_

//#include <AMReX_Print.H>
//#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

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
  amrex::Real pert = (1.0 + ((((rand()%100)/99.) - 0.5) * PeleC::h_prob_parm_device->init_pert_mag));
  amrex::Real u_pert = PeleC::h_prob_parm_device->u_in * pert;

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  //amrex::Real u_profile = (8.0 / (prob_hi[1] * prob_hi[1])) * ProbParm::u_in *
  //  (prob_hi[1] - y) * (y - prob_lo[1]);

  // Set the state
  state(i, j, k, URHO) = PeleC::h_prob_parm_device->rho_in;
  state(i, j, k, UMX) = u_pert * PeleC::h_prob_parm_device->rho_in;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
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

  if (1) {
    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

    if (sgn == 1) {
      // inflow
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();

      // turbulence inflow

      // initialize variables
      const amrex::Real tau = PeleC::h_prob_parm_device->turb_length_scale / PeleC::h_prob_parm_device->turb_velocity;
      const amrex::Real t_hat = 2.0 * constants::PI() * time / tau;
      amrex::Real x_hat[AMREX_SPACEDIM];
      amrex::Real turb_u[AMREX_SPACEDIM] = {0.0};

      amrex::Real isotropic_factor = 1.0 / std::sqrt(1.0/3.0);
      // convert to turbulence terms
      for(int i = 0; i < AMREX_SPACEDIM; i++) {
       x_hat[i] = 2.0 * constants::PI() * x[i] / PeleC::h_prob_parm_device->turb_length_scale;
     }
      // assume inflow is purely x-dimensional.
     x_hat[0] -= (2.0 * constants::PI() * PeleC::h_prob_parm_device->u_in * time)
     / PeleC::h_prob_parm_device->turb_length_scale;
     
     for(int i = 0; i < AMREX_SPACEDIM; i++) {
      for(int n = 0; n < PeleC::h_prob_parm_device->turb_num_modes; n++) {

        amrex::Real arg = isotropic_factor * (PeleC::h_prob_parm_device->d[0][n] * x_hat[0]
          + PeleC::h_prob_parm_device->d[1][n] * x_hat[1]
          + PeleC::h_prob_parm_device->d[2][n] * x_hat[2]);
        arg += PeleC::h_prob_parm_device->omega[n] * t_hat;
        turb_u[i] += PeleC::h_prob_parm_device->p[i][n] * cos(arg)
        + PeleC::h_prob_parm_device->q[i][n] * sin(arg);
      }
	// Again, the sqrt(1/3) is an attempt at isotropy. Change if implementing anisotropy.
      
      turb_u[i] *= std::sqrt(2.0/PeleC::h_prob_parm_device->turb_num_modes);
      turb_u[i] *= PeleC::h_prob_parm_device->turb_velocity * std::sqrt(1.0/3.0); 
    }
    
    rho = PeleC::h_prob_parm_device->rho_in;
      //u = (8.0 / (prob_hi[1] * prob_hi[1])) * ProbParm::u_in * (prob_hi[1] - x[1]) * (x[1] - prob_lo[1]);
    u = PeleC::h_prob_parm_device->u_in + turb_u[0];
    v += turb_u[1];
    if (AMREX_SPACEDIM > 2) {
      w += turb_u[2];
    }
    eint =  PeleC::h_prob_parm_device->e_in;
    T = PeleC::h_prob_parm_device->T_in;
    for (int n = 0; n < NUM_SPECIES; n++)
     s_ext[UFS + n] = PeleC::h_prob_parm_device->rho_in * PeleC::h_prob_parm_device->massfrac[n];

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
  for (int n = 0; n < NUM_SPECIES; n++)
   s_ext[UFS + n] = s_int[UFS+n];
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
