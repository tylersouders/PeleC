#ifndef _PROB_H_
#define _PROB_H_

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
  // Geometric Values
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  
  // Variables to Define Localized Fuel Concentrations
  amrex::Real rho = 0.0;
  amrex::Real eint = 0.0;

  //auto eos = pele::physics::PhysicsType::eos();
  //eos.X2Y(molefrac_alt.begin(), massfrac_alt.begin());
  //eos.PYT2RE(
  //   prob_parm.p_input, massfrac_alt.begin(),
  //   prob_parm.T_input, rho_alt, eint_alt);

  // Add u perturbations
  amrex::Real perturb = (1.0 + ((((rand()%100)/99.0)-0.5)*prob_parm.init_pert_mag));
  amrex::Real u_perturb = prob_parm.u_input * perturb;

  // Initial Pressure Value (is this redundant?)
  const amrex::Real p = prob_parm.p_input;

  // Initialize Energy from Mass composition
  //amrex::Real rho = 0.0, eint = 0.0;
  //amrex::Real massfrac[NUM_SPECIES] = {0.0};
  //for (int n = 0; n < NUM_SPECIES; n++) {
  //	massfrac[n]	= prob_parm.massfrac[n];
  //}

  //auto eos = pele::physics::PhysicsType::eos();
  //eos.PYT2RE(p, massfrac, prob_parm.T_input, rho, eint);
  //amrex::Real rho = prob_parm.rho_input;
  //amrex::Real eint = prob_parm.eint_input;

  // Initialize domain state, using fuel_y_limit
  if (std::abs(y) < prob_parm.fuel_y_limit && x < prob_parm.fuel_x_limit) {
    rho = prob_parm.rho_input;
    state(i, j, k, UEINT) = rho * prob_parm.eint_input;
    state(i, j, k, UEDEN) = rho * (prob_parm.eint_input
      + 0.5 * (u_perturb * u_perturb + prob_parm.v_input * prob_parm.v_input
      + prob_parm.w_input * prob_parm.w_input));
  } else {
    rho = prob_parm.rho_air;
    state(i, j, k, UEINT) = rho * prob_parm.eint_air;
    state(i, j, k, UEDEN) = rho * (prob_parm.eint_air
      + 0.5 * (u_perturb * u_perturb + prob_parm.v_input * prob_parm.v_input
      + prob_parm.w_input * prob_parm.w_input));
  }
  state(i, j, k, URHO) = rho; // was this missing before???
  state(i, j, k, UTEMP) = prob_parm.T_input;
  state(i, j, k, UMX) = rho * u_perturb;
  state(i, j, k, UMY) = rho * prob_parm.v_input;
  state(i, j, k, UMZ) = rho * prob_parm.w_input;
  if (std::abs(y) < prob_parm.fuel_y_limit && x < prob_parm.fuel_x_limit) {
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS + n) = prob_parm.rho_input * prob_parm.massfrac[n];
    }
  } else {
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, UFS+n) = prob_parm.rho_air * prob_parm.massfrac_air[n];
    }
  }
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
  if (1) {
    amrex::Real rho = 0.0, u = 0.0, v = 0.0, w = 0.0, eint = 0.0, T = 0.0;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

    auto eos = pele::physics::PhysicsType::eos();

    if (sgn == 1) {
      const amrex::Real p = prob_parm.p_input;
      T = prob_parm.T_input;
      u = prob_parm.u_input; //cm/s
      v = 0.0;    //cm/s
      w = 0.0;    //cm/s

      //limit propane away from walls
      if (std::abs(x[1]) < prob_parm.fuel_y_limit) {
        rho = prob_parm.rho_input;
        eint = prob_parm.eint_input;
        for (int n=0; n < NUM_SPECIES; n++) {
          s_ext[UFS + n] = prob_parm.rho_input * prob_parm.massfrac[n];
	}
      } else {
        rho = prob_parm.rho_air;
        eint = prob_parm.eint_air;
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[UFS + n] = prob_parm.rho_air * prob_parm.massfrac_air[n];
	}
      }


    } else if (sgn == -1) {

      // Following Blazek p 279, eq. 8.23
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_int = s_int[URHO];
      const amrex::Real u_int = s_int[UMX] / rho_int;
      const amrex::Real v_int = s_int[UMY] / rho_int;
      const amrex::Real w_int = s_int[UMZ] / rho_int;
      const amrex::Real T_int = s_int[UTEMP];
      amrex::Real p_int = 0.0, cs_int = 0.0;

      // Get massfracs from interior points
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_int[n] = s_int[UFS + n] / rho_int;
      }

      eos.RTY2P(rho_int, T_int, massfrac_int.begin(), p_int);
      eos.RTY2Cs(rho_int, T_int, massfrac_int.begin(), cs_int);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p_input;
      const amrex::Real rhob =
        s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
      const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
      const amrex::Real vb = v_int;
      const amrex::Real wb = w_int;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_int) / (xb - xd) * (x[0] - xd) + rho_int;
      const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;
      
      //for (int n = 0; n < NUM_SPECIES; n++)
      //  massfrac[n] = s_inter[UFS + n] / rho_inter;  
 
      eos.RYP2E(rho, massfrac_int.begin(), p, eint);
      eos.EY2T(eint, massfrac_int.begin(), T);

      u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
      v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
      w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;

      // For Outflow
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[UFS + n] = s_int[UFS + n];
      }
    }

    s_ext[URHO] = rho;
    s_ext[UMX] = rho * u;
    s_ext[UMY] = rho * v;
    s_ext[UMZ] = rho * w;
    s_ext[UEINT] = rho * eint;
    s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
    s_ext[UTEMP] = T;
    //for (int n = 0; n < NUM_SPECIES; n++)
    //  s_ext[UFS + n] = rho * prob_parm.massfrac[n];
    //  
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif

