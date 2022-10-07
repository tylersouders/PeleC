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

// EB things
#include "EB.H"
#include <AMReX_EB2_IF_Spline.H>
#include "Utilities.H"
#include "Geometry.H"
class CDNozzle : public pele::pelec::Geometry::Register<CDNozzle>
{
public:
  static std::string identifier() { return "cdnozzle"; }

  void
  build(const amrex::Geometry& geom, const int max_coarsening_level) override;
};

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
  // Initial Condition - Quiescent Flow at STP, All Air
  // Set the state
  state(i, j, k, URHO) = prob_parm.rho_atm;
  state(i, j, k, UMX) = 0.0;
  state(i, j, k, UMY) = 0.0;
  state(i, j, k, UMZ) = 0.0;
  state(i, j, k, UEINT) = prob_parm.rho_atm * prob_parm.eint_atm;
  state(i, j, k, UEDEN) = prob_parm.rho_atm * 
  (prob_parm.eint_atm + 0.5 * (prob_parm.u_atm * prob_parm.u_atm));
  state(i, j, k, UTEMP) = prob_parm.T_atm;
  for (int n = 0; n < NUM_SPECIES; n++)
    state(i, j, k, UFS + n) = prob_parm.rho_atm * prob_parm.massfrac_atm[n];
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
  ProbParmDevice const& prob_parm) {

  // Initialize eos
  auto eos = pele::physics::PhysicsType::eos();

  if (idir == 0) {

    // x-direction BC

    if (sgn == 1) {

      // Simple Dirichlet Inlet
      // Inlet BC - set as dirichlet conditions computed in prob.cpp
      if (time < prob_parm.t_shutoff) {
        s_ext[URHO] = prob_parm.rho_in;
        s_ext[UMX] = prob_parm.rho_in * prob_parm.u_in;
        s_ext[UMY] = 0.0;
        s_ext[UMZ] = 0.0;
        s_ext[UTEMP] = prob_parm.T_in;
        s_ext[UEINT] = prob_parm.rho_in * prob_parm.eint_in;
        s_ext[UEDEN] = prob_parm.rho_in * (prob_parm.eint_in + 0.5 * (prob_parm.u_in*prob_parm.u_in));
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[UFS + n] = prob_parm.rho_in * prob_parm.massfrac_in[n];
        }
      } else {
        // Treat as a wall
        s_ext[URHO] = s_int[URHO];
        s_ext[UMX] = -s_int[UMX];
        s_ext[UMY] = -s_int[UMY];
        s_ext[UMZ] = -s_int[UMZ];
        s_ext[UTEMP] = s_int[UTEMP];
        s_ext[UEINT] = s_int[UEINT];
        s_ext[UEDEN] = s_int[UEDEN];
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[UFS + n] = s_int[UFS + n];
        }
      }

    } else if (sgn == -1) {

      /*
      There are three outlet cases to consider here:
      (1) Interior Flow is supersonic (M >= 1.0), which requires
          an extrapolation condition. No upstream information is
          possible.
      (2) Interior Flow is subsonic (M < 1.0) with positive x-dir
          flow. Use the Blazek non-reflecting boundary condition.
      (3) Interior Flow is subsonic (M < 1.0) and has reversed s.t.
          u < 0.0. I am not 100% sure of the best way to handle
          this, but my best guess is to enforce a Neumann-type 
          boundary condition with p=p_atm at the outlet.
      */

      // Compute the Mach number at the outlet
      // Initialize massfrac array for interior
      amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_int = {0.0};

      // Get species composition at the inlet
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_int[n] = s_int[UFS+n] / s_int[URHO];
      }

      // Compute sound speed and interior Mach number
      amrex::Real cs, Mach;
      eos.RTY2Cs(s_int[URHO], s_int[UTEMP], massfrac_int.begin(),
        cs);
      Mach = std::pow(s_int[UMX]*s_int[UMX] + 
        s_int[UMY]*s_int[UMY] + s_int[UMZ]*s_int[UMZ], 0.5) /
        s_int[URHO] / cs;

      // Get the interior direction to check for flow reversal
      const amrex::Real uchk = s_int[UMX] / s_int[URHO];

      // if ((Mach < prob_parm.outlet_mach_lim) && (uchk > 0.0)) {

      //   // Initialize some variables
      //   amrex::Real rho = 0.0, eint = 0.0, T = 0.0;
      //   amrex::Real u = 0.0, v = 0.0, w = 0.0;

      //   // Apply the Blazek non-reflecting BC
      //   // Interior state (point d)
      //   const amrex::Real* prob_hi = geomdata.ProbHi();
      //   const amrex::Real* dx = geomdata.CellSize();
      //   const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      //   const amrex::Real rho_int = s_int[URHO];
      //   const amrex::Real u_int = s_int[UMX] / rho_int;
      //   const amrex::Real v_int = s_int[UMY] / rho_int;
      //   const amrex::Real w_int = s_int[UMZ] / rho_int;
      //   const amrex::Real T_int = s_int[UTEMP];
      //   amrex::Real p_int = 0.0, cs_int = 0.0;

      //   auto eos = pele::physics::PhysicsType::eos();
      //   eos.RTY2P(rho_int, T_int, massfrac_int.begin(), p_int);
      //   eos.RTY2Cs(rho_int, T_int, massfrac_int.begin(), cs_int);

      //   // Boundary state (point b)
      //   const amrex::Real xb = prob_hi[0];
      //   const amrex::Real pb = prob_parm.p_atm;
      //   const amrex::Real rhob = s_int[URHO] + (pb - p_int) / (cs_int * cs_int);
      //   const amrex::Real ub = u_int + (p_int - pb) / (rho_int * cs_int);
      //   const amrex::Real vb = v_int;
      //   const amrex::Real wb = w_int;

      //   // Ghost state (point a). Linear extrapolation from d and b
      //   rho = (rhob - rho_int) / (xb - xd) * (x[0] - xd) + rho_int;
      //   const amrex::Real p = (pb - p_int) / (xb - xd) * (x[0] - xd) + p_int;

      //   eos.RYP2E(rho, massfrac_int.begin(), p, eint);
      //   eos.EY2T(eint, massfrac_int.begin(), T);

      //   u = (ub - u_int) / (xb - xd) * (x[0] - xd) + u_int;
      //   v = (vb - v_int) / (xb - xd) * (x[0] - xd) + v_int;
      //   w = (wb - w_int) / (xb - xd) * (x[0] - xd) + w_int;

      //   // Set state information
      //   s_ext[URHO] = rho;
      //   s_ext[UMX] = rho * u;
      //   s_ext[UMY] = rho * v;
      //   s_ext[UMZ] = rho * w;
      //   s_ext[UEINT] = rho * eint;
      //   s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
      //   s_ext[UTEMP] = T;
      //   for (int n = 0; n < NUM_SPECIES; n++) {
      //    s_ext[UFS + n] = s_int[UFS+n];
      //   }

      // } else if ((Mach < prob_parm.outlet_mach_lim) && (uchk <= 0.0)) {

      //   // Get interior state information
      //   const amrex::Real u_int = s_int[UMX] / s_int[URHO];
      //   const amrex::Real v_int = s_int[UMY] / s_int[URHO];
      //   const amrex::Real w_int = s_int[UMZ] / s_int[URHO];

      //   // Set state information
      //   s_ext[URHO] = prob_parm.rho_atm;
      //   s_ext[UMX] = prob_parm.rho_atm * u_int;
      //   s_ext[UMY] = prob_parm.rho_atm * v_int;
      //   s_ext[UMZ] = prob_parm.rho_atm * w_int;
      //   s_ext[UEINT] = prob_parm.rho_atm * prob_parm.eint_atm;
      //   s_ext[UEDEN] = prob_parm.rho_atm * (prob_parm.eint_atm + 0.5 * (u_int * u_int 
      //     + v_int * v_int + w_int * w_int));
      //   s_ext[UTEMP] = prob_parm.T_atm;
      //   for (int n = 0; n < NUM_SPECIES; n++) {
      //     s_ext[UFS + n] = prob_parm.rho_atm * prob_parm.massfrac_atm[n];
      //   }

      if ((Mach < prob_parm.outlet_mach_lim) && (uchk > 0.0)) {

      } else {

        // Supersonic outlet - only option is to extrapolate
        for (int n = 0; n < NVAR; n++) {
          s_ext[n] = s_int[n];
        }

      }
    }
  }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
