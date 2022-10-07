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

      // Initialize all quantities with FOExtrap values
      amrex::Real rho = s_int[URHO];
      amrex::Real u = s_int[UMX] / rho;
      amrex::Real v = s_int[UMY] / rho;
      amrex::Real w = s_int[UMZ] / rho;
      amrex::Real eint = s_int[UEINT] / rho;
      amrex::Real T = s_int[UTEMP];
      amrex::Real massfrac[NUM_SPECIES];
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = s_int[UFS + n] / rho;
      }

      // Following Blazek p 279, eq. 8.23 - reimplemented from EB-C10
      // NOTE: it is assumed that mass fractions are constant across the
      // boundary Interior state (point d)
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real xd = prob_hi[0] - 0.5 * dx[0];
      const amrex::Real rho_inter = rho;
      const amrex::Real u_inter = u;
      const amrex::Real v_inter = v;
      const amrex::Real w_inter = w;
      const amrex::Real T_inter = T;
      amrex::Real p_inter = 0.0, cs_inter = 0.0;
      auto eos = pele::physics::PhysicsType::eos();
      eos.RTY2P(rho_inter, T_inter, massfrac, p_inter);
      eos.RTY2Cs(rho_inter, T_inter, massfrac, cs_inter);

      // Boundary state (point b)
      const amrex::Real xb = prob_hi[0];
      const amrex::Real pb = prob_parm.p_atm;
      const amrex::Real rhob =
        s_int[URHO] + (pb - p_inter) / (cs_inter * cs_inter);
      const amrex::Real ub = u_inter + (p_inter - pb) / (rho_inter * cs_inter);
      const amrex::Real vb = v_inter;
      const amrex::Real wb = w_inter;

      // Ghost state (point a). Linear extrapolation from d and b
      rho = (rhob - rho_inter) / (xb - xd) * (x[0] - xd) + rho_inter;
      const amrex::Real p = (pb - p_inter) / (xb - xd) * (x[0] - xd) + p_inter;

      eos.RYP2E(rho, massfrac, p, eint);
      eos.REY2T(rho, eint, massfrac, T);

      u = (ub - u_inter) / (xb - xd) * (x[0] - xd) + u_inter;
      v = (vb - v_inter) / (xb - xd) * (x[0] - xd) + v_inter;
      w = (wb - w_inter) / (xb - xd) * (x[0] - xd) + w_inter;

      // Write to s_ext
      s_ext[URHO] = rho;
      s_ext[UMX] = rho * u;
      s_ext[UMY] = rho * v;
      s_ext[UMZ] = rho * w;
      s_ext[UEINT] = rho * eint;
      s_ext[UEDEN] = rho * (eint + 0.5 * (u * u + v * v + w * w));
      s_ext[UTEMP] = T;
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[UFS + n] = rho * massfrac[n];
      }

    }

  }

}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
