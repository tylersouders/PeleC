#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Vector.H>

struct ProbParmDevice
{
  amrex::Real p_init = 1013250.0;
  amrex::Real T_init = 310.0;
  amrex::Real rho_init = 0.0;
  amrex::Real e_init = 0.0;
  amrex::Real p_in = 1013250.0;
  amrex::Real T_in = 310.0;
  amrex::Real rho_in = 0.0;
  amrex::Real e_in = 0.0;
  amrex::Real phi_in = 0.65;
  amrex::Real u_in = 1590.0;
  amrex::Real init_pert_mag = 0.02;
  amrex::Real turb_length_scale = 0.002;
  amrex::Real turb_intensity = 7;
  amrex::Real Ls_prefactor = 1.5;
  int turb_num_modes = 100;
  static const int sampling_number = 100;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> k;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> p;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> q;
  amrex::Vector<amrex::Vector<amrex::Real>> omega;
  amrex::Vector<amrex::Vector<amrex::Real>> xi;
  amrex::Vector<amrex::Vector<amrex::Real>> zeta;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac;
  amrex::Real xchk[3] = {0.0, 0.10, 0.10};
  amrex::Real divu_check = 0.0;
  amrex::Real timechk = 0.0;
  bool init_turb_fill = 0;
  bool do_turb_inlet = 0;
}; // namespace ProbParm

struct ProbParmHost
{

  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> k;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> p;
  amrex::Vector<amrex::Vector<amrex::Vector<amrex::Real>>> q;
  amrex::Vector<amrex::Vector<amrex::Real>> omega;
  amrex::Vector<amrex::Vector<amrex::Real>> xi;
  amrex::Vector<amrex::Vector<amrex::Real>> zeta;

  ProbParmHost() {}
};

#endif
