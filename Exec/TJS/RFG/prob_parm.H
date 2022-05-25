#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

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
  amrex::Real turb_length_scale;
  amrex::Real turb_velocity;
  static const int turb_num_modes = 50;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> xi;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> zeta;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> d;
  amrex::GpuArray<amrex::Real, turb_num_modes> omega;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> p;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> q;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac;
}; // namespace ProbParm

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
