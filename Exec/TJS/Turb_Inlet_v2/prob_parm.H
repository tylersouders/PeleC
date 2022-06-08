#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

#include "Constants.H"

int compute_turb_num_modes() {

  // Use parmparse to get array of size
  amrex::ParmParse ppgeo("geometry");
  amrex::ParmParse ppamr("amr");

  amrex::Real prob_lo, prob_hi;
  ppgeo.getarr("prob_lo", prob_lo);
  ppgeo.getarr("prob_hi", prob_hi);

  int n_cells[3];
  ppamr.getarr("n_cell", n_cells);

  const amrex::Real Lz = prob_hi[2] - prob_lo[2];
  const amrex::Real dx_base = (probhi[1] - problo[1]) / n_cells[1];

  // Determine the number of elements
  // Compute wavenumber limits (k_lo = largest length scale, k_hi = smallest)
  const amrex::Real k_lo = 2.0 * constants::PI() / Lz;

  amrex::Real Nzmax, intpart, fracpart;
  const amrex::Real kz_max = 1.0 / (2.0 * dx_base);

  Nzmax = kz_max / (2.0 * constants::PI() / PeleC::h_prob_parm_device->turb_length_scale);
  fracpart = std::modf(Nzmax, &intpart);
  return intpart - 1; // Offset to account for 1 index in paper

}

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
  static const int turb_num_modes = compute_turb_num_modes();
  static const int sampling_number = 100;
  amrex::GpuArray<amrex::Array<amrex::Array<amrex::Real, turb_num_modes>, sampling_number>, 3> k;
  amrex::GpuArray<amrex::Array<amrex::Array<amrex::Real, turb_num_modes>, sampling_number>, 3> p;
  amrex::GpuArray<amrex::Array<amrex::Array<amrex::Real, turb_num_modes>, sampling_number>, 3> q;
  amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, sampling_number> omega;
  amrex::GpuArray<amrex::Array<amrex::Real, sampling_number>, 3> xi;
  amrex::GpuArray<amrex::Array<amrex::Real, sampling_number>, 3> zeta;
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
  ProbParmHost() {}
};

#endif
