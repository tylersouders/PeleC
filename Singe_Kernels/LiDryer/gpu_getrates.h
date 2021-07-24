#include "mechanism.H"

AMREX_GPU_GLOBAL
void
gpu_getrates(const double *temperature_array, const double *pressure_array, 
  const double *avmolwt_array, const double *mass_frac_array, const double *diffusion,
  const double dt, const int spec_stride/*NX*NY*NZ in number of doubles*/, double *wdot_array); 
