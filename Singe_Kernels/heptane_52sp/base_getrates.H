#include "mechanism.H"

AMREX_GPU_HOST_DEVICE
void base_getrates(const double pressure, const double temperature, const double 
  avmolwt, const double *mass_frac, const double *diffusion, const double dt, 
  double *wdot); 
