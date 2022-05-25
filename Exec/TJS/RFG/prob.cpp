#include "prob.H"
#include <AMReX_Random.H>

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{

  // Set seed for noise. This makes it deterministic and supports restarts.
  // Change to AMReX random functions. Initialize turbulence inflow.
  const amrex::ULong deterministicSeed = 100;
  amrex::ResetRandomSeed(deterministicSeed);

  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("turb_length_scale", PeleC::h_prob_parm_device->turb_length_scale);
  pp.query("turb_velocity", PeleC::h_prob_parm_device->turb_velocity);

  amrex::Real std_d;
  if (AMREX_SPACEDIM == 3) {
    // 3-dimensional continuity enforcement
    std_d = 0.5;
  } else {
    // 2-dimensional continuity enforcement
    std_d = 1.0 / std::sqrt(3);
  }

  // This needs to be 3-dimensional so that we can take cross products
  // even if AMREX_SPACEDIM = 2.
  for (int i = 0; i < 3; i++) {

    /* 
    AMReX's buit-in amrex::RandomNormal does not allow us to sample from the same distribution.
    However we are testing here. If the seed iterates (which it should do), this may work and be
    more thread-secure.
    */
    for (int n = 0; n < PeleC::h_prob_parm_device->turb_num_modes; n++) {
      PeleC::h_prob_parm_device->xi[i][n]   = amrex::RandomNormal(0.0, 1.0);
      PeleC::h_prob_parm_device->zeta[i][n] = amrex::RandomNormal(0.0, 1.0);
      PeleC::h_prob_parm_device->d[i][n]    = amrex::RandomNormal(0.0, std_d);
    }
  }

  // Following Batten et al. (2004)
  for(int n = 0; n < PeleC::h_prob_parm_device->turb_num_modes; n++) {
    PeleC::h_prob_parm_device->omega[n] = amrex::RandomNormal(1.0, 1.0);
  }

  // p and q are the cross products of populated distributions
  //amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> p;
  //amrex::GpuArray<amrex::Array<amrex::Real, turb_num_modes>, 3> q;

  for(int n = 0; n < PeleC::h_prob_parm_device->turb_num_modes; n++) {
    PeleC::h_prob_parm_device->p[0][n] = (PeleC::h_prob_parm_device->xi[1][n]
                                          * PeleC::h_prob_parm_device->d[2][n])
                                          - (PeleC::h_prob_parm_device->d[1][n]
                                          * PeleC::h_prob_parm_device->xi[2][n]);
    PeleC::h_prob_parm_device->p[1][n] = (PeleC::h_prob_parm_device->xi[2][n]
                                          * PeleC::h_prob_parm_device->d[0][n])
                                          - (PeleC::h_prob_parm_device->xi[0][n]
                                          * PeleC::h_prob_parm_device->d[2][n]);
    PeleC::h_prob_parm_device->p[2][n] = (PeleC::h_prob_parm_device->xi[0][n]
                                          * PeleC::h_prob_parm_device->d[1][n])
                                          - (PeleC::h_prob_parm_device->d[0][n]
                                          * PeleC::h_prob_parm_device->xi[1][n]);
    PeleC::h_prob_parm_device->q[0][n] = (PeleC::h_prob_parm_device->zeta[1][n]
                                          * PeleC::h_prob_parm_device->d[2][n])
                                          - (PeleC::h_prob_parm_device->d[1][n]
                                          * PeleC::h_prob_parm_device->zeta[2][n]);
    PeleC::h_prob_parm_device->q[1][n] = (PeleC::h_prob_parm_device->zeta[2][n]
                                          * PeleC::h_prob_parm_device->d[0][n])
                                          - (PeleC::h_prob_parm_device->zeta[0][n]
                                          * PeleC::h_prob_parm_device->d[2][n]);
    PeleC::h_prob_parm_device->q[2][n] = (PeleC::h_prob_parm_device->zeta[0][n]
                                          * PeleC::h_prob_parm_device->d[1][n])
                                          - (PeleC::h_prob_parm_device->d[0][n]
                                          * PeleC::h_prob_parm_device->zeta[1][n]);
  }

  
  pp.query("p_init", PeleC::h_prob_parm_device->p_init);
  pp.query("T_init", PeleC::h_prob_parm_device->T_init);
  pp.query("phi_init", PeleC::h_prob_parm_device->phi_in);

  pp.query("p_in", PeleC::h_prob_parm_device->p_in);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);
  pp.query("u_in", PeleC::h_prob_parm_device->u_in);

  for (int n = 0; n < NUM_SPECIES; n++)
    PeleC::h_prob_parm_device->molefrac[n] = 0.0;

  // for CH4-air
  const amrex::Real a = 5.0;
  
  PeleC::h_prob_parm_device->molefrac[O2_ID] = 1.0 / (1.0 + (PeleC::h_prob_parm_device->phi_in / a)  + (0.79 / 0.21));
  PeleC::h_prob_parm_device->molefrac[C3H8_ID] = PeleC::h_prob_parm_device->phi_in * PeleC::h_prob_parm_device->molefrac[O2_ID] / a;
  PeleC::h_prob_parm_device->molefrac[N2_ID] = 1.0 - PeleC::h_prob_parm_device->molefrac[C3H8_ID] - PeleC::h_prob_parm_device->molefrac[O2_ID];

  // for initializing the domain with "engineering air"
  // ProbParm::molefrac_init[O2_ID] = 0.21;
  // ProbParm::molefrac_init[N2_ID] = 0.79;

  // Convert mole fracs to mass fracs
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(PeleC::h_prob_parm_device->molefrac.begin(), PeleC::h_prob_parm_device->massfrac.begin());

  // Initialize density and energy from mass fractions, T and P.
  eos.PYT2RE(
	     PeleC::h_prob_parm_device->p_in, PeleC::h_prob_parm_device->massfrac.begin(), PeleC::h_prob_parm_device->T_in,
	     PeleC::h_prob_parm_device->rho_in, PeleC::h_prob_parm_device->e_in);
}
}


void
PeleC::problem_post_timestep()
{
  if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "TIME= " << time << " p[0]  = " << PeleC::h_prob_parm_device->p[0][50] << '\n';
      amrex::Print() << "TIME= " << time << " q[0]  = " << PeleC::h_prob_parm_device->q[0][50] << '\n';
      amrex::Print() << "TIME= " << time << " d[0]  = " << PeleC::h_prob_parm_device->d[0][50] << '\n';
      amrex::Print() << "TIME= " << time << " d[1]  = " << PeleC::h_prob_parm_device->d[1][50] << '\n';
      amrex::Print() << "TIME= " << time << " omega  = " << PeleC::h_prob_parm_device->omega[50] << '\n';
  }
}


void
PeleC::problem_post_init()
{
}


void
PeleC::problem_post_restart()
{
}
