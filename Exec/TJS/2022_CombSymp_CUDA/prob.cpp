#include "prob.H"

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
  // Set seed for random noise in initialization
  srand((unsigned) time(NULL));
  
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_input", PeleC::h_prob_parm_device->p_input);
  pp.query("T_input", PeleC::h_prob_parm_device->T_input);
  pp.query("u_input", PeleC::h_prob_parm_device->u_input);
  pp.query("v_input", PeleC::h_prob_parm_device->v_input);
  pp.query("phi_input", PeleC::h_prob_parm_device->phi_input);
  pp.query("init_pert_mag", PeleC::h_prob_parm_device->init_pert_mag);
  pp.query("fuel_y_limit", PeleC::h_prob_parm_device->fuel_y_limit);
  pp.query("fuel_x_limit", PeleC::h_prob_parm_device->fuel_x_limit);
  

  //Compute random number for pert
  PeleC::h_prob_parm_device->randnum = amrex::Random();    
  const amrex::Real H_n = 8.0;
  const amrex::Real C_m = 3.0;
  const amrex::Real a = C_m + H_n/4.0;

  // Compute mole fractions
  PeleC::h_prob_parm_device->molefrac[O2_ID] = 1.0 / (1.0 + PeleC::h_prob_parm_device->phi_input/a + 3.76);
  PeleC::h_prob_parm_device->molefrac[C3H8_ID] = PeleC::h_prob_parm_device->phi_input / a * 
    PeleC::h_prob_parm_device->molefrac[O2_ID];
  PeleC::h_prob_parm_device->molefrac[N2_ID] = 1.0 - 
    PeleC::h_prob_parm_device->molefrac[O2_ID] - PeleC::h_prob_parm_device->molefrac[C3H8_ID];

  // for "engineering air"
  PeleC::h_prob_parm_device->molefrac_air[O2_ID] = 0.21;
  PeleC::h_prob_parm_device->molefrac_air[N2_ID] = 0.79;
  
  // Convert to mass fraction
  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(PeleC::h_prob_parm_device->molefrac.begin(), PeleC::h_prob_parm_device->massfrac.begin());
  eos.X2Y(PeleC::h_prob_parm_device->molefrac_air.begin(), PeleC::h_prob_parm_device->massfrac_air.begin());

  //eos.RYP2E(
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_input, PeleC::h_prob_parm_device->massfrac.begin(),
    PeleC::h_prob_parm_device->T_input, PeleC::h_prob_parm_device->rho_input,PeleC::h_prob_parm_device->eint_input);
  eos.PYT2RE(
    PeleC::h_prob_parm_device->p_input, PeleC::h_prob_parm_device->massfrac_air.begin(),
    PeleC::h_prob_parm_device->T_input, PeleC::h_prob_parm_device->rho_air,PeleC::h_prob_parm_device->eint_air);
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
