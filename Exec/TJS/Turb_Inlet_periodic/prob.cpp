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

  // Set seed for noise. This makes it deterministic and supports restarts.
  // Change to AMReX random functions. Initialize turbulence inflow.
    const amrex::ULong deterministicSeed = 100;
    amrex::ResetRandomSeed(deterministicSeed);

  // We need to determine cell size via ParmParse here since geom data is not passed.
    amrex::Vector<int> n_cells;
    amrex::ParmParse ppamr("amr");
    ppamr.getarr("n_cell", n_cells);

  // consider using the finest level of AMR if we are using it.
  // For now we will just use the base level.
    int max_level = 0;
    ppamr.query("max_level", max_level);

  // Parse params for turbulence inlet. User can set u_in and I% at runtime
    amrex::ParmParse pp("prob");
    pp.query("turb_intensity", PeleC::h_prob_parm_device->turb_intensity);
    pp.query("u_input", PeleC::h_prob_parm_device->u_in);
    pp.query("init_turb_fill", PeleC::h_prob_parm_device->init_turb_fill);
    pp.query("do_turb_inlet", PeleC::h_prob_parm_device->do_turb_inlet);

  // Compute turbulence length scale based on inlet domain size. Assume x-dimension is much larger than y or z.
  // Currently computing as 1/4th of the domain size, can change if necessary
    PeleC::h_prob_parm_device->turb_length_scale = (probhi[2] - problo[2]);

  // Compute base (no AMR) grid density
    const amrex::Real dx_base = (probhi[1] - problo[1]) / n_cells[1];

  // Compute wavenumber limits (k_lo = largest length scale, k_hi = smallest)
    const amrex::Real k_lo = 2.0 * constants::PI() / PeleC::h_prob_parm_device->turb_length_scale;

  // Adjust for periodic z-dir
    amrex::Real Nzmax, intpart, fracpart;

  // Compute max z wavenumber
    const amrex::Real kz_max = 1.0 / (2.0 * dx_base);
    Nzmax = kz_max / (2.0 * constants::PI() / PeleC::h_prob_parm_device->turb_length_scale);
    fracpart = std::modf(Nzmax, &intpart);
    PeleC::h_prob_parm_device->turb_num_modes = intpart - 1; // Offset to account for 1 index in paper

  // Compute dk based on a discretization of M different discrete wavenumbers
    const amrex::Real dk = 2.0 * constants::PI() / PeleC::h_prob_parm_device->turb_length_scale;

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Checkpoint #1 (before loop) " << '\n';
    }

  // Initialize a few relevant variables.
    amrex::Vector<amrex::Real> zeta_cross_k(3);
    amrex::Vector<amrex::Real> xi_cross_k(3);
    amrex::Vector<amrex::Real> temp_k(3);
    amrex::Real k_mag, theta, phi, f, E, cross_mag, root_arg, 
    std_omega, randa, omega_upper_bound;

  // Loop to initialize and save arrays for zeta and xi to be used in next loop
    // Define xi and zeta
    amrex::Vector<amrex::Vector<amrex::Real>> xi, zeta;
    xi.resize(3);
    zeta.resize(3);
    for (int i = 0; i < 3; i++) {
      xi[i].resize(PeleC::h_prob_parm_device->sampling_number);
      zeta[i].resize(PeleC::h_prob_parm_device->sampling_number);
      for (int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {
        xi[i][n] = amrex::RandomNormal(0.0, 1.0);
        zeta[i][n] = amrex::RandomNormal(0.0, 1.0);
      }
    }

  // Loop to initialize and save arrays for p, q, k, omega
      for (int m = 0; m < PeleC::h_prob_parm_device->turb_num_modes; m++) {

    // wave number magnitude
        k_mag = k_lo + (m * dk);

    // Generate E for this wavenumber
        E = 1.5 * (4.0 * std::pow(PeleC::h_prob_parm_device->turb_intensity * PeleC::h_prob_parm_device->u_in, 2.0)
         * (PeleC::h_prob_parm_device->turb_length_scale / PeleC::h_prob_parm_device->u_in)
         / std::pow(1.0 + 70.8 * std::pow(k_mag * PeleC::h_prob_parm_device->turb_length_scale, 2.0), 5.0/6.0));

        omega_upper_bound = 2.0 * constants::PI() * k_mag * 
        PeleC::h_prob_parm_device->u_in;

        for (int n = 0; n < PeleC::h_prob_parm_device->sampling_number; n++) {

      // Generate random angles for isotropic sphere
          theta = amrex::Random() * 2.0 * constants::PI();
          phi   = amrex::Random() * constants::PI();

      // Convert to cartesian coords
          temp_k[0] = k_mag * cos(theta) * sin(phi);
          temp_k[1] = k_mag * sin(theta) * sin(phi);
          temp_k[2] = k_mag * cos(phi);

      // Generate random values for getting p and q
          // PeleC::prob_parm_host->h_omega[n][m] = amrex::RandomNormal(0.0, omega_upper_bound);
          // Access as 'omega[m*N + n]'
          PeleC::prob_parm_host->h_omega.push_back(amrex::RandomNormal(0.0, omega_upper_bound));
          randa = amrex::Random();

      // Calculate p, q

          xi_cross_k[0] = (xi[1][n] * temp_k[2] - temp_k[1] * xi[2][n]);
          xi_cross_k[1] = (xi[2][n] * temp_k[0] - temp_k[2] * xi[0][n]);
          xi_cross_k[2] = (xi[0][n] * temp_k[1] - temp_k[0] * xi[1][n]);

          cross_mag = std::sqrt(std::pow(xi_cross_k[0],2.0) +
            std::pow(xi_cross_k[1],2.0) + std::pow(xi_cross_k[2],2.0));

          root_arg = randa * 4.0 * E / PeleC::h_prob_parm_device->sampling_number;

          for (int i = 0; i < 3; i++) {
            // PeleC::prob_parm_host->h_p[i][n][m] = (xi_cross_k[i] / cross_mag) * std::sqrt(root_arg);
            PeleC::h_prob_parm_device->h_p.push_back((xi_cross_k[i] / cross_mag) * std::sqrt(root_arg));
          }

          zeta_cross_k[0] = (zeta[1][n] * temp_k[2] - temp_k[1] * zeta[2][n]);
          zeta_cross_k[1] = (zeta[2][n] * temp_k[0] - temp_k[2] * zeta[0][n]);
          zeta_cross_k[2] = (zeta[0][n] * temp_k[1] - temp_k[0] * zeta[1][n]);

          cross_mag = std::sqrt(std::pow(zeta_cross_k[0],2.0) +
            std::pow(zeta_cross_k[1],2.0) + std::pow(zeta_cross_k[2],2.0));

          root_arg = (1.0 - randa) * 4.0 * E / PeleC::h_prob_parm_device->sampling_number;

          for (int i = 0; i < 3; i++) {
            // PeleC::prob_parm_host->h_q[i][n][m] = (zeta_cross_k[i] / cross_mag) * std::sqrt(root_arg);
            PeleC::prob_parm_host->h_q.push_back((zeta_cross_k[i] / cross_mag) * std::sqrt(root_arg));
          }

          // Normalize k to k_tilde using the lowest wave number (equation 21 Huang 2010)
          for (int i = 0; i < 3; i++) {
            // PeleC::prob_parm_host->h_k[i][n][m] /= k_lo;
            PeleC::prob_parm_host->h_k.push_back(temp_k[i] / k_lo);

         }

       }

     }

     if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Checkpoint #4 (after loop) " << '\n';
    }

    pp.query("p_input", PeleC::h_prob_parm_device->p_init);
    pp.query("T_input", PeleC::h_prob_parm_device->T_init);
    pp.query("phi_input", PeleC::h_prob_parm_device->phi_in);

    pp.query("p_input", PeleC::h_prob_parm_device->p_in);
    pp.query("T_input", PeleC::h_prob_parm_device->T_in);

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

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Checkpoint #5 (end of probinit) " << '\n';
    }

  // Get pointer to variables p, q, omega, and k and cast to device
    PeleC::prob_parm_host->hG_k.resize(
      PeleC::prob_parm_host->h_k.size());
    PeleC::prob_parm_host->hG_p.resize(
      PeleC::prob_parm_host->p_k.size());
    PeleC::prob_parm_host->hG_q.resize(
      PeleC::prob_parm_host->q_k.size());
    PeleC::prob_parm_host->hG_omega.resize(
      PeleC::prob_parm_host->omega_k.size());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_k.begin(),
      PeleC::prob_parm_host->h_k.end(),
      PeleC::prob_parm_host->hG_k.begin());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_p.begin(),
      PeleC::prob_parm_host->h_p.end(),
      PeleC::prob_parm_host->hG_p.begin());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_q.begin(),
      PeleC::prob_parm_host->h_q.end(),
      PeleC::prob_parm_host->hG_q.begin());
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_omega.begin(),
      PeleC::prob_parm_host->h_omega.end(),
      PeleC::prob_parm_host->hG_omega.begin());

    PeleC::h_prob_parm_device->k = PeleC::prob_parm_host->hG_k.data();
    PeleC::h_prob_parm_device->p = PeleC::prob_parm_host->hG_p.data();
    PeleC::h_prob_parm_device->q = PeleC::prob_parm_host->hG_q.data();
    PeleC::h_prob_parm_device->omega = PeleC::prob_parm_host->hG_omega.data();

  }
}


void
PeleC::problem_post_timestep()
{
 if (amrex::ParallelDescriptor::IOProcessor()) {
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " p[0]  = " << PeleC::h_prob_parm_device->p[0][5][10] << '\n';
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " q[0]  = " << PeleC::h_prob_parm_device->q[0][5][10] << '\n';
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " k[0]  = " << PeleC::h_prob_parm_device->k[0][5][10] << '\n';
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " k[1]  = " << PeleC::h_prob_parm_device->k[1][5][10] << '\n';
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " omega  = " << PeleC::h_prob_parm_device->omega[5][10] << '\n';
   // amrex::Print() << "TIME= " << PeleC::h_prob_parm_device->timechk << " divu_check = " << PeleC::h_prob_parm_device->divu_check << '\n';
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
