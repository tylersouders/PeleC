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

    // Set useful variables
    auto eos = pele::physics::PhysicsType::eos();

    // Query user-changable values
    amrex::ParmParse pp("prob");
    pp.query("p_atm", PeleC::h_prob_parm_device->p_atm);
    pp.query("T_atm", PeleC::h_prob_parm_device->T_atm);
    pp.query("p0", PeleC::h_prob_parm_device->p0);
    pp.query("T0", PeleC::h_prob_parm_device->T0);
    pp.query("t_shutoff", PeleC::h_prob_parm_device->t_shutoff);

    // Convert the tank pressure (in psi) to Ba
    PeleC::h_prob_parm_device->p0 *= 68947.6;

    // Initialize species fractions for tank and atm
    for (int n = 0; n < NUM_SPECIES; n++) {
      PeleC::h_prob_parm_device->molefrac_atm[n] = 0.0;
      PeleC::h_prob_parm_device->molefrac_in[n] = 0.0;
    }

    // for initializing the domain with "engineering air"
    PeleC::h_prob_parm_device->molefrac_atm[O2_ID] = 0.21;
    PeleC::h_prob_parm_device->molefrac_atm[N2_ID] = 0.79;

    // Set inlet to pure CO2
    PeleC::h_prob_parm_device->molefrac_in[CO2_ID] = 1.0;

    eos.X2Y(PeleC::h_prob_parm_device->molefrac_atm.begin(), 
      PeleC::h_prob_parm_device->massfrac_atm.begin());
    eos.X2Y(PeleC::h_prob_parm_device->molefrac_in.begin(), 
      PeleC::h_prob_parm_device->massfrac_in.begin());

    // Initialize density and energy from mass fractions, T and P.
    amrex::Real rho0, eint0;
    eos.PYT2RE(
      PeleC::h_prob_parm_device->p_atm, PeleC::h_prob_parm_device->massfrac_atm.begin(), PeleC::h_prob_parm_device->T_atm,
      PeleC::h_prob_parm_device->rho_atm, PeleC::h_prob_parm_device->eint_atm);
    eos.PYT2RE(
      PeleC::h_prob_parm_device->p0, PeleC::h_prob_parm_device->massfrac_in.begin(), PeleC::h_prob_parm_device->T0,
      rho0, eint0);

    // Use EB limits to compute throat area
    amrex::ParmParse pp2("tailored_bb");
    amrex::Vector<amrex::Real> throatvec, boxhi;
    pp2.getarr("tri_0_point_2", throatvec);
    pp2.getarr("top_box_lo", boxhi);
    amrex::Real rt, r0, At, A0;
    rt = throatvec[1];
    r0 = boxhi[1];
    At = constants::PI() * std::pow(rt, 2.0);
    A0 = constants::PI() * std::pow(r0, 2.0);

    // Get spec heat ratio (using CO2 values but close enough)
    amrex::Real cp, cv, gamma;
    eos.TY2Cp(PeleC::h_prob_parm_device->T0, 
      PeleC::h_prob_parm_device->massfrac_in.begin(), cp);
    eos.RTY2Cv(rho0,
      PeleC::h_prob_parm_device->T0,
      PeleC::h_prob_parm_device->massfrac_in.begin(), cv);
    gamma = cp/cv;

    // Compute mass flow rate from At (everything in cgs units, should be OKAY)
    amrex::Real mdot = At * std::pow(
      gamma * rho0 * PeleC::h_prob_parm_device->p0 *
      std::pow((2.0 / (gamma + 1.0)), (gamma+1.0)/(gamma-1.0)),
      0.5);

    // Employ a bisection method search to compute Mach
    amrex::Real Mlo, M, Mhi, fa, fb, fc;
    amrex::Real err = 9999;
    amrex::Real errlim = 1.0e-8;
    int it = 0;

    Mlo = 0.01;
    Mhi = 0.999;
    while ((err > errlim) && (it < 100)) {

      // Evaluate guesses
      M = 0.5 * (Mlo + Mhi);

      fa = 1 / (std::pow(Mlo, 2.0)) * 
      std::pow(2.0/(gamma+1.0)*(1+(gamma-1.0)/2.0*std::pow(Mlo, 2.0)), (gamma+1.0)/(gamma-1.0)) -
      std::pow(A0/At, 2.0);
      fb = 1 / (std::pow(M, 2.0)) * 
      std::pow(2.0/(gamma+1.0)*(1+(gamma-1.0)/2.0*std::pow(M, 2.0)), (gamma+1.0)/(gamma-1.0)) -
      std::pow(A0/At, 2.0);
      fc = 1 / (std::pow(Mhi, 2.0)) * 
      std::pow(2.0/(gamma+1.0)*(1+(gamma-1.0)/2.0*std::pow(Mhi, 2.0)), (gamma+1.0)/(gamma-1.0)) -
      std::pow(A0/At, 2.0);

      if ((fb == 0.0) || ((Mhi - Mlo)/2.0 <- errlim)) {
        break;
      } else if ((fb/fa)>0.0) {
        Mlo = M;
      } else {
        Mhi = M;
      }
      it++;
    }

    // Save this Mach number
    prob_parm.Mach_inlet = M;

    // Use Total-to-Static to get actual inlet conditions!
    PeleC::h_prob_parm_device->p_in = PeleC::h_prob_parm_device->p0 / 
      std::pow(1.0 + (gamma-1.0)/2.0*std::pow(M, 2.0), gamma/(gamma-1.0));
    PeleC::h_prob_parm_device->T_in = PeleC::h_prob_parm_device->T0 /
      (1.0 + (gamma-1.0)/2.0*std::pow(M, 2.0));
    eos.PYT2RE(PeleC::h_prob_parm_device->p_in, PeleC::h_prob_parm_device->massfrac_in.begin(),
      PeleC::h_prob_parm_device->T_in, PeleC::h_prob_parm_device->rho_in, PeleC::h_prob_parm_device->eint_in);

    // Use mass flow rate to compute inlet velocity
    PeleC::h_prob_parm_device->u_in = mdot / (A0 * PeleC::h_prob_parm_device->rho_in);

    // Debug
    amrex::Print() << "\n\n";
    amrex::Print() << "----> ECHO CHOKING CONDITIONS <----\n";
    amrex::Print() << "Total Density = " << rho0 << " [g/cm^3]\n";
    amrex::Print() << "Total Pressure = " << PeleC::h_prob_parm_device->p0 << " [g/cm/s^2]\n";
    amrex::Print() << "Total Temperature = " << PeleC::h_prob_parm_device->T0 << " [K]\n";
    amrex::Print() << "Area Ratio = " << (A0/At) << "\n";
    amrex::Print() << "Mass Flow Rate = " << mdot << " [g/s]\n";
    amrex::Print() << "cp/cv = " << gamma << "\n\n";

    amrex::Print() << "Inlet Mach Number = " << M << "\n";
    amrex::Print() << "Inlet Pressure Ratio = " << (PeleC::h_prob_parm_device->p_in/PeleC::h_prob_parm_device->p0) << "\n";
    amrex::Print() << "Inlet Temperature Ratio = " << (PeleC::h_prob_parm_device->T_in/PeleC::h_prob_parm_device->T0) << "\n";
    amrex::Print() << "Inlet Density Ratio = " << (PeleC::h_prob_parm_device->rho_in/rho0) << "\n";
    amrex::Print() << "Inlet Velocity = " << (PeleC::h_prob_parm_device->u_in) << " [cm/s]\n";
    amrex::Print() << "\n\n";

  }

}

void
PeleC::problem_post_timestep()
{

  // if (amrex::ParallelDescriptor::IOProcessor()) {
  //   // Debug prints
  //   amrex::Print() << "prob_parm.p_in = " << PeleC::h_prob_parm_device->p_in << " [g/cm/s^2]\n";
  //   amrex::Print() << "prob_parm.T_in = " << PeleC::h_prob_parm_device->T_in << " [K]\n";
  //   amrex::Print() << "prob_parm.rho_in = " << PeleC::h_prob_parm_device->rho_in << " [g/cm^3]\n";
  //   amrex::Print() << "prob_parm.eint_in = " << PeleC::h_prob_parm_device->eint_in << " [units]\n";
  //   amrex::Print() << "prob_parm.u_in = " << PeleC::h_prob_parm_device->u_in << " [cm/s]\n";
  // }
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

// void
// CDNozzle::build(
//   const amrex::Geometry& geom, const int max_coarsening_level) {

//   const int npts_in_tri = 3;
//   const int max_tri = 3;

//   int num_tri;
//   amrex::ParmParse pp("tailored_bb");
//   // hard-code AMREX_SPACEDIM=3; this is for a triangle after all.
//   amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM>> alltri(
//     npts_in_tri * max_tri);

//   // initalize all triangles with some dummy values
//   // that fall outside of the domain
//   const amrex::Real* problo;
//   const amrex::Real* probhi;
//   amrex::Real maxlen;

//   problo = geom.ProbLo();
//   probhi = geom.ProbHi();

//   amrex::Vector<amrex::Real> topboxlo;
//   amrex::Vector<amrex::Real> topboxhi;
//   amrex::Vector<amrex::Real> botboxlo;
//   amrex::Vector<amrex::Real> botboxhi;
//   amrex::Vector<amrex::Real> bluffbodyhi;
//   amrex::Vector<amrex::Real> bluffbodylo;    

//   pp.getarr("top_box_lo", topboxlo);
//   pp.getarr("top_box_hi", topboxhi);
//   pp.getarr("bot_box_lo", botboxlo);
//   pp.getarr("bot_box_hi", botboxhi);

//   // Generate Top Box
//   amrex::EB2::BoxIF topbox(
//     {AMREX_D_DECL(topboxlo[0], topboxlo[1], problo[2]-1.0)},
//     {AMREX_D_DECL(topboxhi[0], topboxhi[1], probhi[2]+1.0)}, false);

//   // Generate Bottom Box
//   amrex::EB2::BoxIF bottombox(
//     {AMREX_D_DECL(botboxlo[0], botboxlo[1], problo[2]-1.0)},
//     {AMREX_D_DECL(botboxhi[0], botboxhi[1], probhi[2]+1.0)}, false);

//   maxlen = std::max(
//     std::max(geom.ProbLength(0), geom.ProbLength(1)), geom.ProbLength(2));

//   // // setting all triangles to be waaay outside the domain initially
//   // for (int itri = 0; itri < max_tri; itri++) {
//   //   alltri[npts_in_tri * itri + 0][0] = problo[0] + 100.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 0][1] = problo[1] + 100.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 0][2] = 0.0;

//   //   alltri[npts_in_tri * itri + 1][0] = probhi[0] + 101.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 1][1] = problo[1] + 100.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 1][2] = 0.0;

//   //   alltri[npts_in_tri * itri + 2][0] = probhi[0] + 101.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 2][1] = problo[1] + 101.0 * maxlen;
//   //   alltri[npts_in_tri * itri + 2][2] = 0.0;
//   // }

//   // tailored_bb.tri_2_point_0 = -0.34524  0.0    0.0
//   // tailored_bb.tri_2_point_1 =  0.001   -0.1999 0.0
//   // tailored_bb.tri_2_point_2 =  0.001    0.1999 0.0

//   // HARD CODING BAD BUT SOMETIMES NECESSARY
//   for (int itri = 0; itri < max_tri; itri++) {
//     alltri[npts_in_tri * itri + 0][0] = -0.34;
//     alltri[npts_in_tri * itri + 0][1] = 0.0;
//     alltri[npts_in_tri * itri + 0][2] = 0.0;

//     alltri[npts_in_tri * itri + 1][0] = 0.0;
//     alltri[npts_in_tri * itri + 1][1] = -0.199;
//     alltri[npts_in_tri * itri + 1][2] = 0.0;

//     alltri[npts_in_tri * itri + 2][0] = 0.0;
//     alltri[npts_in_tri * itri + 2][1] = 0.199;
//     alltri[npts_in_tri * itri + 2][2] = 0.0;
//   }

//   // get user defined number of triangles
//   pp.get("num_tri", num_tri);

//   for (int itri = 0; itri < num_tri; itri++) {
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> point{
//       AMREX_D_DECL(0.0, 0.0, 0.0)};

//     for (int ipt = 0; ipt < npts_in_tri; ipt++) {
//       std::string pointstr =
//         "tri_" + convertIntGG(itri) + "_point_" + convertIntGG(ipt);
//       amrex::Vector<amrex::Real> vecpt;
//       pp.getarr(pointstr.c_str(), vecpt, 0, AMREX_SPACEDIM);
//       for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
//         point[dir] = vecpt[dir];
//       }
//       alltri[npts_in_tri * itri + ipt] = point;
//     }
//   }

//   // intersection of the 3 planes in a triangle for all triangles
//   amrex::Vector<std::unique_ptr<amrex::EB2::IntersectionIF<
//     amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>>
//     impfunc_triangles(max_tri);

//   for (int itri = 0; itri < max_tri; itri++) {
//     // make sure points are in anti clockwise direction to set the inside of
//     // the triangle as solid phase correctly
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> norm0;
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> norm1;
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> norm2;

//     amrex::Array<amrex::Real, AMREX_SPACEDIM> point0;
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> point1;
//     amrex::Array<amrex::Real, AMREX_SPACEDIM> point2;

//     point0 = alltri[npts_in_tri * itri + 0];
//     point1 = alltri[npts_in_tri * itri + 1];
//     point2 = alltri[npts_in_tri * itri + 2];

//     norm0[0] = -(point1[1] - point0[1]);
//     norm0[1] = (point1[0] - point0[0]);
//     norm0[2] = 0.0;

//     norm1[0] = -(point2[1] - point1[1]);
//     norm1[1] = (point2[0] - point1[0]);
//     norm1[2] = 0.0;

//     norm2[0] = -(point0[1] - point2[1]);
//     norm2[1] = (point0[0] - point2[0]);
//     norm2[2] = 0.0;

//     // normalize so that magnitude is 1
//     amrex::Real norm = sqrt(norm0[0] * norm0[0] + norm0[1] * norm0[1]);
//     norm0[0] = norm0[0] / norm;
//     norm0[1] = norm0[1] / norm;

//     // normalize so that magnitude is 1
//     norm = sqrt(norm1[0] * norm1[0] + norm1[1] * norm1[1]);
//     norm1[0] = norm1[0] / norm;
//     norm1[1] = norm1[1] / norm;

//     // normalize so that magnitude is 1
//     norm = sqrt(norm2[0] * norm2[0] + norm2[1] * norm2[1]);
//     norm2[0] = norm2[0] / norm;
//     norm2[1] = norm2[1] / norm;

//     amrex::EB2::PlaneIF plane0(point0, norm0);
//     amrex::EB2::PlaneIF plane1(point1, norm1);
//     amrex::EB2::PlaneIF plane2(point2, norm2);

//     impfunc_triangles[itri] = std::make_unique<amrex::EB2::IntersectionIF<
//       amrex::EB2::PlaneIF, amrex::EB2::PlaneIF, amrex::EB2::PlaneIF>>(

//       plane0, plane1, plane2);
//   }

//   auto alltri_IF = amrex::EB2::makeUnion(
//     *impfunc_triangles[0], *impfunc_triangles[1], *impfunc_triangles[2]);

//   // Combine Geometry
//   auto comb_shape = amrex::EB2::makeUnion(alltri_IF, topbox, bottombox);

//   // auto comb_shape_extrude_IF = amrex::EB2::extrude(comb_shape, 2); // along z

//   auto gshop = amrex::EB2::makeShop(comb_shape);
//   amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);

// }
