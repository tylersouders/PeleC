#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

// #include "mechanism.H"

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Forcing.H"
#include "Utilities.H"
#include "turb_inflow.H"

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
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, energy, T, pres, e;

  amrex::Real kernel_location[3];
  const amrex::Real kernel_diameter = 100.0e-04;
  const amrex::Real kernel_y_loc    = 1.92;
  const amrex::Real kernel_height   = 1.01e-01;
  double r;

  kernel_location[0] = 0.0;
  kernel_location[1] = kernel_y_loc-kernel_height/2.0;
  kernel_location[2] = 0.0;

    for (int n = 0; n < NVAR; n++)
      state(i, j, k, n) = prob_parm.fuel_state[n];    

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
  ProbParmDevice const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  amrex::Real velx, vely, velz;
  amrex::Real rho, T, pres;
  
  if(idir == 0){
    if(PeleC::h_prob_parm_device->turbulence){
      // int i = (x[0] + 0.3+dx[0]*4)/(dx[0]);
      // int k = (x[2] + 0.0125)/dx[2];
      // printf("i = %i x = %f  j = %i  k = %i \n",i,x[0],j,k);

      // int n_cell[3] = {0};
      // amrex::ParmParse parse1("amr");
      // parse1.query("n_cell", n_cell);
      // printf("nx = %i ny = %i nz = %i \n",n_cell[0],n_cell[1],n_cell[2] );

      // amrex::ParmParse parse2("geometry");
      // amrex::Real prob_lo_global[3] = {0};
      // parse2.query("prob_lo", prob_lo_global);
      // printf("n_cell = %i %i %i\n", n_cell[0],n_cell[1],n_cell[2]);
      // printf("prob_lo_global = %f %f %f\n", prob_lo_global[0],prob_lo_global[1],prob_lo_global[2]);
      // amrex::Abort();


      //I am starting in 2D so I have to match only the cell index for y 
      int j = (x[1] - 1.02)/dx[1];
      if(j > 255)j=255;
      if(j < 0)j=0;
      // j = 0;
      
      //  int k = 0;
      int k = (x[2] + 0.8)/dx[2];
      // printf("k = %i \n", k);
      if(k > 255)k=255;
      if(k < 0)k=0;


      amrex::Real u[3] = {0.0};
      amrex::Real molefrac[NUM_SPECIES] = {0.0};
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      amrex::Real t = time + PeleC::h_prob_parm_device->time_init_turb;
      // amrex::Print() << "Overall time = " << t << "time = " << time << "time from plt = " << PeleC::h_prob_parm_device->time_init_turb << std::endl;
      int increment = PeleC::h_prob_parm_device->nt; //number of cells to read in x every time read_turbulence_bc() is called

      
      if(PeleC::h_prob_parm_device->first_time){
        // printf("Reading turbulence for the first time \n");
        // amrex::Real t = PeleC::h_prob_parm_device->time_init_turb; //0.00314002107022254; //time;
        
        //find index for time direction based on current time
        int idx = (t * PeleC::h_prob_parm_device->vn_in / PeleC::h_prob_parm_device->dx_turb) - 2;
        amrex::Print() << "Reading BC data for the first time" <<std::endl;
        
        if(idx < 0)idx = 0;
        // read file for the first time based on current simulation time
        PeleC::h_prob_parm_device->x_index[0] = idx;
        PeleC::h_prob_parm_device->x_index[1] = idx + increment;

        // printf("target time = %e \n", t);
        // printf("Index range: %i and %i \n",PeleC::h_prob_parm_device->x_index[0],PeleC::h_prob_parm_device->x_index[1]);
        read_turbulent_bc(PeleC::h_prob_parm_device->x_index,PeleC::h_prob_parm_device->vn_in);
        PeleC::h_prob_parm_device->first_time = false;
      }

      if(t > PeleC::h_prob_parm_device->d_xarray_turb[increment-1] + PeleC::h_prob_parm_device->recycled_turb_file){

        amrex::Print() << "Reading new turbulence data" << std::endl;
        
        // if current time is greater than the time read from file, read it again with the new data
        PeleC::h_prob_parm_device->x_index[0] = PeleC::h_prob_parm_device->x_index[1]-2;
        PeleC::h_prob_parm_device->x_index[1] = PeleC::h_prob_parm_device->x_index[1] + increment - 2;
        //start reading from time zero if the end of HIT BC was reached
        if(t > 1024 * PeleC::h_prob_parm_device->dx_turb / PeleC::h_prob_parm_device->vn_in + PeleC::h_prob_parm_device->recycled_turb_file){
          PeleC::h_prob_parm_device->recycled_turb_file += 1024 * PeleC::h_prob_parm_device->dx_turb / PeleC::h_prob_parm_device->vn_in ;
          amrex::Print() << "Restarting HIT BC..." <<std::endl;
          PeleC::h_prob_parm_device->x_index[0] = 0;
          PeleC::h_prob_parm_device->x_index[1] = increment;        
        }
        read_turbulent_bc(PeleC::h_prob_parm_device->x_index,PeleC::h_prob_parm_device->vn_in);
      }


      // if(prob_parm.turbulence){

      amrex::Real uinterp[3] = {0.0};

      // Interpolation factors
      amrex::Real mod[3] = {0.0};
      int idx[3] = {0};
      int idxp1[3] = {0};
      amrex::Real slp[3] = {0.0};


      mod[0] = std::fmod(t, PeleC::h_prob_parm_device->Lxturb);
      locate(PeleC::h_prob_parm_device->d_xarray_turb, PeleC::h_prob_parm_device->ires_x, mod[0], idx[0]);
      idxp1[0] = (idx[0] + 1) % PeleC::h_prob_parm_device->ires_x;
      slp[0] = (mod[0] - PeleC::h_prob_parm_device->d_xarray_turb[idx[0]]) / PeleC::h_prob_parm_device->d_xdiff_turb[idx[0]];

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);

      int nx = PeleC::h_prob_parm_device->ires_x;
      int ny = PeleC::h_prob_parm_device->ires_y;
      int nz = PeleC::h_prob_parm_device->ires_z;

      // Interpolate data in time direction
      for (int iscal=0;iscal<3;iscal++){
        uinterp[iscal] = PeleC::h_prob_parm_device->d_data_ic_turb[idx[0]   + nx*j   + nx*ny*k   + iscal*nx*ny*nz]*f0 + 
                         PeleC::h_prob_parm_device->d_data_ic_turb[idxp1[0] + nx*j   + nx*ny*k   + iscal*nx*ny*nz]*f1;
      }


      velx = uinterp[0] + PeleC::h_prob_parm_device->vn_in;
      vely = uinterp[1];
      velz = uinterp[2];      
    }
    else{
      velx = 0.0 + PeleC::h_prob_parm_device->vn_in;
      vely = 0.0;
      velz = 0.0;
    }

    rho = prob_parm.fuel_state[URHO];

    s_ext[URHO]  = rho;
    s_ext[UMX]   = rho*velx;
    s_ext[UMY]   = rho*vely;
    s_ext[UMZ]   = rho*velz;
    s_ext[UEINT] = prob_parm.fuel_state[UEINT];
    s_ext[UEDEN] = prob_parm.fuel_state[UEINT] + rho * (0.5 * (velx * velx + vely * vely + velz * velz));
    s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];

    // printf("End of turbulent BC \n");
  }
  else if(idir == 1){
    amrex::Real velx = 0.0;
    amrex::Real vely = 0.0;
    amrex::Real velz = 0.0;
    s_ext[URHO]  = prob_parm.fuel_state[URHO];
    s_ext[UMX]   = prob_parm.fuel_state[URHO]*velx;
    s_ext[UMY]   = prob_parm.fuel_state[URHO]*vely;
    s_ext[UMZ]   = prob_parm.fuel_state[URHO]*velz;
    s_ext[UEINT] = prob_parm.fuel_state[UEINT];
    s_ext[UEDEN] = prob_parm.fuel_state[UEINT] + prob_parm.fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
    s_ext[UTEMP] = prob_parm.fuel_state[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++)
      s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
      // s_ext[UFS + n] = prob_parm.fuel_state[UFS + n];
  }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void problem_post_init(//){};  
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real rho, energy, T, pres, e;

  amrex::Real kernel_location[3];
  const amrex::Real kernel_diameter = 100.0e-04;
  const amrex::Real kernel_y_loc    = 1.92;
  const amrex::Real kernel_height   = 1.01e-01;
  double r;

  kernel_location[0] = 0.0;
  kernel_location[1] = kernel_y_loc-kernel_height/2.0;
  kernel_location[2] = 0.0;

  if(PeleC::h_prob_parm_device->init_kernel){
    
    r = sqrt(pow((kernel_location[0]-x),2)+pow((kernel_location[2]-z),2));
    if(y > PeleC::h_prob_parm_device->d_yarray[0] and y < PeleC::h_prob_parm_device->d_yarray[PeleC::h_prob_parm_device->ny-1] 
      and r < PeleC::h_prob_parm_device->d_xarray[PeleC::h_prob_parm_device->nx-1]){

      // Fill in the velocities and interpolated quantities.
      amrex::Real u[3] = {0.0};
      // amrex::Real uinterp[PeleC::h_prob_parm_device->nscal] = {0.0};
      amrex::Real uinterp[10] = {0.0};

      // Interpolation factors
      amrex::Real mod[3] = {0.0};
      int idx   = 0;
      int idxp1 = 0;
      int idy   = 0;
      int idyp1 = 0;
      amrex::Real slp[3] = {0.0};
      
      mod[0] = std::fmod(r, PeleC::h_prob_parm_device->Lxinput);
      locate(PeleC::h_prob_parm_device->d_xarray, PeleC::h_prob_parm_device->nx, mod[0], idx);
      idxp1 = (idx + 1) % PeleC::h_prob_parm_device->nx;
      //idxp1 = (idx + 1);
      slp[0] = (mod[0] - PeleC::h_prob_parm_device->d_xarray[idx]) / PeleC::h_prob_parm_device->d_xdiff[idx];

      mod[1] = std::fmod(y, PeleC::h_prob_parm_device->Lyinput);
      locate(PeleC::h_prob_parm_device->d_yarray, PeleC::h_prob_parm_device->ny, mod[1], idy);
      idyp1 = (idy + 1) % PeleC::h_prob_parm_device->ny;
      //idyp1 = (idy + 1);
      slp[1] = (mod[1] - PeleC::h_prob_parm_device->d_yarray[idy]) / PeleC::h_prob_parm_device->d_ydiff[idy];

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]);
      const amrex::Real f2 = (1 - slp[0]) * slp[1];
      const amrex::Real f3 = slp[0] * slp[1];

      // Interpolate data
      for (int iscal=0;iscal<PeleC::h_prob_parm_device->nscal;iscal++){
        uinterp[iscal] = PeleC::h_prob_parm_device->d_data_ic[idx  +idy  *PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f0 +
                         PeleC::h_prob_parm_device->d_data_ic[idxp1+idy  *PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f1 +
                         PeleC::h_prob_parm_device->d_data_ic[idx  +idyp1*PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f2 +
                         PeleC::h_prob_parm_device->d_data_ic[idxp1+idyp1*PeleC::h_prob_parm_device->nx+iscal*PeleC::h_prob_parm_device->nx*PeleC::h_prob_parm_device->ny]*f3;
      
      }

      // // temperature, pressure, internal energy and composition
      double p = uinterp[1]*10.0; //Pa to dyn/cm2
      // rho = uinterp[0]/1.0e+3;
      T = uinterp[2];
      double massfrac[NUM_SPECIES]={0.0};

      for (int n = 0; n < NUM_SPECIES; n++)
        massfrac[n] = uinterp[n+6];

      double sum = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++)
        sum = sum + massfrac[n];
      // printf("Sum of all species = % e \n", sum);
      massfrac[N2_ID] = massfrac[N2_ID] + 1.0 - sum; //making sure that sum of ys is 1


      auto eos = pele::physics::PhysicsType::eos();
      // eos.RTY2P(rho, T, massfrac, p);
      eos.PYT2RE(p, massfrac, T, rho, e);

      double beta  = atan(z/x);
      double alpha = atan(x/z);

      //Velocity in each direction
      if(x >= 0.0){
        u[0] = uinterp[5]*cos(beta)*100.0;
      }
      else{
        u[0] = -uinterp[5]*cos(beta)*100.0;
      }

      u[1] = uinterp[4]*100.0;

      if(z >= 0){
        u[2] = uinterp[5]*cos(alpha)*100.0;
      }
      else{
        u[2] = -uinterp[5]*cos(alpha)*100.0;
      }


      // amrex::Real rho_init = state(i, j, k, URHO);
      amrex::Real velx = u[0] + state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
      amrex::Real vely = u[1] + state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
      amrex::Real velz = u[2] + state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum

      for (int n = 0; n < NUM_SPECIES; n++)
        state(i, j, k, UFS + n) = rho * massfrac[n];

      state(i, j, k, URHO)  = rho;
      state(i, j, k, UTEMP) = T;
      state(i, j, k, UMX)   = rho * velx;
      state(i, j, k, UMY)   = rho * vely;
      state(i, j, k, UMZ)   = rho * velz;
      state(i, j, k, UEINT) = rho * e;
      state(i, j, k, UEDEN) = rho * (e + 0.5 * (velx * velx + vely * vely + velz * velz)); 

      amrex::Real spec_sum = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        spec_sum = spec_sum + state(i, j, k, UFS + n);
      }
      if (
        amrex::Math::abs(state(i, j, k, URHO) - spec_sum) >
        1.e-6 * state(i, j, k, URHO)) {
        printf("In problem_post_restart: Sum of (rho X)_i = %e vs rho = %e at (%i,%i,%i)\n",spec_sum,state(i,j,k,URHO),i,j,k);
      }

    }  
    else{//initialize premixed state outside the kernel region
      // amrex::Real rho  = state(i, j, k, URHO);
      amrex::Real velx = state(i, j, k, UMX); //state() at this point is storing velocity, not momentum
      amrex::Real vely = state(i, j, k, UMY); //state() at this point is storing velocity, not momentum
      amrex::Real velz = state(i, j, k, UMZ); //state() at this point is storing velocity, not momentum
      state(i, j, k,URHO)  = PeleC::h_prob_parm_device->fuel_state[URHO];
      state(i, j, k,UMX)   = PeleC::h_prob_parm_device->fuel_state[URHO]*velx;
      state(i, j, k,UMY)   = PeleC::h_prob_parm_device->fuel_state[URHO]*vely;
      state(i, j, k,UMZ)   = PeleC::h_prob_parm_device->fuel_state[URHO]*velz;
      state(i, j, k,UEINT) = PeleC::h_prob_parm_device->fuel_state[UEINT];
      state(i, j, k,UEDEN) = PeleC::h_prob_parm_device->fuel_state[UEINT] + PeleC::h_prob_parm_device->fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
      state(i, j, k,UTEMP) = PeleC::h_prob_parm_device->fuel_state[UTEMP];
      for (int n = 0; n < NUM_SPECIES; n++)
        state(i, j, k,UFS + n) = PeleC::h_prob_parm_device->fuel_state[UFS + n];


    }
  }
  else{ // init_kernel = false: initialize premixed state with velocity field coming from plt file
    // amrex::Real rho  = state(i, j, k, URHO)
    amrex::Real velx = state(i, j, k, UMX);
    amrex::Real vely = state(i, j, k, UMY);
    amrex::Real velz = state(i, j, k, UMZ);
    state(i, j, k,URHO)  = PeleC::h_prob_parm_device->fuel_state[URHO];
    state(i, j, k,UMX)   = PeleC::h_prob_parm_device->fuel_state[URHO]*velx;
    state(i, j, k,UMY)   = PeleC::h_prob_parm_device->fuel_state[URHO]*vely;
    state(i, j, k,UMZ)   = PeleC::h_prob_parm_device->fuel_state[URHO]*velz;
    state(i, j, k,UEINT) = PeleC::h_prob_parm_device->fuel_state[UEINT];
    state(i, j, k,UEDEN) = PeleC::h_prob_parm_device->fuel_state[UEINT] + PeleC::h_prob_parm_device->fuel_state[URHO] * (0.5 * (velx * velx + vely * vely + velz * velz));
    state(i, j, k,UTEMP) = PeleC::h_prob_parm_device->fuel_state[UTEMP];
    for (int n = 0; n < NUM_SPECIES; n++)
      state(i, j, k,UFS + n) = PeleC::h_prob_parm_device->fuel_state[UFS + n];
  }

}

void read_input_file(
  const std::string iname,
  int z_coord,
  int nx,
  int ny,
  int nz,
  int nscal,
  amrex::Vector<amrex::Real>& gridx_input,
  amrex::Vector<amrex::Real>& gridy_input,
  amrex::Vector<amrex::Real>& gridz_input,
  amrex::Vector<amrex::Real>& data);

void get_inputs(
  const std::string& iname,
   int z_coord);

void pc_prob_close();

void
init_composition(amrex::Real molefrac[NUM_SPECIES]);

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
