#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Vector.H>

struct ProbParmDevice
{
  // Ambient Settings
  amrex::Real p_atm = 1013250.0;
  amrex::Real T_atm = 298.0;
  amrex::Real rho_atm = 0.0;
  amrex::Real eint_atm = 0.0;
  amrex::Real u_atm = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_atm;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac_atm;
  // Total Pressure and Temperature
  amrex::Real p0 = 1013250.0;
  amrex::Real T0 = 298.0;
  // Inlet conditions computed from choking requirement 
  amrex::Real p_in = 0.0;
  amrex::Real T_in = 0.0;
  amrex::Real rho_in = 0.0;
  amrex::Real eint_in = 0.0;
  amrex::Real u_in = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_in;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac_in;
  // Time to "shut off" the flow
  amrex::Real t_shutoff = 999999.9;
  amrex::Real outlet_mach_lim = 0.95;
};

struct ProbParmHost
{
  // amrex::Vector<amrex::Real> h_k;
  // amrex::Vector<amrex::Real> h_p;
  // amrex::Vector<amrex::Real> h_q;
  // amrex::Vector<amrex::Real> h_omega;
  // amrex::Gpu::DeviceVector<amrex::Real> hG_k;
  // amrex::Gpu::DeviceVector<amrex::Real> hG_p;
  // amrex::Gpu::DeviceVector<amrex::Real> hG_q;
  // amrex::Gpu::DeviceVector<amrex::Real> hG_omega;
  // ProbParmHost() 
  // : hG_k(0), hG_p(0), hG_q(0), hG_omega(0)
  // {
  // }
};

#endif
