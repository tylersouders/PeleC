#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  // Initialize values to possibly be changed via input file
  amrex::Real p_input       = 1013250.0;      //Pressure Pa
  amrex::Real T_input       = 0.0;            //Temperature K
  amrex::Real rho_input     = 0.00116;      //Density g/cm^3 - I believe this is redundant
  amrex::Real eint_input    = 0.0;         //Internal Energy input - I believe this is redundant
  amrex::Real rho_air       = 0.00116;     //Second set of variables for air regions 
  amrex::Real eint_air      = 0.0;  
  amrex::Real u_input       = 1500.0;         //u-direction velocity cm/s
  amrex::Real v_input       = 0.0;            //y-direction velocity cm/s
  amrex::Real w_input       = 0.0;            //z-direction velocity cm/s
  amrex::Real phi_input     = 0.0;          //Equivalence Ratio - Default to nonreacting

  amrex::Real fuel_y_limit  = 100.0;          //y-bounds for propane
  amrex::Real fuel_x_limit  = 100.0;          //x-max for propane initialization
  amrex::Real init_pert_mag = 0.0; // Decimal for initial perturbation percent

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac_air = {0.0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES> molefrac_air = {0.0};
};

struct ProbParmHost
{
  ProbParmHost() {}
};

#endif
