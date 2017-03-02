// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_PM_HELESHAWDOMAIN_HH
#define DUNE_PM_HELESHAWDOMAIN_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>


#include<src/physics/hydraulicparameters.hh>
#include<src/physics/physical_chemistry.hh>


//! \brief physical constants

//24.789 598(42) dm3/mol at 25 °C

//! atmospheric pressure [Pa]
const double patm = Dune::PM::PhysicalChemistry<double>::STANDARD_AIR_PRESSURE;

//!  ideal gas constant [N m / (mol * K)] = [J / (mol * K)]
const double R = Dune::PM::PhysicalChemistry<double>::IDEAL_GAS_CONSTANT;

//! specific gas constant [J / (kg * K)], \f$ R_s = \frac{R}{M}]\f$, where \f$M\f$ is molar mass of the gas
//! molar_mass_air =   0.026
//! ideal_gass_constant = 8.314472
const double R_S =
  /* Dune::PM::PhysicalChemistry<double>::IDEAL_GAS_CONSTANT / Dune::PM::PhysicalChemistry<double>::MOLAR_MASS_AIR */
  287.58; // wikipedia

const double M_g =  0.028911858;

//! water density [kg/m^3]
const double Rho_l = Dune::PM::PhysicalChemistry<double>::DENSITY_WATER ;

//! molar mass (water) [kg/mol], 0.018015
const double M_l = Dune::PM::PhysicalChemistry<double>::MOLAR_MASS_WATER;

//! \brief Heleshaw Domain description
template<int dim>
struct HeleshawDomain
{
  HeleshawDomain (const Dune::ParameterTree & param)
  {

    // read parameters
    height = param.get<double>("height");
    width  = param.get<double>("width");
    depth  = param.get<double>("depth");

    lense_width_min  = param.get<double>("lensexmin",0.0);
    lense_width_max  = param.get<double>("lensexmax",0.0);
    lense_height_min  = param.get<double>("lenseymin",0.0);
    lense_height_max  = param.get<double>("lenseymax",0.0);
    refine = param.get<int>("refine");
    nx     = param.get<int>("nx");
    ny     = param.get<int>("ny");
    if (param.get<bool>("1d",false))
      {
        ny =1;
        height=width/nx/refine;
      }
    nz = 1;
    if (dim == 3)
      nz     = param.get<int>("nz");
  }
  double height;
  double width;
  double depth;
  double lense_width_min;
  double lense_width_max;
  double lense_height_min;
  double lense_height_max;
  int refine;
  int nx;
  int ny;
  int nz;
};
#endif
