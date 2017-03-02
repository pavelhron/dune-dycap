// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARAMETERESTIMATION_UNITS_HH
#define DUNE_DYCAP_PARAMETERESTIMATION_UNITS_HH


#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <stdlib.h>
#include <map>

#include"compiler.h"

class Units
{
public:
  typedef enum {OD,ZZ} UNITS;
private:
  double temperature;
  double weight;
  double odinit;
  double cellinit;
  UNITS system;



public:

  virtual double Temperature()
  {
    return temperature;
  }

  // return weight of ecolis [mg cells/ml water]
  virtual double getWeight(double value)
  {
    if (temperature<0. || weight<0.0)
      {
        std::cout << "Temperature and weight need to be set" << std::endl;
        throw error_abort();
      }

    if (system==OD)
      {
        if (odinit<0.0 || cellinit<0.0)
          {
            std::cout << "ODinit and cellinit need to be set" << std::endl;
            throw error_abort();
          }
      }


    if (system==OD)
      return cellinit/odinit*value*weight*1000.;
    else
      return value*weight*1000.;
  }


  void setTemperature(double value)
  {
    temperature=value;
  }

  void setWeight(double value)
  {
    weight=value;
  }

  void setODinit(double value)
  {
    odinit=value;
  }

  void setCellInit(double value)
  {
    cellinit=value;
  }

  void setSystem(std::string unitSystem)
  {

    if (unitSystem=="OD")
      system = OD;
    else if (unitSystem=="ZZ")
      system = ZZ;
    else
      {
        std::cout << "Unknown unit system " << unitSystem << ". Valid unit systems are ZZ and OD" << std::endl;
        throw error_abort();
      }

    std::cout << "new initsystem is " << system << std::endl;
  }


  UNITS GetSystem()
  {
    return system;
  }

  Units() : temperature(-1.0),
            weight(-1.0),
            odinit(-1.0),
            cellinit(-1.0)
  {
  }

  virtual ~Units()
  {};
};

Units UNITS;

#endif
