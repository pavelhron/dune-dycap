#pragma once
#include <dune/pdelab/instationary/onestepparameter.hh>
#include <dune/pdelab/instationary/onestep.hh>


namespace Dune {
  namespace PDELab {

    template<typename RF>
    class TimeSteppingMethods
    {
    public:
      TimeSteppingMethods()
      {
        methods.clear();
        // explicit methods
        methods.insert(std::make_pair("ExplicitEuler",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::ExplicitEulerParameter<RF>())));
        methods.insert(std::make_pair("Heun",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::HeunParameter<RF>())));
        methods.insert(std::make_pair("Shu3",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Shu3Parameter<RF>())));
        methods.insert(std::make_pair("RK4",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::RK4Parameter<RF>())));
        // implicit methods
        methods.insert(std::make_pair("ImplicitEuler",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::ImplicitEulerParameter<RF>())));
        methods.insert(std::make_pair("OneStepTheta",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::OneStepThetaParameter<RF>(0.5))));
        methods.insert(std::make_pair("Alexander2",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Alexander2Parameter<RF>())));
        methods.insert(std::make_pair("FractionalStep",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::FractionalStepParameter<RF>())));
        methods.insert(std::make_pair("Alexander3",std::shared_ptr<TimeSteppingParameterInterface<RF> >(new Dune::PDELab::Alexander3Parameter<RF>())));
      }

      template<typename OSM>
      void setTimestepMethod(OSM& osm, std::string method_name)
      {

        if (methods.count(method_name) )
          {
            std::shared_ptr<TimeSteppingParameterInterface<RF> > method = methods.find(method_name)->second;
            osm.setMethod(*method);
          }
        else  DUNE_THROW(Exception,"method " << method_name << "can NOT be found");
      }

    private:
      std::map<std::string, std::shared_ptr<TimeSteppingParameterInterface<RF> > > methods;
    };
  }
}
