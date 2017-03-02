#pragma once

#include <dune/pdelab/instationary/imexonestepparameter.hh>

namespace Dune {
  namespace PDELab {

    template<typename RF>
    class IMEXTimeSteppingMethods
    {
    public:
      IMEXTimeSteppingMethods()
      {
        methods.clear();
        // explicit methods
        methods.insert(std::make_pair("IMEXEuler",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXEulerParameter<RF>())));
        methods.insert(std::make_pair("IMEXTrapez",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXTrapezParameter<RF>())));
                methods.insert(std::make_pair("IMEXTheta",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXThetaParameter<RF>(0.5))));
                methods.insert(std::make_pair("IMEXAscher3",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXAscher3Parameter<RF>())));
                methods.insert(std::make_pair("IMEXAlexander2",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXAlexander2Parameter<RF>())));
                methods.insert(std::make_pair("IMEXPareschi2",std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> >(new Dune::PDELab::IMEXPareschi2Parameter<RF>())));
      }

      template<typename OSM>
      void setTimestepMethod(OSM& osm, std::string method_name)
      {

        if (methods.count(method_name) )
          {
            std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> > method = methods.find(method_name)->second;
            osm.setMethod(*method);
          }
        else  DUNE_THROW(Exception,"method " << method_name << "can NOT be found");
      }

    private:
      std::map<std::string, std::shared_ptr<IMEXTimeSteppingParameterInterface<RF> > > methods;
    };

  }
}
