// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifndef DUNE_DYCAP_ADHESION_ESTIMATION_FITMODEL_HH
#define DUNE_DYCAP_ADHESION_ESTIMATION_FITMODEL_HH

#include <src/estimation/fitclass.hh>
#include <dune/common/exceptions.hh>
#include <src/estimation/utilities.hh>
#include <src/estimation/fit_ecoli_model_base.hh>



template<typename RModel>
class FitModel : public FitEcoliModelBase
{

private:
  typedef FitEcoliModelBase Base;

  /** class containint the anaerobic model
   *  and all the functionality:
   *  set initial solution
   *  set initial parameters
   *  compute solution
   *  update parameters
   */

  RModel& model;

  //explicit using base class members (not neccessary)
  using Base::curParameter_;
  using Base::curStringParameter_;
  using Base::fitParameter_;



private:
  std::string fname;

  //replace string in str from to
bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

public:

  //! constructor
  FitModel(RModel& model_, string parameterFileName, string dataFile)
    :Base(parameterFileName, dataFile),
     model(model_)
  {
    fname = parameterFileName;
    replace(fname,"parameters","results");
    fname.resize(fname.length()-4);
    fname.append(".sol");


    std::cout << "Construct class FitModel" << std::endl;
    try {

      for (map<string, parameterType>::iterator param = curParameter_.begin(); param
             != curParameter_.end(); ++param){
        model.setParamValue(param->second.name, param->second.value);
        std::cout << "parameter " <<  param->second.name << " value " << param->second.value << std::endl;}
    }

    catch (Dune::Exception e) {
      std::cout << e << std::endl;
    }
  }

  ~FitModel()
  {
  }

  void SolutionOutput(FitDataClass &data, const FLOAT timesteps=100.)
  {
    std::ofstream s;
    std::cout << "Writting solution to file " << fname << std::endl;
    s.open(fname.c_str());
    if (!s.is_open())
      DUNE_THROW(Dune::IOError, "File " << fname <<  " was not opened!");
    // set all parameters to the growth model
    for(size_t i=0; i<fitParameter_.size(); ++i)
      model.setParamValue(fitParameter_[i].name, fitParameter_[i].value);

    // vector with time steps
    std::vector<FLOAT> tt;
    // fill time vector
    for(auto it=data.GetTimes().begin();it!=data.GetTimes().end();++it){
      tt.push_back(*it);
    }

    int numberofexperiments = data.getNumberOfExperiments();
    // int numberofactiveexperiments = data.GetNumValues() /  tt.size();

    // loop over all experiments
    for(int k = 0; k<numberofexperiments; k++)
      {
        model.printSolution(s,timesteps);
      }
    s.close();
  }

private:


  /** \brief function describing the Model
   *
   * \param data data to be fitted (all the set of data)
   * \param yVal vector with computed values, should be computed
   *
   *
   */
  void Model(FitDataClass &data, vector<FLOAT> &yVal)
  {

    // set all parameters to the growth model
    for(size_t i=0; i<fitParameter_.size(); ++i)
      model.setParamValue(fitParameter_[i].name, fitParameter_[i].value);

    model.reset();

    // vector with time steps
    std::vector<FLOAT> tt;
    // fill time vector
    for(auto it=data.GetTimes().begin();it!=data.GetTimes().end();++it){
      tt.push_back(*it);
    }

    int numberofexperiments = data.getNumberOfExperiments();
    // int numberofactiveexperiments = data.GetNumValues() /  tt.size();

    // loop over all experiments
    for(int k = 0; k<numberofexperiments; k++)
      {
        // set parameters c
        // setParameters(data,k);

        // initial values
        auto measurement = data.GetMeasurements(*tt.begin());
        for(auto element = measurement.begin(); element!=measurement.end(); element++)
          {
            if (element->active && element->localindex==k)
              {
#if TWOCOMPONENTS
                // flow is 1,2
                yVal[element->index] = model.getOutflow(k+1);
#else
                yVal[element->index] = model.getOutflow();
#endif
              }
          }

        // time loop
        for(auto tit=tt.begin();tit!=tt.end()-1;++tit)
          {
            // solve model
            auto dt = *(tit+1)-*tit;
#if TWOCOMPONENTS
                     model.apply(*tit,dt,k+1);
#else
                     model.apply(*tit,dt);
#endif


            auto measurement = data.GetMeasurements(*(tit+1));
            for(auto element = measurement.begin(); element!=measurement.end(); element++)
              {
                if (element->active && element->localindex==k)
                  {
#if TWOCOMPONENTS
                    yVal[element->index] = model.getOutflow(k+1);
#else
                    yVal[element->index] = model.getOutflow();
#endif
                  }
              } // measurements and model evaluation
          } // time loop

      } // experiments loop
  }

};

#endif /* FITECOLIMODEL_HH_ */
