// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_MULTISPAN_UTILITIES_HH
#define DUNE_DYCAP_MULTISPAN_UTILITIES_HH

#include<iostream>
#include<vector>
#include<map>
#include<list>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>
#include"time_step_manager.hh"



template <typename GV, typename Real>
class MultiSpanTwoPhase : public ConvergenceAdaptiveTimeStepper<Real>
{
public:
  typedef ConvergenceAdaptiveTimeStepper<Real> Base;

private:

  class SpanComputeParameter
  {
  public:
    const Real start;
    const Real end;
    const Real dt_min;
    const Real dt_max;
    std::map<std::string,bool> computeMap;
    int number;

    SpanComputeParameter(Real start_, Real end_,  Real dt_min_, Real dt_max_, std::map<std::string,bool> computeMap_, int number_)
      : start(start_), end(end_),dt_min(dt_min_), dt_max(dt_max_),
        computeMap(computeMap_), number(number_)
    {}
  };

  typedef std::list<SpanComputeParameter> SpanComputeList;
  SpanComputeList span_compute;
  typename SpanComputeList::iterator compute_iterator;

public:
  enum {dim=GV::Grid::dimension};

  MultiSpanTwoPhase(const GV& gv_, const Dune::ParameterTree & param)
    : Base(),  verbosity_level(1)  {



    typedef typename Dune::ParameterTree::KeyVector KeyVector;
    const Dune::ParameterTree intervals = param.Dune::ParameterTree::sub("intervals");
    const KeyVector keyvector = intervals.getSubKeys();
    KeyVector::const_iterator it = keyvector.begin();
    KeyVector::const_iterator eit = keyvector.end();

    if (gv_.comm().rank()>0)
      verbosity_level = 0;

    if (verbosity_level)
      std::cout << "creating MultiSpanTwoPhase " << std::endl;

       Real dt_min_default = param.get<Real>("InletsDefault.dt_min", 0.);
       Real dt_max_default = param.get<Real>("InletsDefault.dt_max", 1.e100);

    for(; it!=eit; ++it){
      const Dune::ParameterTree sub = intervals.Dune::ParameterTree::sub(*it);

      KeyVector subkey = sub.getValueKeys();
      KeyVector::iterator sit = subkey.begin();
      KeyVector::iterator seit = subkey.end();

      std::map<std::string,bool> computeMap;

      for(; sit!=seit; ++sit){
        if ( (*sit).find("compute")!= std::string::npos)
          {
            if (computeMap.count(*sit)<1)
              computeMap.insert(make_pair(*sit,sub.get<bool>(*sit)));
            else
              DUNE_THROW(Dune::RangeError, "Key '" << *sit << "' appears in the this interval 2x");
          }
      }


      SpanComputeParameter pc(sub.get<Real>(std::string("start")),
                              sub.get<Real>(std::string("end")),
                              sub.get<Real>(std::string("dt_min"), dt_min_default),
                              sub.get<Real>(std::string("dt_max"), dt_max_default),
                              computeMap,std::distance(keyvector.begin(),it));

      span_compute.push_back(pc);

      if (verbosity_level)
        {
          std::cout << "\nTime Span Interval from " << pc.start << " to " << pc.end << "\ncompute:"<< std::endl;

          std::string toReplace = "compute";

          for (auto iter = computeMap.begin(); iter != computeMap.end(); iter++)
            {
              std::string compute = iter->first;
              std::string replace = "compute";
              compute.replace(compute.find(replace),replace.length(),"");
            std::cout << compute << " " << iter->second << std::endl;
            }
          std::cout << "\n";
        }
      //else std::cout << "verbosity level is " << verbosity_level << " with rank " << gv_.comm().rank() << std::endl;

    }
    compute_iterator = span_compute.begin();

    Base::final_time = span_compute.front().end;
    Base::time = span_compute.front().start;
    Base::dt_max = ( span_compute.front().end -span_compute.front().start );
    Base::dt_eps = param.get<Real>("InletsDefault.dt_eps");
    Base::dt_plot = param.get<Real>("InletsDefault.dt_plot");
    Base::next_plot_time = Base::time;
    Base::increase_rate = param.get<Real>("InletsDefault.increase_rate");
    Base::success_time_factor = param.get<Real>("Timeloop.success_time_factor");
    Base::failure_time_factor = param.get<Real>("Timeloop.failure_time_factor");

    if (Base::success_time_factor < 1 || std::abs(Base::failure_time_factor)>=1.)
        DUNE_THROW(Dune::Exception,"Succes or failure time factor set to bad value!");
    //   Base::dt = param.get<Real>("propagation.dt");

    init = true;
    }

  /**
     \brief Provides the current time step size.
  */
  Real getTimeStepSize()
  {
    if(!Base::notified)
      return Base::dt;

    while(Base::time < compute_iterator->start){
      if (verbosity_level)
        std::cout << "Time " << Base::time << " Span start " << compute_iterator->start<< std::endl;
      if(compute_iterator == span_compute.begin())
        {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
      --compute_iterator;
    }

    while(Base::time >= compute_iterator->end){
      if (verbosity_level)
        std::cout << "compute timespan end " << compute_iterator->end << std::endl;
      if(compute_iterator == span_compute.end())
        {
          if (verbosity_level)
            std::cout << "Last Time Span, time " << Base::time << " timespan end " << compute_iterator->end << std::endl;
          {DUNE_THROW(Dune::Exception,"Time is out of any given interval " << Base::time);}
        }
      else
        {
          ++compute_iterator;
        }
      init = true;

    }
    Base::dt_min = compute_iterator->dt_min;
    Base::dt_max = std::min(compute_iterator->dt_max, compute_iterator->end - compute_iterator->start);
    Base::final_time = compute_iterator->end;


    /*
    if (std::abs(compute_iterator->sideflux) < 1.e-7 && std::abs(compute_iterator->bottomflux)<1.e-7)
      noflow = true;
    else
      noflow = false;
    */

    if (verbosity_level)
      std::cout << "Current time step limits: dt_min : " << Base::dt_min
                << ",  dt_max : " << Base::dt_max << " Time " << Base::time << " dt " << Base::getTimeStepSize() << " time_end " << Base::final_time << " final_span_time " << compute_iterator->end  << std::endl;


    return Base::getTimeStepSize();
  }


  bool finalize()
  {
    return Base::time>=span_compute.back().end-1.e-6;
  }


  const bool compute(std::string part)
  {
    std::string part2;
    if (part.find("compute") == std::string::npos)
      part="compute"+part;

    if (part == "computeInitial")
      {
        if (!init)
          return false;
        else
          init = false;
      }

    if ( (compute_iterator->computeMap.count(part)<1))
      return false;
    else return compute_iterator->computeMap.find(part)->second;
  }

  void set_dt(Real dt_)
  {
    Base::dt = dt_;
  }

  bool isNoFlow()
  {
    return noflow;
  }

  Real getDtMax()
  {
    return Base::dt_max;
  }

  Real getDtMin()
  {
    return Base::dt_min;
  }

  int getNumberOfSpans()
  {
    return span_compute.size();
  }

  size_t getSpanNumber()
  {
    return compute_iterator->number;
  }




private:
  bool init;
  bool noflow;
  unsigned int verbosity_level;

};

template <typename Real>
class MultiSpanTransport : public ConvergenceAdaptiveTimeStepper<Real>
{
public:
  typedef ConvergenceAdaptiveTimeStepper<Real> Base;

private:

  class SpanComputeParameter
  {
  public:
    const Real start;
    const Real end;
    const Real dt_min;
    const Real dt_max;
    std::map<std::string,bool> computeMap;

    SpanComputeParameter(Real start_, Real end_,  Real dt_min_, Real dt_max_, std::map<std::string,bool> computeMap_)
      : start(start_), end(end_),dt_min(dt_min_), dt_max(dt_max_),
        computeMap(computeMap_)
    {}
  };

  typedef std::list<SpanComputeParameter> SpanComputeList;
  SpanComputeList span_compute;
  typename SpanComputeList::const_iterator compute_iterator;

public:

  MultiSpanTransport(const Dune::ParameterTree & param, const size_t rank)
    : ConvergenceAdaptiveTimeStepper<Real>(),compute_iterator(span_compute.begin()), verbosity_level(1)
  {
    typedef typename Dune::ParameterTree::KeyVector KeyVector;
    const Dune::ParameterTree intervals = param.Dune::ParameterTree::sub("intervals");
    const KeyVector keyvector = intervals.getSubKeys();
    KeyVector::const_iterator it = keyvector.begin();
    KeyVector::const_iterator eit = keyvector.end();



    if (rank>0)
         verbosity_level = 0;

       Real dt_min_default = param.get<Real>("SpanDefault.dt_min", 0.);
       Real dt_max_default = param.get<Real>("SpanDefault.dt_max", 1.e100);

    for(; it!=eit; ++it){
      const Dune::ParameterTree sub = intervals.Dune::ParameterTree::sub(*it);

      const KeyVector subkey = sub.getValueKeys();
      auto sit = subkey.begin();
      auto seit = subkey.end();

      std::map<std::string,bool> computeMap;

      for(; sit!=seit; ++sit){
        if ( (*sit).find("compute")!= std::string::npos)
          {

            if (computeMap.count(*sit)<1)
              computeMap.insert(make_pair(*sit,sub.get<bool>(*sit)));
            else
              DUNE_THROW(Dune::RangeError, "Key '" << *sit << "' appears in the this interval 2x");
          }
      }


      SpanComputeParameter pc(sub.get<Real>(std::string("start")),
                              sub.get<Real>(std::string("end")),
                              sub.get<Real>(std::string("dt_min"), dt_min_default),
                              sub.get<Real>(std::string("dt_max"), dt_max_default),
                              computeMap);

      span_compute.push_back(pc);

      if (verbosity_level)
        {
          std::cout << "\nTime Span Interval from " << pc.start << " to " << pc.end << std::endl;
          for (auto iter = computeMap.begin(); iter != computeMap.end(); iter++)
            std::cout << iter->first << " " << iter->second << std::endl;

        }

    }
    compute_iterator = span_compute.begin();

    Base::final_time = span_compute.front().end;
    Base::time = span_compute.front().start;
    Base::dt_max = ( span_compute.front().end -span_compute.front().start );
    Base::dt_eps = param.get<Real>("SpanDefault.dt_eps");
    Base::dt_plot = param.get<Real>("SpanDefault.dt_plot");
    Base::next_plot_time = Base::time;
    init = true;
    }

  /**
     \brief Provides the current time step size.
  */
  Real getTimeStepSize()
  {
    if(!Base::notified)
      return Base::dt;

    while(Base::time < compute_iterator->start){
      if (verbosity_level)
        std::cout << "Time " << Base::time << " Span start " << compute_iterator->start<< std::endl;
      if(compute_iterator == span_compute.begin())
        {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
      --compute_iterator;
    }

    while(Base::time >= compute_iterator->end){
      if (verbosity_level)
        std::cout << "compute timespan end " << compute_iterator->end << std::endl;
      if(compute_iterator == span_compute.end())
        {
          if (verbosity_level)
            std::cout << "Last Time Span, time " << Base::time << " timespan end " << compute_iterator->end << std::endl;
          {DUNE_THROW(Dune::Exception,"Time is out of any given interval " << Base::time);}
        }
      else
        {
          ++compute_iterator;
        }
      init = true;

    }
    Base::dt_min = compute_iterator->dt_min;
    Base::dt_max = std::min(compute_iterator->dt_max, compute_iterator->end - compute_iterator->start);
    Base::final_time = compute_iterator->end;


    if (verbosity_level)
      std::cout << "Current time step limits: dt_min : " << Base::dt_min
                << ",  dt_max : " << Base::dt_max << " Time " << Base::time << " dt " << Base::getTimeStepSize() << " time_end " << Base::final_time << " final_span_time " << compute_iterator->end  << std::endl;


    return Base::getTimeStepSize();
  }


  bool finalize()
  {
    return Base::time>=span_compute.back().end-1.e-6;
  }


  void set_dt(Real dt)
  {
    Base::dt = dt;
  }

  Real getDtMax()
  {
    return Base::dt_max;
  }

  Real getDtMin()
  {
    return Base::dt_min;
  }

  Real getFinalTime()
  {
    return span_compute.back().end;;
  }



private:
  bool init;
  bool noflow;
  unsigned int verbosity_level;

};

#endif
