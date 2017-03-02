// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifdef DUNE_DYCAP_MULTISPAN_UTILITIES_HH
#warning ***** WARNING ***** multispantwophase.hh was already included ******
#endif

#ifndef DUNE_DYCAP_MULTISPANCF_UTILITIES_HH
#define DUNE_DYCAP_MULTISPANCF_UTILITIES_HH

#include<iostream>
#include<vector>
#include<map>
#include<list>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>

#include"time_step_manager.hh"
#include"../parameters/heleshawdomain.hh"


template<class GV, class Real>
class TwoPhaseParameter;

template <typename GV, typename Real>
class MultiSpanCFReaction : public ConvergenceAdaptiveTimeStepper<Real>, public TwoPhaseParameter<GV,Real>
{
public:
  typedef ConvergenceAdaptiveTimeStepper<Real> Base;
  typedef TwoPhaseParameter<GV,Real> Base2;



private:
  class SpanTwoPhaseParameter
  {
  public:
    const Real start;
    const Real end;
    const Real dt_min;
    const Real dt_max;
    const std::vector<Real> leftinlets;
    const Real leftinletsize;
    const Real leftflux;
    const std::vector<Real> rightinlets;
    const Real rightinletsize;
    const Real rightflux;

    const std::vector<Real> bottominlets;
    const Real bottominletsize;
    const Real bottomflux;
    const std::vector<Real> topinlets;
    const Real topinletsize;
    const Real topflux;


    SpanTwoPhaseParameter(Real start_,Real end_, Real dt_min_, Real dt_max_, std::vector<Real> leftinlets_, Real leftinletsize_, Real leftflux_, std::vector<Real> rightinlets_, Real rightinletsize_, Real rightflux_, std::vector<Real> bottominlets_, Real bottominletsize_, Real bottomflux_,std::vector<Real> topinlets_, Real topinletsize_, Real topflux_) DUNE_DEPRECATED
      : start(start_), end(end_), dt_min(dt_min_), dt_max(dt_max_), leftinlets(leftinlets_), leftinletsize(leftinletsize_), leftflux(leftflux_),rightinlets(rightinlets_), rightinletsize(rightinletsize_), rightflux(rightflux_), bottominlets(bottominlets_), bottominletsize(bottominletsize_), bottomflux(bottomflux_), topinlets(topinlets_), topinletsize(topinletsize_), topflux(topflux_)
    {}
  };

  typedef std::list<SpanTwoPhaseParameter> SpanTwoPhaseList;
  SpanTwoPhaseList span_twophase;
  typename SpanTwoPhaseList::const_iterator twophase_iterator;


  class SpanComputeParameter
  {
  public:
    const bool initial;
    const bool twophase;
    const bool O2liquid;
    const bool DOCliquid;
    const bool Ecoliliquid;
    const bool O2gas;
    const bool reaction;

    SpanComputeParameter(bool initial_, bool twophase_, bool O2liquid_, bool DOCliquid_, bool Ecoliliquid_, bool O2gas_, bool reaction_)
      : initial(initial_), twophase(twophase_), O2liquid(O2liquid_), DOCliquid(DOCliquid_), Ecoliliquid(Ecoliliquid_), O2gas(O2gas_), reaction(reaction_)
    {}
  };
  typedef std::list<SpanComputeParameter> SpanComputeList;
  SpanComputeList span_compute;
  typename SpanComputeList::const_iterator compute_iterator;

public:
  enum {dim=GV::Grid::dimension};

  MultiSpanCFReaction(const GV& gv_, const HeleshawDomain<dim>& dom_, const Dune::ParameterTree & param) DUNE_DEPRECATED
    : ConvergenceAdaptiveTimeStepper<Real>(), TwoPhaseParameter<GV,Real>(gv_, dom_, param), twophase_iterator(span_twophase.begin()), compute_iterator(span_compute.begin()), noflow(true), verbosity_level(1)
  {
    //  Dune::ParameterTreeParser::readINITree(filename,param,false);
    typedef typename Dune::ParameterTree::KeyVector KeyVector;
    const Dune::ParameterTree intervals = param.Dune::ParameterTree::sub("intervals");
    const KeyVector keyvector = intervals.getSubKeys();
    KeyVector::const_iterator it = keyvector.begin();
    KeyVector::const_iterator eit = keyvector.end();

    std::vector<Real> leftinlets_default =  param.get<std::vector<Real> >("default.sideinlets",std::vector<Real>(1, 0.0));
    Real leftinletsize_default =  param.get<Real>("default.sideinletsize",0.);
    std::vector<Real> rightinlets_default =  param.get<std::vector<Real> >("default.sideinlets",std::vector<Real>(1, 0.0));
    Real rightinletsize_default =  param.get<Real>("default.sideinletsize",0.);
    std::vector<Real> bottominlets_default =  param.get<std::vector<Real> >("default.bottominlets",std::vector<Real>(1, 0.0));
    Real bottominletsize_default =  param.get<Real>("default.bottominletsize",0.);
    std::vector<Real> topinlets_default =  param.get<std::vector<Real> >("default.topinlets",std::vector<Real>(1, 0.0));
    Real topinletsize_default =  param.get<Real>("default.topinletsize",0.);

    if (gv_.comm().rank()>0)
      verbosity_level = 0;

    for(; it!=eit; ++it){
      const Dune::ParameterTree sub = intervals.Dune::ParameterTree::sub(*it);

      SpanTwoPhaseParameter p(sub.get<Real>(std::string("start")),
                              sub.get<Real>(std::string("end")),
                              sub.get<Real>(std::string("dt_min")),
                              sub.get<Real>(std::string("dt_max")),
                              sub.get<std::vector<Real> >(std::string("leftinlets"), sub.get<std::vector<Real> >(std::string("sideinlets"),leftinlets_default)),
                              sub.get<Real>(std::string("leftinletsize"),sub.get<Real>(std::string("sideinletsize"),leftinletsize_default)),
                              sub.get<Real>(std::string("leftflux"), sub.get<Real>(std::string("sideflux"),0.)),
                              sub.get<std::vector<Real> >(std::string("rightinlets"), sub.get<std::vector<Real> >(std::string("sideinlets"),rightinlets_default)),
                              sub.get<Real>(std::string("rightinletsize"),sub.get<Real>(std::string("sideinletsize"),rightinletsize_default)),
                              sub.get<Real>(std::string("rightflux"), sub.get<Real>(std::string("sideflux"),0.)),

                              sub.get<std::vector<Real> >(std::string("bottominlets"),bottominlets_default),
                              sub.get<Real>(std::string("bottominletsize"),bottominletsize_default),
                              sub.get<Real>(std::string("bottomflux"),0.),
                              sub.get<std::vector<Real> >(std::string("topinlets"),topinlets_default),
                              sub.get<Real>(std::string("topinletsize"),topinletsize_default),
                              sub.get<Real>(std::string("topflux"),0.));

      SpanComputeParameter pc(sub.get<bool>(std::string("computeInitial" ),false),
                              sub.get<bool>(std::string("computeTwoPhase"),false),
                              sub.get<bool>(std::string("computeO2Liquid"),false),
                              sub.get<bool>(std::string("computeDOCLiquid"),false),
                              sub.get<bool>(std::string("computeEcoliLiquid"),false),
                              sub.get<bool>(std::string("computeO2Gas"),false),
                              sub.get<bool>(std::string("computeReaction"),false));


      if(span_twophase.size() && p.start != span_twophase.back().end)
        {DUNE_THROW(Dune::Exception,"Found inconsistent time spans. One ends at "
                    << span_twophase.back().end << " while the next begins at "
                    << p.start);}

      span_twophase.push_back(p);
      span_compute.push_back(pc);

      if (verbosity_level)
        {
          std::cout << "\nTime stepping span (" <<*it<< "): "<< p.start<<"-"<<p.end
                    <<":\n  dt_min = "<<p.dt_min
                    <<" dt_max = " <<p.dt_max
                    <<"\n  Side:"
                    << "\n  inlets = ";
          for(typename std::vector<Real>::size_type n = 0; n < p.leftinlets.size(); ++n)
            std::cout  << p.leftinlets[n] << " ";
          std::cout << "\n  flux = " << p.leftflux
                    << "  leftinletsize = " << p.leftinletsize
                    <<"\n  inlets = ";
          for(typename std::vector<Real>::size_type n = 0; n < p.rightinlets.size(); ++n)
            std::cout  << p.rightinlets[n] << " ";
          std::cout << "\n  flux = " << p.rightflux
                    << "  rightinletsize = " << p.rightinletsize
                    <<"\n  Bottom:"
                    << "\n  inlets = ";
          for(typename std::vector<Real>::size_type n = 0; n < p.bottominlets.size(); ++n)
            std::cout  << p.bottominlets[n] << " ";
          std::cout << "\n  flux = " << p.bottomflux
                    << "  inletsize = " << p.bottominletsize
                    <<"\n  Top:"
                    << "\n  inlets = ";
          for(typename std::vector<Real>::size_type n = 0; n < p.topinlets.size(); ++n)
            std::cout  << p.topinlets[n] << " ";
          std::cout << "\n  flux = " << p.topflux
                    << "  inletsize = " << p.topinletsize
                    << "\n";

          std::cout << "\n  computeInitial     = " << pc.initial
                    <<   "  computeTwoPhase    = " << pc.twophase
                    << "\n  computeO2Liquid    = " << pc.O2liquid
                    <<   "  computeDOCLiquid   = " << pc.DOCliquid
                    << "\n  computeEcoliLiquid = " << pc.Ecoliliquid
                    <<   "  computeO2Gas       = " << pc.O2gas
                    << "\n  ComputeReaction  = " << pc.reaction
                    << std::endl;
        }

    }
    twophase_iterator = span_twophase.begin();
    compute_iterator = span_compute.begin();

    Base::final_time = span_twophase.front().end;
    Base::time = span_twophase.front().start;
    Base::dt_eps = param.get<Real>("default.dt_eps");
    Base::dt_plot = param.get<Real>("default.dt_plot");
    Base::next_plot_time = Base::time;
    Base::dt_increase_ratio = param.get<Real>("default.dt_increase_rate");
    //   Base::dt = param.get<Real>("propagation.dt");

    Base2::incompressible = param.get<bool>("Setup.incompressible",false);

    Base2::setInlets(Base2::left,span_twophase.front().leftinlets,span_twophase.front().leftinletsize,span_twophase.front().leftflux);
Base2::setInlets(Base2::right,span_twophase.front().rightinlets,span_twophase.front().rightinletsize,span_twophase.front().rightflux);
    Base2::setInlets(Base2::top,span_twophase.front().topinlets,span_twophase.front().topinletsize,span_twophase.front().topflux);
    Base2::setInlets(Base2::bottom,span_twophase.front().bottominlets,span_twophase.front().bottominletsize,span_twophase.front().bottomflux);

  init = true;
  }

  /**
     \brief Provides the current time step size.
  */
  Real getTimeStepSize()
  {
    if(!Base::notified)
      return Base::dt;

    while(Base::time < twophase_iterator->start){
      if(twophase_iterator == span_twophase.begin())
        {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
      --twophase_iterator;
      --compute_iterator;
    }

    while(Base::time >= twophase_iterator->end){
      if (verbosity_level)
        std::cout << " timespan end " << twophase_iterator->end << std::endl;
      if(twophase_iterator == span_twophase.end())
        {
          if (verbosity_level)
            std::cout << "Last Time Span, time " << Base::time << " timespan end " << twophase_iterator->end << std::endl;
          {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
        }
      else
        {
          ++twophase_iterator;
          ++compute_iterator;
        }
      init = true;
    }
    Base::dt_min = twophase_iterator->dt_min;
    Base::dt_max = twophase_iterator->dt_max;
    Base::final_time = twophase_iterator->end;

    Base2::setInlets(Base2::left,twophase_iterator->leftinlets,twophase_iterator->leftinletsize,twophase_iterator->leftflux);
    Base2::setInlets(Base2::right,twophase_iterator->rightinlets,twophase_iterator->rightinletsize,twophase_iterator->rightflux);
    Base2::setInlets(Base2::top,twophase_iterator->topinlets,twophase_iterator->topinletsize,twophase_iterator->topflux);
    Base2::setInlets(Base2::bottom,twophase_iterator->bottominlets,twophase_iterator->bottominletsize,twophase_iterator->bottomflux);



    if (std::abs(twophase_iterator->leftflux) < 1.e-7 && std::abs(twophase_iterator->rightflux)<1.e-7 && std::abs(twophase_iterator->bottomflux) < 1.e-7 && std::abs(twophase_iterator->topflux)<1.e-7)
      noflow = true;
    else
      noflow = false;

    if (verbosity_level)
      std::cout << "Current time step limits: dt_min : " << Base::dt_min
                << ",  dt_max : " << Base::dt_max << " Time " << Base::time << " dt " << Base::getTimeStepSize() << " time_end " << Base::final_time << " final_span_time " << twophase_iterator->end  << std::endl;


    return Base::getTimeStepSize();
  }


  bool finalize()
  {
    return Base::time>=span_twophase.back().end-1.e-6;
  }

  const bool computeInitial()
  {
    if (!init)
      return false;
    else
      init = false;

    return compute_iterator->initial;
  }

  const bool computeTwoPhase()
  {
    return compute_iterator->twophase;
  }

  const bool computeTransportO2Liquid()
  {
    return compute_iterator->O2liquid;
  }

  const bool computeTransportDOCLiquid()
  {
    return compute_iterator->DOCliquid;
  }

  const bool computeTransportEcoliLiquid()
  {
    return compute_iterator->Ecoliliquid;
  }

  const bool computeTransportO2Gas()
  {
    return compute_iterator->O2gas;
  }

  const bool computeReaction()
  {
    return compute_iterator->reaction;
  }

  void set_dt(Real dt)
  {
    Base::dt = dt;
  }

  bool isNoFlow()
  {
    return noflow;
  }


private:
  bool init;
  bool noflow;
  unsigned int verbosity_level;

};
#endif


