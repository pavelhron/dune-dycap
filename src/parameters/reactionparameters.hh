// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_REACTIONPARAMETER_HH
#define DUNE_DYCAP_REACTIONPARAMETER_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>

#include"twophaseparameters.hh"
#include<dune/dycap/models/reactionop.hh>


template<typename GV, typename RF, typename TP,  typename SLDGF, typename VEL>
class ReactionParameter
  : public Dune::PDELab::ReactionParameterInterface<Dune::PDELab::ReactionParameterTraits<GV,RF>,
                                                    ReactionParameter<GV,RF, TP, SLDGF, VEL> >
{

public:

  //! Types related to the Traits and TwoPhaseParameter
  //! @{
  typedef Dune::PDELab::ReactionParameterTraits<GV,RF> Traits;
  //typedef TwoPhaseParameter<GV,RF> TP;
  //! @}

  enum {dim=GV::Grid::dimension};

  //! constructor
  ReactionParameter(const GV& gv_, const TP& tp_, const SLDGF& solddgf_, const SLDGF& snewdgf_, const VEL& vel_) DUNE_DEPRECATED:
    gv(gv_),
    tp(tp_),
    param(tp.getParam()),  // parameter class
    solddgf(solddgf_), snewdgf(snewdgf_), // old and new liquid saturation
    vel(vel_),

    growthmodeltype(param.sub("Microorganism").get<std::string>("model")),
    exchangetype(param.sub("Microorganism").get<std::string>("exchangetype")),
    xmax(param.sub("Microorganism").get<RF>("xmax",1e3)),
    mumax(param.sub("Microorganism").get<RF>("mumax")/3600.),
    ks(param.sub("Microorganism").get<RF>("Ks")),
    ko(param.sub("Microorganism").get<RF>("Ko")),
    ys(param.sub("Microorganism").get<RF>("Ys")),
    yo(param.sub("Microorganism").get<RF>("Yo")),
    mo(param.sub("Microorganism").get<RF>("mo")/3600.),
    rd(param.sub("Microorganism").get<RF>("Rd")/3600.),
    rb(param.sub("Microorganism").get<RF>("rb")),
    kappa(param.sub("Microorganism").get<RF>("kappa",1.)),



    mumax_an(param.sub("Microorganism").get<RF>("mumax_an")/3600.),
    ks_an(param.sub("Microorganism").get<RF>("Ks_an")),
    ys_an(param.sub("Microorganism").get<RF>("Ys_an")),
    smin(param.sub("Microorganism").get<RF>("smin",-1.0)),

    eta(param.sub("Adhesion").get<RF>("eta",-1.0)),
    dg(param.sub("Adhesion").get<RF>("dg",-1.0)),
    kdet(param.sub("Adhesion").get<RF>("Kdet",-1.0)),
    rhobulk(param.sub("Adhesion").get<RF>("rhobulk",-1.0)),



    k_h(param.sub("Microorganism").get<RF>("K_H")),
    lambda(param.sub("Microorganism").get<RF>("lambda",0)/3600),
    time(0.), t0(0), dt(1e6)

  {
  }

  //! porosity
  typename Traits::RangeFieldType
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return tp.phi(e,x);
  }

  //! new liquid saturation
  typename Traits::RangeFieldType
  slnew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SLDGF::Traits::RangeType y;
    snewdgf.evaluate(e,x,y);
    return y;
  }

  //! old liquid saturation
  typename Traits::RangeFieldType
  slold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SLDGF::Traits::RangeType y;
    solddgf.evaluate(e,x,y);
    return y;
  }

  //! old liquid saturation
  typename Traits::RangeFieldType
  velocity (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename VEL::Traits::RangeType y;
    vel.evaluate(e,x,y);
    return y.two_norm();
  }


  //! liquid saturation (time dependent average between new and old values)
  typename Traits::RangeFieldType
  sl (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    RF factor = (time-t0)/dt;
    return (1.-factor)*slold(e,x) + factor*slnew(e,x);
  }

   //! maximum specific growth rate
  typename Traits::RangeFieldType
  Xmax (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return xmax;
  }

  //! maximum specific growth rate
  typename Traits::RangeFieldType
  Mumax (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mumax;
  }

  //! maximum specific growth rate
  typename Traits::RangeFieldType
  Smin (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (smin==-1)
      DUNE_THROW(Dune::Exception, "smin must be set in parameter file");
    return smin;
  }

  //! DOC half saturation constant
  typename Traits::RangeFieldType
  Ks (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return ks;
  }


  //! oxygen half saturation constant
  typename Traits::RangeFieldType
  Ko (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return ko;
  }

  //! death rate
  typename Traits::RangeFieldType
  Rd (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return rd;
  }

  //! growth yield coefficient based on DOC utilized
  typename Traits::RangeFieldType
  Ys (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return ys;
  }

  //! growth yield coefficient based on oxygen utilized
  typename Traits::RangeFieldType
  Yo (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return yo;
  }

  //! growth yield coefficient based on oxygen utilized
  typename Traits::RangeFieldType
  Mo (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mo;
  }


  //! rbonent in model
  typename Traits::RangeFieldType
  Rb (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return rb;
  }

  //! rbonent in model
  typename Traits::RangeFieldType
  Kappa (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return kappa;
  }

    //! maximum specific growth rate
  typename Traits::RangeFieldType
  Mumax_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return mumax_an;
  }

  //! DOC half saturation constant
  typename Traits::RangeFieldType
  Ks_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return ks_an;
  }


  //! growth yield coefficient based on DOC utilized
  typename Traits::RangeFieldType
  Ys_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return ys_an;
  }


  //! Henry's Constant
  typename Traits::RangeFieldType
  K_H(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return k_h;
  }

  //! Henry's Constant
  typename Traits::RangeFieldType
  K_H() const
  {
    return k_h;
  }


  typename Traits::RangeFieldType
  Eta() const {
    if (eta<0)
      DUNE_THROW(Dune::Exception, "eta smaller than zero");
  return eta;}


  typename Traits::RangeFieldType
  Dg() const
  {
    if (dg<0)
      DUNE_THROW(Dune::Exception, "dg smaller than zero");
    return dg;
  }


  typename Traits::RangeFieldType
  Kdet() const
  {
    if (kdet<0)
      DUNE_THROW(Dune::Exception, "kdet smaller than zero");
    return kdet;}


  typename Traits::RangeFieldType
  Rhobulk() const {
    if (rhobulk<0)
      DUNE_THROW(Dune::Exception, "rhobulk smaller than zero");
  return rhobulk;}


  //! oxygen exchange Constant
  typename Traits::RangeFieldType
  k_l(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return lambda;
  }

  //! oxygen concentration in water
  typename Traits::RangeFieldType
  c_O2_liquid (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
               typename Traits::RangeFieldType c) const
  {
    RF factor = (time-t0)/dt;
    return c*phi(e,x)*( (1.-factor)*slold(e,x) + factor*slnew(e,x));
  }

  //! oxygen concentration in air
  typename Traits::RangeFieldType
  c_O2_gas (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType c) const
  {
    RF factor = (time-t0)/dt;
    return  c*phi(e,x)*( (1.-factor)*(1.-slold(e,x)) + factor*(1.-slnew(e,x)));
  }

  //! DOC concentration
  typename Traits::RangeFieldType
  c_DOC (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
	 typename Traits::RangeFieldType c) const
  {
    RF factor = (time-t0)/dt;
    return  c*phi(e,x)*( (1.-factor)*slold(e,x) + factor*slnew(e,x));
  }

  //! Ecoli concentration
  typename Traits::RangeFieldType
  c_Ecoli (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
	   typename Traits::RangeFieldType c) const
  {
    RF factor = (time-t0)/dt;
    return c*phi(e,x)*( (1.-factor)*slold(e,x) + factor*slnew(e,x));;
  }


  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

  //! to be called once before each time step
  void preStep (RF time_, RF dt_, int stages)
  {
    t0 = time_;
    dt = dt_;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

  const std::string GrowthModelType()
  {
    return growthmodeltype;
  }

  const std::string ExchangeType()
  {
    return exchangetype;
  }


  void setExchangeType(std::string type_)
  {
    exchangetype = type_;
  }

  const TP& getTP()
  {
    return tp;
  }

private:
  const GV& gv;
  const TP& tp;

  // parameter class
  const Dune::ParameterTree & param;

  const SLDGF& solddgf;
  const SLDGF& snewdgf;
  const VEL& vel;

  const std::string growthmodeltype;
  std::string exchangetype;

  // DDM parameters
  const RF xmax;
  const RF mumax;
  const RF ks;
  const RF ko;
  const RF ys;
  const RF yo;
  const RF mo;
  const RF rd;

  const RF rb;
  const RF kappa;

  const RF mumax_an;
  const RF ks_an;
  const RF ys_an;
  const RF smin;

  const RF eta;
  const RF dg;
  const RF kdet;
  const RF rhobulk;

  const RF k_h;
  const RF lambda;



  // time values
  RF time, t0, dt;

};

template<typename GV, typename RF>
class IronReactionParameter
  : public Dune::PDELab::ReactionParameterInterface<Dune::PDELab::ReactionParameterTraits<GV,RF>,
                                                    IronReactionParameter<GV,RF> >
{

public:

  //! Types related to the Traits and TwoPhaseParameter
  //! @{
  typedef Dune::PDELab::ReactionParameterTraits<GV,RF> Traits;
  //! @}

  enum {dim=GV::Grid::dimension};

  //! constructor
  IronReactionParameter(const GV& gv_, const Dune::ParameterTree & param_):
    gv(gv_),
    param(param_),  // parameter class


    time(0.), t0(0), dt(1e6)

  {
  }

  //! porosity
  /* typename Traits::RangeFieldType
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return tp.phi(e,x);
    }*/



  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

  //! to be called once before each time step
  void preStep (RF time_, RF dt_, int stages)
  {
    t0 = time_;
    dt = dt_;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }



private:
  const GV& gv;
  // parameter class
  const Dune::ParameterTree & param;

  // DDM parameters
  const RF khom;
  const RF Ks;
  const RF Ns;
  const RF Sa;
  const RF A;

  // time values
  RF time, t0, dt;

};


template<typename GV, typename RF, typename TP,  typename SLDGF, typename VEL>
class ReactionParameterEcoli
  : public Dune::PDELab::ReactionParameterInterface<Dune::PDELab::ReactionParameterTraits<GV,RF>,
                                                    ReactionParameterEcoli<GV,RF, TP, SLDGF, VEL> >
{

public:

  //! Types related to the Traits and TwoPhaseParameter
  //! @{
  typedef Dune::PDELab::ReactionParameterTraits<GV,RF> Traits;
  //typedef TwoPhaseParameter<GV,RF> TP;
  //! @}

  enum {dim=GV::Grid::dimension};

  //! constructor
  ReactionParameterEcoli(const GV& gv_, const TP& tp_, const SLDGF& solddgf_, const SLDGF& snewdgf_, const VEL& vel_) :
    gv(gv_),
    tp(tp_),
    param(tp.getParam()),  // parameter class
    solddgf(solddgf_), snewdgf(snewdgf_), // old and new liquid saturation
    vel(vel_),
    time(0.), t0(0), dt(1e6)

  {
  }

  //! porosity
  typename Traits::RangeFieldType
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return tp.phi(e,x);
  }

  //! new liquid saturation
  typename Traits::RangeFieldType
  slnew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SLDGF::Traits::RangeType y;
    snewdgf.evaluate(e,x,y);
    return y;
  }

  //! old liquid saturation
  typename Traits::RangeFieldType
  slold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SLDGF::Traits::RangeType y;
    solddgf.evaluate(e,x,y);
    return y;
  }

  //! old liquid saturation
  typename Traits::RangeFieldType
  velocity (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename VEL::Traits::RangeType y;
    vel.evaluate(e,x,y);
    return y.two_norm();
  }


  //! liquid saturation (time dependent average between new and old values)
  typename Traits::RangeFieldType
  sl (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    RF factor = (time-t0)/dt;
    return (1.-factor)*slold(e,x) + factor*slnew(e,x);
  }


  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

  //! to be called once before each time step
  void preStep (RF time_, RF dt_, int stages)
  {
    t0 = time_;
    dt = dt_;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

  const TP& getTP()
  {
    return tp;
  }

private:
  const GV& gv;
  const TP& tp;

  // parameter class
  const Dune::ParameterTree & param;

  const SLDGF& solddgf;
  const SLDGF& snewdgf;
  const VEL& vel;

  // time values
  RF time, t0, dt;

};


#endif
