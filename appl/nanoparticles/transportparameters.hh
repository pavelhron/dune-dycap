// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifdef DUNE_DYCAP_TRANSPORTPARAMETERS_HH
#warning ***** WARNING ***** transportparameters.hh was already included ******
#endif

#ifndef DUNE_DYCAP_TRANSPORTPARAMETERS_HH
#define DUNE_DYCAP_TRANSPORTPARAMETERS_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>

#include<src/physics/hydraulicparameters.hh>
#include<src/physics/physical_chemistry.hh>
#include<src/models/componenttransportop.hh>

#include <src/utilities/inlets_utilities.hh>

/** \brief Class to define the boundary condition types
 */
struct ConvectionDiffusionBoundaryConditions
{
  enum Type { Dirichlet=1, Neumann=-1, Outflow=-2, None=-3 }; // BC requiring constraints must be >0 if
  // constraints assembler coming with PDELab is used
};

//==============================================================================
// Problem definition : component transport
//==============================================================================

//! Transport in water phase
template<typename GV, typename RF>
class TransportLiquid :
  public Dune::PDELab::ModifiedTransportSpatialParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>,
                                                                  TransportLiquid<GV,RF> >,
  public Dune::PDELab::ModifiedTransportTemporalParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>,
                                                                   TransportLiquid<GV,RF> >
{
  enum {dim=GV::Grid::dimension};

public:
  typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;
  typedef Dune::Dycap::TransportInlets<Traits,RF> Inlets;
  typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;

  TransportLiquid (const GV& gv_, const Dune::ParameterTree& param,  const std::string cname, const std::string subintervals="intervals")
    : gv(gv_),
      time(0.),
      phi(param.sub("Setup").template get<RF>("phi")),
      velocityx(param.sub(cname).template get<RF>("velocityx")),
      velocityy(param.sub(cname).template get<RF>("velocityy")),
      Dt(param.sub(cname).template get<RF>("D")),
      bg_concentration(param.sub(cname).template get<RF>("input")),
      inlets(gv,param,cname,subintervals),
      advection(true),
      diffusion(true)

  {}

  /* phi .. porosity
     velocityx ... velocity in x direction
     velocityy ... velocity in y direction (in 2D)
     Dt ... diffusion
     bg_concentration ... concentration in background (e.g. for oxygen)
  */

  //! capacity function
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (dt<=0)
      DUNE_THROW(Dune::Exception, "dt for component transport was not specified in transportparameters, current dt is " << dt);
    RF s_new = snew(e,x);
    RF s_old = sold(e,x);
    RF factor = (tend-time)/dt;
    return phi*(factor*s_old + (1.-factor)*s_new )*(1+riso);
  }

  //! saturation at new time level
  typename Traits::RangeFieldType
  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  // set Retardation
  void setRiso(RF iso_) {riso=iso_;}

  //! saturation at old time level
  typename Traits::RangeFieldType
  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 1.0;
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType velo(0.0);
    if (!advection) return velo;

    velo[0]=velocityx;
    if (dim>1)
      velo[1]=velocityy;

    return velo;
  }

  //! tensor permeability
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (!diffusion)
      return 0.0;

    return Dt*sold(e,x)*phi;
  }

  //! source/reaction term
  typename Traits::RangeFieldType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 0.0;
  }

  //! boundary condition type function
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x)
  {
    //typename Traits::RangeType global = is.geometry().global(x);

    typename Traits::RangeType velo(0.0);
    if (!advection)
      {
        advection=true;
        velo = v(is.inside(), is.geometryInInside().center());
        advection=false;
      }
    else
      velo=v(is.inside(), is.geometryInInside().center());


    RF vn = velo*is.centerUnitOuterNormal();

    if (vn>0.)
      return BCType::Outflow;

    if (std::abs(vn)<1.e-10)
      return BCType::Neumann;

    if (vn<1e-10)
      return BCType::Dirichlet;

    return BCType::Neumann;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return inlets.evaluate(is,x);
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

  void setTimeTarget(RF time_, RF dt_)
  {
    tend = time_;
    dt = dt_;
    setTime(time_-dt);
  }

  //! to be called once before each time step
  void preStep (RF time_, RF dt_, int stages)
  {
    inlets.setTime(time_+dt_);
  }


  void onlyAdvection()
  {
    reset();
    diffusion=false;
  }

  void onlyDiffusion()
  {
    reset();
    advection=false;
  }

  void reset()
  {
    advection=true;
    diffusion=true;
  }

  void setDiffusion(RF Dt_)
  {
    Dt=Dt_;
  }

  void setVelocity(RF velocity_)
  {
    velocityx=velocity_;
  }

  const RF getBgConcentration() const
  {
    return bg_concentration;
  }

  inline const bool computeReaction()
  {
    return true;
  }


private:
  const GV& gv; // store access to two phase parameters
  RF time, tend, dt;
  RF phi;
  RF velocityx, velocityy;
  RF Dt;
  const RF bg_concentration;
  Inlets inlets;
  bool advection, diffusion;
  RF riso;
};




// initial conditions for component concentration
template<typename GV, typename RF>
class ConcInitial
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  ConcInitial<GV,RF> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ConcInitial<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  ConcInitial (const GV& gv, const Dune::ParameterTree& param, const std::string cname) : BaseT(gv)
  {
    if (!param.sub(cname).hasKey("initial"))
      DUNE_THROW(Dune::Exception, "There is no initial concentration for component " << cname);

    v = param.sub(cname).template get<RF>("initial");
  }

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = v;
  }
private:
  typename Traits::RangeFieldType v;
};

#endif
