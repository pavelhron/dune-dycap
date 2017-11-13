// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TRANSPORTPARAMETERS_HH
#define DUNE_DYCAP_TRANSPORTPARAMETERS_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/shared_ptr.hh>
#include<dune/common/deprecated.hh>

#include<src/physics/hydraulicparameters.hh>
#include<src/physics/physical_chemistry.hh>
#include"twophaseparameters.hh"
#include <src/utilities/inlets_utilities.hh>


namespace Dune {
  namespace PDELab {
    template<class T, class Imp>
    class ModifiedTransportSpatialParameterInterface;

    template<class T, class Imp>
    class ModifiedTransportTemporalParameterInterface;

    template<typename GV, typename RF>
    struct TransportParameterTraits;
  }
}



//==============================================================================
// Problem definition : component transport
//==============================================================================
//! Transport in water phase
template<typename GV, typename RF>
class TransportParameterBase :
  public Dune::PDELab::ModifiedTransportSpatialParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>,
                                                                  TransportParameterBase<GV,RF> >,
  public Dune::PDELab::ModifiedTransportTemporalParameterInterface<Dune::PDELab::TransportParameterTraits<GV,RF>,
                                                                   TransportParameterBase<GV,RF> >

{
  enum {dim=GV::Grid::dimension};

public:
  typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;
  typedef typename Traits::BCType BCType;


  TransportParameterBase ()
    : time(0.),tend(0.),dt(0.),diffusion_type(0), cname("not defined")
  {
  }



  /* c0 ... initial amount of oxygen
     vg  ... input
     vD  ... diffusivity */

  //! capacity function
  virtual typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own c(e,x)");
  }

  //! saturation at new time level
  virtual typename Traits::RangeFieldType
  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own snew");
  }

  //! saturation at old time level
  virtual typename Traits::RangeFieldType
  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own snew");
  }

  //! velocityvector
  virtual typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own velocity");
  }

  //! tensor permeability
  virtual typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own D");
  }

  //! source/reaction term

  virtual typename Traits::RangeFieldType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own q");
  }

  //! boundary condition type function
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  virtual BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own bctype");
  }

  //! Dirichlet boundary condition value
  virtual typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own g");
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  virtual typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    DUNE_THROW(Dune::NotImplemented, "You have to write your own j");
  }

  virtual void setRetardation(RF retardation_)
  {DUNE_THROW(Dune::NotImplemented, "setRetardation");}


  //! set time target for step in multiphase flow
  void setTimeTarget(RF time_, RF dt_)
  {
    tend = time_;
    dt = dt_;
    setTime(time_-dt);
  }

  //! to be called once before each time step
  virtual void preStep (RF time_, RF dt_, int stages)=0;

  virtual void setTime(RF time)=0;

  std::string getName() const
  {
    return cname;
  }


  void set_diffusion_type(std::string type)
  {
    if (type == "default") // both, but more for liquid
      diffusion_type = 0;
    else if (type == "Chiogna") // liquid
      {
        diffusion_type = 1;
      }
    else if (type == "MQ60") //gas
      {
        diffusion_type = 2;
      }
    else if (type == "MQ61") //gas
      {
        diffusion_type = 3;
      }
    else if (type == "EN") // both
      {
        diffusion_type = 4;
      }
    else if (type == "MMS") // both
      {
        diffusion_type = 5;
      }
    else if (type == "JJ") // both
      {
        diffusion_type = 6;
      }
    else
      DUNE_THROW(Dune::Exception, "Diffusion type " << type << " is not known.");
  }

  virtual bool computeReaction()
  {
    return true;
  }

protected:
  RF time, tend, dt;
  short diffusion_type;
  std::string cname;
};



//! Transport in water phase
template<typename GV, typename RF, typename TP, typename SOLDDGF, typename SNEWDGF, typename UDGF>
class TransportLiquid :
  public TransportParameterBase<GV,RF>
{
  enum {dim=GV::Grid::dimension};

public:
  typedef typename TransportParameterBase<GV,RF>::Traits Traits;
  typedef Dune::Dycap::TransportInlets<Traits,RF> Inlets;
  typedef Dune::LocalUniversalMapper<typename GV::Traits::Grid> IntersectionMapper;
  typedef typename GV::Grid::ctype DF;
  typedef typename TransportParameterBase<GV,RF>::BCType BCType;


  TransportLiquid (TP& tp_, const HeleshawDomain<dim>& dom_, const SOLDDGF& solddgf_, const SNEWDGF& snewdgf_, const UDGF& udgfold_, const UDGF& udgf_,  const std::string cname_)
    : tp(tp_), domain(dom_), solddgf(solddgf_), snewdgf(snewdgf_), udgfold(udgfold_), udgf(udgf_),
      // c0(tp.param.sub(cname).template get<RF>("initial")),
      Dt(tp.param.sub(cname_).template get<RF>("D")),
      grain_diameter(0.),
      csource(tp.param.sub(cname_).template get<RF>("source",0.0)),
      retardation(tp.param.sub(cname_).template get<RF>("retardation",1.0)),
      flowontop(tp.param.sub(cname_).template get<bool>("flowontop",false)),
      robin(tp.param.sub(cname_).template get<bool>("robin",false)),
      inlets(tp.getGridView(),tp.getParam(),cname_),
      velocityvector()
  {
    cname = cname_;
    this->set_diffusion_type(tp.param.sub(cname).template get<std::string>("diffusiontype"));

    if (diffusion_type == 1){
      if ( tp.param.sub(cname).template hasKey("graindiameter"))
        grain_diameter =  tp.param.sub(cname).template get<RF>("graindiameter"); // grain diameter in m
      else DUNE_THROW(Dune::Exception, "Diffusion model by Chiogna need a graindiameter parameter in .conf file ");
    }
    if (tp.getGridView().comm().rank()==0)
      std::cout << "Component " << cname  << " uses diffusion model " << diffusion_type << " with modified Robing BC " << robin << std::endl;
    createIntersectionMapper();
    velocityvector.resize(intersectionMapper->size());
    update();
  }

  //! intersection mapper
  /**
   * Returns the intersection mapper.  If there is no intersection mapper yet,
   * initialize it first.
   */
  virtual const Dune::shared_ptr<IntersectionMapper> &createIntersectionMapper() {
    if(!intersectionMapper) {
      typedef typename GV::template Codim<0>::Iterator EIterator;
      typedef typename GV::IntersectionIterator IIterator;

      intersectionMapper.reset(new IntersectionMapper(tp.getGridView().grid()));

      const EIterator &eend = tp.getGridView().template end<0>();
      for(EIterator eit = tp.getGridView().template begin<0>(); eit != eend; ++eit) {
        const IIterator &iend = tp.getGridView().iend(*eit);
        for(IIterator iit = tp.getGridView().ibegin(*eit); iit != iend; ++iit)
          // for a universal mapper, this will create a new map entry
          intersectionMapper->subIndex(*eit, iit->indexInInside(), dim-1);
      }
    }
    return intersectionMapper;
  }

  //! update the construction values for a given timestep
  void update()
  {

    Dune::Timer watch;
    watch.reset();
    // iterate over all cells
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename GV::Intersection Intersection;

    // loop over cels
    for (ElementIterator it = tp.getGridView().template begin<0>(); it!=tp.getGridView().template end<0>(); ++it)
      {
        unsigned int intersection_index = 0;
        IntersectionIterator iit = tp.getGridView().ibegin(*it);
        IntersectionIterator eiit = tp.getGridView().iend(*it);

        // loop over sides
        for(; iit!=eiit; ++iit, ++intersection_index)
          {
            // intersection
            typedef typename Dune::PDELab::IntersectionGeometry<Intersection> IG;
            IG ig(*iit,intersection_index);

            // face geometry
            const Dune::FieldVector<DF,IG::dimension-1>&
              face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

            Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
              global = ig.intersection().geometry().global(face_local);

            typename Traits::RangeType velo;
            udgf.evaluate((ig.inside()),global,velo);
            velocityvector[intersectionMapper->subIndex(ig.inside(), ig.indexInInside(), dim-1)] = velo;
          }

      }
  }

  /* c0 ... initial amount of oxygen
     vg  ... input
     vD  ... diffusivity */

  //! capacity function
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType s_new = snew(e,x);
    typename SOLDDGF::Traits::RangeType s_old = sold(e,x);
    RF factor = (tend-time)/dt;
    if (factor>(1+1.e-6) || factor<0)
      DUNE_THROW(Dune::Exception, "Bad time in component transport, time " << time << " tend " << tend << " dt " << dt);
    if (factor>1)
      factor=1.0;

    //  RF factor = 0.0;
#ifdef SCONTROL_DEBUG
    if (s_new != s_old)
      std::cout << "snew " << s_new << " sold " << s_old << " sactive "  <<(factor*s_old + (1.-factor)*s_new )
                << " tend " << tend << " time " << time << " dt " << dt << " rank "  << tp.getGridView().comm().rank()<< std::endl;
#endif
    return tp.phi(e,x)*(factor*s_old + (1.-factor)*s_new )*retardation; // * tp.nu_l(e,x,0.0);
  }

  //! saturation at new time level
  typename Traits::RangeFieldType
  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType y;
    snewdgf.evaluate(e,x,y);
    return y;
  }

  //! saturation at old time level
  typename Traits::RangeFieldType
  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SOLDDGF::Traits::RangeType y;
    solddgf.evaluate(e,x,y);
    return y;
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return v(is.inside(), is.geometryInInside().center());
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType velo;
    udgf.evaluate(e,x,velo);
    return velo;
  }

  //! tensor permeability
  /**
   * type 0   default  D*phi*S  liquid
   * type 1   Chiogna           liquid
   * type 2   MQ60     (D*phi*S)*S*phi^{1/3} gas
   * type 3   MQ61     (D*phi*S)*(S*phi)^{7/3}/(phi*phi) gas
   * type 4   EN       (D*phi*S)* tau*(S*phi)^{xl-1}/(phi)^{xl} gas
   * type 4   EN       (D*phi*S)* tau/phi liquid
   * type 5   MMS      (D*phi*S)* (S*phi)^{2*xl+1}/(phi*phi) both
   * type 6   JJ       (D*phi*S)* (S*phi)/(phi)^{2/3} both
   **/
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    RF sat = c(e,x)/tp.phi(e,x);
    RF thetaa = Dt* sat*tp.phi(e,x);

    switch (diffusion_type) {
    case 0:
      return thetaa;
      break;
    case 1:
      {
        RF velocity = v(e,x).two_norm();
        return thetaa + velocity * grain_diameter / std::sqrt(velocity*grain_diameter / Dt + 123.);
        break;
      }
    case 2:
      {
        return thetaa * sat*std::pow(tp.phi(e,x),1./3);
        break;
      }
    case 3:
      {
        return thetaa * std::pow(tp.phi(e,x)*sat,7./3)/(tp.phi(e,x)*tp.phi(e,x));
        break;
      }
    case 4:
      {
        RF tau = 0.273; // pm 0.08
        //        RF xl = 3.28; // pm 0.4

        return thetaa * tau / tp.phi(e,x); // liquid
        break;
      }
    case 5:
      {
        RF xl = 0.6*sat*sat*sat-0.758*sat*sat+0.493*sat+0.56 ;
        return thetaa * std::pow(tp.phi(e,x)*sat, 2*xl+1)/(tp.phi(e,x)*tp.phi(e,x)); // both
        break;
      }
    case 6:
      {
        return thetaa * tp.phi(e,x)*sat / std::pow(tp.phi(e,x),2./3); // both
        break;
      }

    default:
      std::cout << "diffusion model not known" << std::endl;
      break;

      return Dt*sat*tp.phi(e,x);
    }
    return -1;
  }

  //! source/reaction term
  typename Traits::RangeFieldType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {

    //    return 0.0;
    RF source=tp.q_l(e,x,time);
    if (source>0){
      return source*csource;
    }
    return source;


  }

  //! boundary condition type function
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::RangeType global = is.geometry().global(x);

    typename Traits::RangeType velo=v(is.inside(), is.geometryInInside().center());
    RF vn = velo*is.centerUnitOuterNormal();

    //outflow
    if (vn>0.)
      {
        if (robin)
          return BCType::Robin;
        return BCType::Outflow;
      }

    // no flow at the top of the domain
    if (global[dim-1]>domain.height-1e-6 && !flowontop)
      {
        if (robin)
          return BCType::Robin;
        return BCType::Neumann;
      }
    // inflow
    if (vn<0)
      {
        if (robin)
          return BCType::Robin;
      return BCType::Dirichlet;
      }

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

  // set timetarget for inlets
  void preStep (RF time_, RF dt_, int stages)
  {
    inlets.setTime(time_+dt_);
  }

  void setRetardation(RF retardation_)
  {retardation=retardation_;}


private:
  const TP& tp; // store access to two phase parameters
  const HeleshawDomain<dim>& domain;
  const SOLDDGF& solddgf;
  const SNEWDGF& snewdgf;
  const UDGF& udgfold;
  const UDGF& udgf;
  using TransportParameterBase<GV,RF>::cname;
  using TransportParameterBase<GV,RF>::time;
  using TransportParameterBase<GV,RF>::tend;
  using TransportParameterBase<GV,RF>::dt;
  using TransportParameterBase<GV,RF>::diffusion_type;
  // RF c0;
  RF Dt;
  RF grain_diameter;
  RF csource;
  RF retardation;
  const bool flowontop;
  const bool robin;
  Inlets inlets;
  std::vector<typename Traits::RangeType> velocityvector;
  Dune::shared_ptr<IntersectionMapper> intersectionMapper;
};


//! Transport in gas phase
template<typename GV, typename RF, typename TP, typename SOLDDGF, typename SNEWDGF, typename UDGF, typename POLDDGF, typename PNEWDGF>
class TransportGas :
  public TransportParameterBase<GV,RF>
{
  // typedef TwoPhaseParameter<GV,RF> TP;
  enum {dim=GV::Grid::dimension};

public:
  typedef typename TransportParameterBase<GV,RF>::Traits Traits;
  typedef typename TransportParameterBase<GV,RF>::BCType BCType;

  TransportGas (const TP& tp_, const HeleshawDomain<dim>& dom_, const SOLDDGF& solddgf_, const SNEWDGF& snewdgf_, const UDGF& udgf_, const POLDDGF& polddgf_, const PNEWDGF& pnewdgf_, const std::string cname_ = "TransportOxygenGas")
    : tp(tp_), domain(dom_), solddgf(solddgf_), snewdgf(snewdgf_), udgf(udgf_), polddgf(polddgf_), pnewdgf(pnewdgf_),
      c0(tp.param.sub(cname_).template get<RF>("initial")),
      Dt(tp.param.sub(cname_).template get<RF>("D")),
      robin(tp.param.sub(cname_).template get<RF>("robin",false))
  {
    cname = cname_;
    this->set_diffusion_type(tp.param.sub(cname).template get<std::string>("diffusiontype"));
    if (tp.getGridView().comm().rank()==0)
      std::cout << "Component " << cname  << " uses diffusion model " << diffusion_type << std::endl;
  }

  //! capacity function
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType s_new = snew(e,x);
    typename SOLDDGF::Traits::RangeType s_old = sold(e,x);
    //  typename POLDDGF::Traits::RangeType p_old;
    // typename PNEWDGF::Traits::RangeType p_new;
    //  polddgf.evaluate(e,x,p_old);
    //  pnewdgf.evaluate(e,x,p_new);
    RF factor = (tend-time)/dt;
    if (factor>(1+1.e-6) || factor<0)
      DUNE_THROW(Dune::Exception, "Bad time in component transport, time " << time << " tend " << tend << " dt " << dt);
    if (factor>1)
      factor=1.0;

    return tp.phi(e,x)*(factor*s_old + (1.-factor)*s_new ); // * tp.nu_l(e,x,0.0);
  }

  //! saturation at new time level
  typename Traits::RangeFieldType
  snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SNEWDGF::Traits::RangeType y;
    snewdgf.evaluate(e,x,y);
    return y;
  }

  //! saturation at old time level
  typename Traits::RangeFieldType
  sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename SOLDDGF::Traits::RangeType y;
    solddgf.evaluate(e,x,y);
    return y;
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return v(is.inside(), is.geometryInInside().center());
  }

  //! velocityvector
  typename Traits::RangeType
  v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::RangeType velo;
    udgf.evaluate(e,x,velo);
    return velo;
  }

  //! tensor permeability
  /**
   * type 0   default  D*phi*S  liquid
   * type 2   MQ60     (D*phi*S)*S*phi^{1/3} gas
   * type 3   MQ61     (D*phi*S)*(S*phi)^{7/3}/(phi*phi) gas
   * type 4   EN       (D*phi*S)* tau*(S*phi)^{xl-1}/(phi)^{xl} gas
   * type 5   MMS      (D*phi*S)* (S*phi)^{2*xl+1}/(phi*phi) both
   * type 6   JJ       (D*phi*S)* (S*phi)/(phi)^{2/3} both
   **/
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    RF thetaa = Dt* sold(e,x)*tp.phi(e,x);

    switch (diffusion_type) {
    case 0:
      return thetaa;
      break;
    case 2:
      {
        return thetaa * sold(e,x)*std::pow(tp.phi(e,x),1./3);
        break;
      }
    case 3:
      {
        return thetaa * std::pow(tp.phi(e,x)*sold(e,x),7./3)/(tp.phi(e,x)*tp.phi(e,x));
        break;
      }
    case 4:
      {
        RF tau = 0.273; // pm 0.08
        RF xl = 3.28; // pm 0.4

        return thetaa * tau * std::pow(tp.phi(e,x)*sold(e,x), xl-1.)/std::pow(tp.phi(e,x),xl); // gas
        break;
      }
    case 5:
      {
        RF xl = 0.6*sold(e,x)*sold(e,x)*sold(e,x)-0.758*sold(e,x)*sold(e,x)+0.493*sold(e,x)+0.56 ;
        return thetaa * std::pow(tp.phi(e,x)*sold(e,x), 2*xl+1)/(tp.phi(e,x)*tp.phi(e,x)); // gas
        break;
      }
    case 6:
      {
        return thetaa * tp.phi(e,x)*sold(e,x) / std::pow(tp.phi(e,x),2./3); // both
        break;
      }

    default:
      std::cout << "diffusion model not known" << std::endl;
      break;

      return Dt*sold(e,x)*tp.phi(e,x);
    }
    return -1;
  }

  //! source/reaction term
  typename Traits::RangeFieldType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& ) const
  {
    return 0.0;
  }

  //! boundary condition type function
  // 0 means Neumann
  // 1 means Dirichlet
  // 2 means Outflow (zero diffusive flux)
  // 4 means Modified Robin, upw c - D grad c = 0 usefull for gas transport
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    typename Traits::RangeType global = is.geometry().global(x);

    typename Traits::RangeType velo=v(is.inside(), is.geometryInInside().center());
    RF vn = velo*is.centerUnitOuterNormal();

    // if (global[1]<1.e-3)
    //    std::cout << global[0] << " " << vn << std::endl;

    //   if (std::abs(vn)<1.e-9)
    //    return 0;

    if (global[dim-1]<domain.height-1e-6)
      return BCType::Neumann;

    if (robin)
      return BCType::Robin;

     if (vn<1e-8)
      return BCType::Dirichlet;


    if (vn>1e-8)
      return BCType::Outflow;

    return BCType::Neumann;

    /*
      typename Traits::RangeType global = is.geometry().global(x);
      if (global[dim-1]>domain.height-1e-6)
      return 1;
      return 0;
    */
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    //  std::cout << vg << std::endl;
    return c0;//inlets.evaluate(is,x);
  }

  //! Neumann boundary condition
  // Good: The dependence on u allows us to implement Robin type boundary conditions.
  // Bad: This interface cannot be used for mixed finite elements where the flux is the essential b.c.
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  void setTime (RF t)
  {
    time = t;
    //inlets.setTime(t);
  }

  // set timetarget for inlets (no inlets for air)
  void preStep (RF time_, RF dt_, int stages)
  {
    //inlets.setTime(time_+dt_);
  }


private:
  const TP& tp; // store access to two phase parameters
  const HeleshawDomain<dim>& domain;
  const SOLDDGF& solddgf;
  const SNEWDGF& snewdgf;
  const UDGF& udgf;
  const POLDDGF& polddgf;
  const PNEWDGF& pnewdgf;
  using TransportParameterBase<GV,RF>::cname;
  using TransportParameterBase<GV,RF>::time;
  using TransportParameterBase<GV,RF>::tend;
  using TransportParameterBase<GV,RF>::dt;
  using TransportParameterBase<GV,RF>::diffusion_type;
  RF c0;
  RF Dt;
  const bool robin;
};



//! Transport in water phase
template<typename GV, typename RF>
class TransportLiquidSimple :
  public TransportParameterBase<GV,RF>
{
  enum {dim=GV::Grid::dimension};

public:
  typedef typename TransportParameterBase<GV,RF>::Traits Traits;
  typedef Dune::Dycap::TransportInlets<Traits,RF> Inlets;
  typedef Dune::LocalUniversalMapper<typename GV::Traits::Grid> IntersectionMapper;
  typedef typename GV::Grid::ctype DF;
  typedef typename TransportParameterBase<GV,RF>::BCType BCType;


  TransportLiquidSimple (GV& gv_, const Dune::ParameterTree & param_, const std::string cname_)
    :param(param_),
     gv(gv_),
     Dt(param.sub(cname_).template get<RF>("D")),
     grain_diameter(0.),
      csource(param.sub(cname_).template get<RF>("source",0.0)),
      flowontop(param.sub(cname_).template get<bool>("flowontop",false)),
     inlets(gv,param,cname_)

  {
    cname = cname_;
    this->set_diffusion_type(param.sub(cname).template get<std::string>("diffusiontype"));

    if (diffusion_type == 1){
      if ( param.sub(cname).template hasKey("graindiameter"))
        grain_diameter =  param.sub(cname).template get<RF>("graindiameter"); // grain diameter in m
      else DUNE_THROW(Dune::Exception, "Diffusion model by Chiogna need a graindiameter parameter in .conf file ");
    }
    if (gv.comm().rank()==0)
      std::cout << "Component " << cname  << " uses diffusion model " << diffusion_type << std::endl;
  }



  //! tensor permeability
  /**
   * type 0   default  D*phi*S  liquid
   * type 1   Chiogna           liquid
   * type 2   MQ60     (D*phi*S)*S*phi^{1/3} gas
   * type 3   MQ61     (D*phi*S)*(S*phi)^{7/3}/(phi*phi) gas
   * type 4   EN       (D*phi*S)* tau*(S*phi)^{xl-1}/(phi)^{xl} gas
   * type 4   EN       (D*phi*S)* tau/phi liquid
   * type 5   MMS      (D*phi*S)* (S*phi)^{2*xl+1}/(phi*phi) both
   * type 6   JJ       (D*phi*S)* (S*phi)/(phi)^{2/3} both
   **/
  typename Traits::RangeFieldType
  D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    RF thetaa = Dt;

    switch (diffusion_type) {
    case 0:
      return thetaa;
      break;
    default:
      std::cout << "diffusion model not known" << std::endl;
      break;
    }
    return -1;
  }

  //! source/reaction term
  typename Traits::RangeFieldType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return csource;
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

  // set timetarget for inlets
  void preStep (RF time_, RF dt_, int stages)
  {
    inlets.setTime(time_+dt_);
  }

 const Dune::ParameterTree & param;

private:
  GV& gv;
  using TransportParameterBase<GV,RF>::cname;
  using TransportParameterBase<GV,RF>::time;
  using TransportParameterBase<GV,RF>::tend;
  using TransportParameterBase<GV,RF>::dt;
  using TransportParameterBase<GV,RF>::diffusion_type;
  // RF c0;
  RF Dt;
  RF grain_diameter;
  RF csource;
  const bool flowontop;
  Inlets inlets;
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

    std::cout << "Initial concentration for " << cname << " created with initial value " << v << std::endl;
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
