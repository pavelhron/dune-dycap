// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifdef DUNE_DYCAP_REACTION_MODEL_HH
#warning ***** WARNING ***** reactionmodel.hh was already included ******
#endif

#ifndef DUNE_DYCAP_REACTION_MODEL_HH
#define DUNE_DYCAP_REACTION_MODEL_HH

#include <vector>
#include<dune/common/parametertreeparser.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>


//! traits class for reaction parameter class
template<typename GV, typename RF>
struct ReactionParameterTraits
{
  //! \brief the grid view
  typedef GV GridViewType;

  //! \brief Enum for domain dimension
  enum {
    //! \brief dimension of the domain
    dimDomain = GV::dimension
  };

  //! \brief Export type for domain field
  typedef typename GV::Grid::ctype DomainFieldType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

  //! \brief Export type for range field
  typedef RF RangeFieldType;

  //! \brief range type
  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

  //! grid types
  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
  typedef typename GV::Intersection IntersectionType;
};


template<typename GV, typename RF>
class ReactionParameter
{

public:

  //! Types related to the Traits and TwoPhaseParameter
  //! @{
  typedef ReactionParameterTraits<GV,RF> Traits;
  //! @}

  enum {dim=GV::Grid::dimension};

  //! constructor
  ReactionParameter(const GV& gv_, Dune::ParameterTree & param_) :
    gv(gv_),
    param(param_),  // parameter class
    katt(param.sub("Reaction").get<RF>("katt")),
    smax(param.sub("Reaction").get<RF>("smax")),
    beta(param.sub("Reaction").get<RF>("beta")),
    kdet(0.0),
    riso(0.0),
    phi(param.sub("Setup").get<RF>("phi")),
    bulk_density(param.sub("Setup").get<RF>("bulk")),
    grain_size(param.sub("Setup").get<RF>("grain")), // median grain size
    velocity(param.sub("Reaction").get<RF>("velocityx",0.0)),
    time(0.), t0(0), dt(1e6)

  {
  }

  void setBeta(RF beta_) {beta=beta_;}
  void setSmax(RF smax_) {smax=smax_; }
  void setKatt(RF katt_) { katt= katt_; }
  void setRiso(RF riso_){ riso= riso_;  }
  void setKdet(RF kdet_){ kdet= kdet_;  }

  void setVelocity(RF velocity_)
  {
    velocity= velocity_;
  }


  //! power parameter
  typename Traits::RangeFieldType
  getBeta(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return beta;
  }

   //! maximal adhesion
  typename Traits::RangeFieldType
  getRiso(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return riso;
  }
 //! detachment coefficient
  typename Traits::RangeFieldType
  getKdet(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return kdet;
  }

  //! adhesion parameter
  typename Traits::RangeFieldType
  getKatt(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return katt;
  }

  //! capacity parameter
  typename Traits::RangeFieldType
  getSmax(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return smax;
  }

  //! material porosity
  typename Traits::RangeFieldType
  getPhi(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return phi;
  }

  //! velocity in x direction
  typename Traits::RangeFieldType
  getVelocity(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return velocity;
  }

  //! material bulk density
  typename Traits::RangeFieldType
  getBulkDensity(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return bulk_density;
  }

    //! material grain size
  typename Traits::RangeFieldType
  getGrainSize(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return grain_size;
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



private:
  const GV& gv;
  // parameter class
  const Dune::ParameterTree & param;
  RF katt, smax, beta, kdet, riso;
  const RF phi;
  const RF bulk_density, grain_size;
  RF velocity;
  // time values
  RF time, t0, dt;

};



//! Adapter for system of equations
/*!
  \tparam M model type
  This class is used for adapting system of ODE's (right side)
  to system of PDE's (as source term)
  What need to be implemented is only the evaluate function
*/
template<class M>
class ReactionAdapter
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

  typedef typename M::V V;

  //! constructor stores reference to the model
  ReactionAdapter(M& model_)
    : model(model_), u(model.size()), result(model.size())
  {
    model.initialize(t,u);
    dt = 0.1;
  }


  //! get current time
  time_type get_time () const
  {
    return t;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_dt () const
  {
    return dt;
  }

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dt_;
  }

  //! set time step for subsequent steps
  void set_time (time_type time_)
  {
    time = time_;
  }

  //! evaluate model in entity eg with local values x
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
  {

    // model setup
    model.setup(eg);
    // copy local DOF to a vector
    for(int i = 0;i<M::Nr;i++)
      u[i]=x(lfsu,i);

    // evaluate model
    result = 0.0;
    model.f_adapter(t,u,result);   // evaluate model

    // accumulate it to the residuum
    for(int i = 0;i<M::Nr;i++){
      r.accumulate(lfsu,i,-eg.geometry().volume()*result[i]);
    }
  }

protected:
  M& model;
  time_type t, dt;
  V u;
  V result;
};




template<typename TP>
class AdhesionModel
{
  //! type for reaction parameters traits
  typedef typename TP::Traits Traits;
  //! type for range field
  typedef typename Traits::RangeFieldType RF;
  //! type for number type
  typedef RF N;
  //! type for time tipe
  typedef RF T;

public:

  //! \enum number of component
  enum { Nr = 2};

  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! type for vector of unknownts
  typedef Dune::FieldVector<RF,Nr> V;
  //! type for matrix
  typedef Dune::FieldMatrix<RF,Nr,Nr> Matrix;


  //! constructor
  /** stores reaction parameters and verbosity level
   */
  AdhesionModel (TP& tp_, int verbosityLevel_ = 0) :
    tp(tp_),
    verbosityLevel(verbosityLevel_)
  {
    verbosityLevel = 0;
  }


  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center in global coordinates
    Dune::FieldVector<typename Entity::Geometry::ctype, Entity::dimension>
      x = e.geometry().center();

    z = x[0]; // z-coordinate

    sl = 1.;
    phi = tp.getPhi(e,x);
    katt = tp.getKatt(e,x);
    smax = tp.getSmax(e,x);
    beta = tp.getBeta(e,x);
    kdet = tp.getKdet(e,x);
    riso = tp.getRiso(e,x);
    bulk_density = tp.getBulkDensity(e,x);
    grain_size = tp.getGrainSize(e,x);


  }


  //! suggested time for the problem
  T suggest_dt (T dt, const V& x)
  {
    /*
      this would be suitable only for explicit methods without timestep adaptation
      RF dtsafe = 0.5;
      dt = std::min(dt,1./(lambda));
      dt = std::min(dt,1./(lambda/(maxadhesion*1.e8)));
      return  dtsafe*dt;
    */
    return dt;
  }

  //! return number of componentes for the model
  std::size_t size () const
  {
    return Nr;
  }

  //! set initial state including time value
  void initialize (T& t0, V& x0) const
  {
    t0 = 0;
    x0[0] = 0;
    x0[1] = 0;
    assert(Nr==x0.size());
  }

  //! model evaluation \n
  /*!
    x[0] ... cell mass concentracion \n
  */
  void f (const T& t, const V& x, V& result) const
  {


    //RF d = -katt*(1-x[1]/smax)*std::pow((grain_size+z)/grain_size,-beta)*x[0];
    //std::cout << maxadhesion << std::endl;
    //  if (x[1]>smax)
    //   std::cout << smax  << " " << x[1] << std::endl;

    RF d =  -katt*std::max(1-x[1]/smax,0.0)*std::pow((grain_size+z)/grain_size,-beta)*x[0] + kdet * bulk_density/(phi*sl)*std::max(x[1],0.);
    result[0]=d/(1.+riso);
    result[1]=-(d/bulk_density*phi*sl);
  }

  void f_adapter (const T& t, const V& x, V& result) const
  {
    f(t,x,result);
    result[0]*=sl*phi;
  }

  //! model derivative
  void f_x (const T& t, const V& x, Matrix& result) const
  {
    std::cerr << "f_x not implemented!" << std::endl;
  }

  //! change verbosity level; 0 means completely quiet
  void setVerbosityLevel (int level)
  {
    verbosityLevel = level;
  }

  inline TP& getParameters()
  {
    return tp;
  }

private:
  TP& tp;
  RF sl,phi;
  RF katt, smax,beta, riso, kdet;
  RF bulk_density, grain_size, z;
  RF dtsafe;
  int verbosityLevel;
};

#endif
