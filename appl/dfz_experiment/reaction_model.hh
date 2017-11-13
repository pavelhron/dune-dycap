// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifdef DUNE_DYCAP_REACTION_MODEL_HH
#warning ***** WARNING ***** reactionmodel.hh was already included ******
#endif

#ifndef DUNE_DYCAP_REACTION_MODEL_HH
#define DUNE_DYCAP_REACTION_MODEL_HH

#include <vector>
#include<gsl/gsl_sf_gamma.h>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>


//! Ecoli Growth Model for Ecoli
/*!
  \tparam TP Reaction parameters
  This class is used for solving aerobic growth model based on double monod kinetics
*/
template<typename TP>
class EcoliGrowthAdhesionModel
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
  enum { Nr = 5};

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
  EcoliGrowthAdhesionModel (TP& tp_, int verbosityLevel_ = 0) :
    tp(tp_), param(tp.getTP().getParam()),
    growthmodel(param.sub("Microorganism").get<std::string>("model")),
    oxmumax(param.sub("Microorganism").get<RF>("mumax")/3600.),
    oxks(param.sub("Microorganism").get<RF>("Ks")),
    ko(param.sub("Microorganism").get<RF>("Ko")),
    oxys(param.sub("Microorganism").get<RF>("Ys")),
    yo(param.sub("Microorganism").get<RF>("Yo")),
    mo(param.sub("Microorganism").get<RF>("mo")/3600.),
    rd(param.sub("Microorganism").get<RF>("Rd")/3600.),
    mumax(param.sub("Microorganism").get<RF>("mumax_an")/3600.),
    ks(param.sub("Microorganism").get<RF>("Ks_an")),
    ys(param.sub("Microorganism").get<RF>("Ys_an")),
    smin(param.sub("Microorganism").get<RF>("smin")),
    xmax(param.sub("Microorganism").get<RF>("xmax")),

    exchangetype(param.sub("PhaseExchange").get<std::string>("model")),
    k_h(param.sub("PhaseExchange").get<RF>("K_H")),
    md(param.sub("PhaseExchange").get<RF>("md")),
    kappa(param.sub("PhaseExchange").get<RF>("kappa")),

    katt(param.sub("Adhesion").get<RF>("katt")),
    maxadhesion(param.sub("Adhesion").get<RF>("maxadhesion")*param.sub("Microorganism").get<RF>("weight")),

    Dw(param.sub("TransportOxygenWater").get<RF>("D")),
    // alpha = alpha' / 100 * rho * g
    alpha(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("alpha")),
    n(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("n")),
    m(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("m", 1.-1./n)),
    sgmin(param.sub("Setup").get<RF>("sgmin")),

    kdet(param.sub("Adhesion").get<RF>("kdet")),
    velocity_dependence(param.sub("Adhesion").get<bool>("velocitydependence")),
    sat_correction(param.sub("Adhesion").get<bool>("saturation")),
    verbosityLevel(verbosityLevel_)
  {
    if (tp.getGridView().comm().rank()>0)
      verbosityLevel = 0;

    sl=sg=phi=0;
  }


  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    typedef typename Traits::DomainFieldType ct;
    // cell center local
    auto x  = Dune::ReferenceElements<ct, Entity::dimension>::general(e.geometry().type()).position(0,0);

    sl = tp.slnew(e,x);
    sg = 1.-sl;
    phi = tp.phi(e,x);
    RF velocity = tp.velocity(e,x);
    RF sigma = tp.getTP().sigma_l(e,x);
    RF pc = tp.getTP().pc(sl);

    beta = Dw*(2/md + sqrt(velocity/(md*M_PI*Dw)));

    if (exchangetype == "Gvirtzmann")
      agw = 6*(1-phi)/md*sg*kappa;
    else if (exchangetype == "Niemet")
      {
        const RF w = m-1./n;
        const RF z = 1.+1./n;
        const RF u = std::pow(sl,1/m);
        agw = 3*phi/(2.*sigma)*m/alpha *gsl_sf_beta(w,z)*(1.-gsl_sf_beta_inc(w,z,u));
      }
    else if (exchangetype == "Miller")
      agw = pc*phi/sigma*sg;
    else if (exchangetype == "quadratic")
      agw = (1699.4 - 1969.8*sl - 1.4*pc +274*sl*sl + 1.219*sl*pc +0.000289799*pc*pc);
    else if (exchangetype =="Geistlinger")
      agw = 6.*phi/md*sg*kappa;
    else
      DUNE_THROW(Dune::Exception, "exchange model " << exchangetype << " is not known.");



    kattv = katt;
    kdetv = kdet;
    if (velocity_dependence)
      {
        kattv*=velocity;
        kdetv*=velocity;
      }


    if (verbosityLevel>2)
      {
        std::cout  << "\noxmumax " << oxmumax
                   << "\noxks " << oxks
                   << "\nko " << ko
                   << "\noxys " << oxys
                   << "\nyo " << yo
                   << "\nalpha " << alpha
                   << "\nmumax " << mumax
                   << "\nsmin " << smin
                   << "\nks " << ks
                   << "\nys " << ys
                   << "\nRd " << rd
                   <<"\nsl " << sl
                   <<"\nsg " << sg
                   << std::endl;
      }
  }

  //! suggested time for the problem
  T suggest_dt (T dt, const V& x)
  {
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
    x0[2] = 0;
    x0[3] = 0;
    x0[4] = 0;
  }

  //! model evaluation \n
  /*!
    x[0] ... cell mass concentracion \n
    x[1] ... dissolved carbon concentration \n
    x[2] ... dissolved oxygen concentration \n
    x[3] ... oxygen concentration in air \n
  */
  void f (const T& t, const V& x, V& result) const
  {
    for (auto &it:result)
      it=0;

    f_adapter(t,x,result);
    result[0] /= (sl*phi);
    result[1] /= (sl*phi);
    result[2] /= (sl*phi);
    result[3] /= (sg*phi);
    result[4] /= (sl*phi);

    for (auto it=result.begin();it!=result.end();++it)
      if (!std::isfinite(*it)) std::cout << "element "<< it - result.begin()<< "  is infinite" << std::endl;
  }

  void f_adapter (const T& t, const V& x, V& result) const
  {


    Dune::FieldMatrix<RF,2,2> A;
    Dune::FieldVector<RF,2> b;
    Dune::FieldVector<RF,2> xx(0.);

    // const RF Mn = 28; // g/mol
    const RF Mo = 32.; // g/mol
    //  const RF Mc = 44; // g/mol

    const RF oxygen_air = x[3]*Mo; // [g/m3 air]=[mg/dm3]

    const RF oxygen_water = std::max(x[2],0.0); // [mg/dm3] water
    //  const RF o2weight_water = oxygen_water * sl * phi; //mg
    // const RF o2weight_air = oxygen_air * sg * phi; //mg

    // compute oxygen in water under equilibrium
    RF oxygen_equilibrium =  oxygen_water;
    // skip cells with small gas content
    if (sg>sgmin)
      {
        xx[0] =  oxygen_water;
        xx[1] = oxygen_air;

        A[0][0] =sl * phi;
        A[0][1] =sg * phi;

        A[1][0] = 0.0;
        A[1][1] = 0.0;
        b = 0.0;
        A.umv(xx,b);
        A[1][0] = 1.0;
        A[1][1] = -k_h;
        A.solve(xx,b);
        oxygen_equilibrium = xx[0];
      }

    RF cells = x[0]+x[4];
    RF oxft=oxmumax * S(x[1])/ (S(x[1])+oxks*cells+eps)*(oxygen_water/(ko*cells+oxygen_water+eps));
    RF ft=mumax * S(x[1])/ (S(x[1])+ks*cells+eps);

    // maximal growth rate
    ft = std::max(ft-oxft,0.0);

    // this makes the maximal adhesion proportional to water saturation
    RF alimit = std::max((1-x[4]/maxadhesion),0.0);
    // maximal adhesion is anti-proportional to water content, it does not make sence
    if (sat_correction)
      alimit = (1-x[4]*(2.-sl)/maxadhesion);


    RF adhesion = kattv * std::max(x[0],0.0)*sl*phi*alimit - sl*phi*kdetv*std::max(x[4],0.0);

    /*
      this makes difference in inflow in lower part of the domain
    if (!adh_correction)
      result[0] = (oxft + ft)*cells*sl*phi - rd* x[0]*sl*phi;// -adhesion;
    else
    */
    result[0] = (oxft + ft)*cells*sl*phi - rd* x[0]*sl*phi -adhesion;
    result[1] = - ((oxft/oxys) + (ft/ys))*cells*sl*phi;
    result[4] =  adhesion - rd*x[4]*sl*phi;

    // new version
    RF uptake = beta * agw * (oxygen_equilibrium - oxygen_water);
    RF uptake_water = - (1./yo)*oxft*cells*1000.*sl*phi+ uptake;
    RF uptake_air = - uptake;

    uptake_water-=mo*cells*1000.*oxygen_water*sl*phi;

    /*
    if (verbosityLevel>1){
      std::cout << "total oxygen uptake " << uptake << " mg"  << std::endl;
      std::cout << "oxygen uptake in water " << uptake_water << " [mg/dm3 water]"  << std::endl;
      std::cout << "oxygen uptake in air   " << uptake_air << " [mg/dm3 air]"  << std::endl;
    }
    */

    result[2]=uptake_water;
    result[3] = uptake_air/Mo;

  }

  inline void growth_function(const RF & x0, const RF & x1,const RF & oxygen_water, RF &ft, RF &oxft) const
  {
    if (growthmodel == "monod")
      {
        ft = mumax * S(x1)/ (S(x1)+ks);
        oxft = oxmumax * S(x1)/ (S(x1)+oxks) * oxygen_water/(ko+oxygen_water);
      }
    else if (growthmodel == "tessier")
      {
        ft = mumax * (1-std::exp(-S(x1)/ks));
        oxft = oxmumax *  (1-std::exp(-S(x1)/oxks)) * (1-std::exp(-oxygen_water/ko));
      }
    else if (growthmodel == "contois")
      {
        ft = mumax * S(x1)/ (S(x1)+ks*x0+eps);
        oxft = oxmumax * S(x1)/ (S(x1)+oxks*x0+eps)*(oxygen_water/(ko*x0+oxygen_water+eps));

      }
    else
      DUNE_THROW(Dune::Exception, "bacterial model " << growthmodel << " is not known.");
  }

  inline  RF S(RF s) const
  {
    return std::max(s-smin,0.0);
  }


  //! model derivative
  void f_x (const T& t, const V& x, Matrix& result) const
  {
    /* need only for implicit methods */
  }

  //! change verbosity level; 0 means completely quiet
  void setVerbosityLevel (int level)
  {
    verbosityLevel = level;
  }

  void setExchangeType (std::string type_)
  {
    tp.setExchangeType(type_);
  }

private:
  TP& tp;
  const Dune::ParameterTree & param;
  const RF eps=1.e-14;
  RF sl,sg,phi;
  const std::string growthmodel;
  const RF oxmumax, oxks, ko, oxys, yo, mo, rd;
  const RF mumax, ks, ys, smin,xmax;
  const std::string exchangetype;
  const RF k_h, md, kappa;
  const RF katt, maxadhesion;
  const RF Dw;
  const RF alpha,n,m;
  const RF sgmin;
  RF kdet, beta, agw, kattv, kdetv;
  const bool velocity_dependence;
  const bool sat_correction;
  RF dtsafe;
  int verbosityLevel;
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


//! Ecoli Growth Model for Ecoli
/*!
  \tparam TP Reaction parameters
  This class is used for solving aerobic growth model based on double monod kinetics
*/
template<typename TP>
class EcoliGrowthModel
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
  enum { Nr = 4};

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
  EcoliGrowthModel (TP& tp_, int verbosityLevel_ = 0) :
    tp(tp_), param(tp.getTP().getParam()),
    growthmodel(param.sub("Microorganism").get<std::string>("model")),
    oxmumax(param.sub("Microorganism").get<RF>("mumax")/3600.),
    oxks(param.sub("Microorganism").get<RF>("Ks")),
    ko(param.sub("Microorganism").get<RF>("Ko")),
    oxys(param.sub("Microorganism").get<RF>("Ys")),
    yo(param.sub("Microorganism").get<RF>("Yo")),
    mo(param.sub("Microorganism").get<RF>("mo")/3600.),
    rd(param.sub("Microorganism").get<RF>("Rd")/3600.),
    mumax(param.sub("Microorganism").get<RF>("mumax_an")/3600.),
    ks(param.sub("Microorganism").get<RF>("Ks_an")),
    ys(param.sub("Microorganism").get<RF>("Ys_an")),
    smin(param.sub("Microorganism").get<RF>("smin")),
    xmax(param.sub("Microorganism").get<RF>("xmax",1.e9)),

    exchangetype(param.sub("PhaseExchange").get<std::string>("model")),
    k_h(param.sub("PhaseExchange").get<RF>("K_H")),
    md(param.sub("PhaseExchange").get<RF>("md")),
    kappa(param.sub("PhaseExchange").get<RF>("kappa")),

    Dw(param.sub("TransportOxygenWater").get<RF>("D")),
    // alpha = alpha' / 100 * rho * g
    alpha(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("alpha")),
    n(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("n")),
    m(param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("m", 1.-1./n)),
    sgmin(param.sub("Setup").get<RF>("sgmin")),
    verbosityLevel(verbosityLevel_)
  {
    if (tp.getGridView().comm().rank()>0)
      verbosityLevel = 0;

    sl=sg=phi=0;
  }


  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center
    Dune::FieldVector<typename Entity::Geometry::ctype, Entity::dimension>
      x = e.geometry().center();

    if (verbosityLevel)
      std::cout << "\ncoordinates "<<x[0] << " " << x[1] << std::endl;

    sl = tp.slnew(e,x);
    sg = 1.-sl;
    phi = tp.phi(e,x);
    RF sigma = tp.getTP().sigma_l(e,x);
    RF pc = tp.getTP().pc(sl);

    beta = 2*Dw/md;

    if (exchangetype == "Gvirtzmann")
      agw = 6.*(1.-phi)/md*sg*kappa;
    else if (exchangetype == "Niemet")
      {
        const RF w = m-1./n;
        const RF z = 1.+1./n;
        const RF u = std::pow(sl,1/m);
        agw = 3*phi/(2.*sigma)*m/alpha *gsl_sf_beta(w,z)*(1.-gsl_sf_beta_inc(w,z,u));
      }
    else if (exchangetype == "Miller")
      agw = pc*phi/sigma*sg;
    else if (exchangetype == "quadratic")
      agw = (1699.4 - 1969.8*sl - 1.4*pc +274*sl*sl + 1.219*sl*pc +0.000289799*pc*pc);
    else if (exchangetype =="Geistlinger")
      agw = 6.*phi/md*sg*kappa;
    else
      DUNE_THROW(Dune::Exception, "exchange model " << exchangetype << " is not known.");

    agw=std::max(0.,agw);

    if (verbosityLevel>2)
      {
        std::cout  << "\noxmumax " << oxmumax
                   << "\noxks " << oxks
                   << "\nko " << ko
                   << "\noxys " << oxys
                   << "\nyo " << yo
                   << "\nalpha " << alpha
                   << "\nmumax " << mumax
                   << "\nsmin " << smin
                   << "\nks " << ks
                   << "\nys " << ys
                   << "\nRd " << rd
                   <<"\nsl " << sl
                   <<"\nsg " << sg
                   << std::endl;
      }
  }

  //! suggested time for the problem
  T suggest_dt (T dt, const V& x)
  {
    DUNE_THROW(Dune::Exception,"function suggest_dt for ODE model was not implemented");
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
    x0[2] = 0;
    x0[3] = 0;
  }

  //! model evaluation \n
  /*!
    x[0] ... cell mass concentracion \n
    x[1] ... dissolved carbon concentration \n
    x[2] ... dissolved oxygen concentration \n
    x[3] ... oxygen concentration in air \n
  */

  void f (const T& t, const V& x, V& result) const
  {
    result[0] = 0;
    result[1] = 0;
    result[2] = 0;
    result[3] = 0;

    /*  for (auto it:x)
      if (it<0){
        for (auto itt:x)
          std::cout << itt << " ";
        std::cout << " time " << t << std::endl;
        break;
      }
    */
    f_adapter(t,x,result);
    result[0] /= (sl*phi);
    result[1] /= (sl*phi);
    result[2] /= (sl*phi);
    result[3] /= (sg*phi);
    if (verbosityLevel)
      std::cout << "f finished" << std::endl;

    if (!std::isfinite(result[0])) std::cout << "0 is infinite" << std::endl;
    if (!std::isfinite(result[1])) std::cout << "1 is infinite" << std::endl;
    if (!std::isfinite(result[2])) std::cout << "2 is infinite" << std::endl;
    if (!std::isfinite(result[3])) std::cout << "3 is infinite" << std::endl;
  }

  void f_adapter (const T& t, const V& x, V& result) const
  {

    Dune::FieldMatrix<RF,2,2> A;
    Dune::FieldVector<RF,2> b;
    Dune::FieldVector<RF,2> xx(0.);

    if (verbosityLevel>1)
      {
        std::cout.precision(5);
        std::cout << "ecoli             " << x[0] << " [g/dm3] water"
                  << "\ndoc             " << x[1] << " [g/dm3] water"
                  << "\noxygen in water " << x[2] << " [mg/dm3] water"
                  << "\noxygen in air   " << x[3] << " [mol/m3] air" << std::endl;
      }
    const RF Mo = 32; // g/mol

    const RF oxygen_air = x[3]*Mo; // [g/m3 air]=[mg/dm3]
    RF oxygen_water = x[2];//std::max(x[2],0.0); // [mg/dm3] water


    const RF o2weight_water = oxygen_water * sl * phi; //mg
    const RF o2weight_air = oxygen_air * sg * phi; //mg

    if (verbosityLevel>1){
      std::cout << std::fixed << "overall weight of oxygen:"
                <<"\nwater " << o2weight_water<< " mg"
                <<"\nair   " << o2weight_air<< " mg"
                << "\noverall " << o2weight_water + o2weight_air << " mg" <<std::endl;
    }

    /*   if (x[2]<0)
      {
        std::cout << "oxygen " << x[2] << " " << x[3]/Mo  << " overall " <<  o2weight_water + o2weight_air << " mg" << std::endl;
        oxygen_water=0.0;
      }
    */
     RF oxygen_equilibrium =  oxygen_water;

    if (sg>sgmin)
      {
        xx[0] =  oxygen_water;
        xx[1] = oxygen_air;

        A[0][0] =sl * phi;
        A[0][1] =sg * phi;

        A[1][0] = 0.0;
        A[1][1] = 0.0;
        b = 0.0;
        A.umv(xx,b);
        A[1][0] = 1.0;
        A[1][1] = -k_h;
        A.solve(xx,b);
        oxygen_equilibrium = xx[0];
      }

    if (verbosityLevel>1){
      std::cout << std::setw(8) << "oxygen concentration in momentan equilibrium:"
                << "\nwater " << xx[0] << " mg/dm3"
                << "\nair " << xx[1] << " mg/dm3\n";
      std::cout << std::setw(8) << "oxygen concentration in equilibrium:"
                << "\nwater " << 21*10.*1.2*k_h << " mg/dm3"
                << "\nair " << 21/100.*1.2*1000 << " mg/dm3\n\n";
    }

    RF oxft;
    RF ft;
    growth_function(x[0],x[1],oxygen_water,ft,oxft);

    ft = std::max(ft-oxft,0.0);
    result[0] = (oxft + ft -rd)*x[0]*sl*phi;
    result[1] =( - (1./oxys)*oxft*x[0] - (1./ys)*ft*x[0])*sl*phi;

    if (verbosityLevel>1)
      std::cout << "ft " << ft << " oxft " << oxft << " 1/oxys " << (1./oxys)  << " " << (1./oxys)*oxft*x[0] << std::endl;

    // new version
    RF uptake = beta * agw * (oxygen_equilibrium - oxygen_water);
    // smaller interface for small water saturation
    RF uptake_water = - (1./yo)*oxft*x[0]*1000.*sl*phi+ uptake;
    RF uptake_air = - uptake;

    if (verbosityLevel)
      std::cout << "mo " << mo << " uptake " << uptake_water << " mo " <<  mo*x[0]*1000.*oxygen_water*sl*phi << std::endl;

    uptake_water-=mo*x[0]*1000.*oxygen_water*sl*phi;

    if (verbosityLevel>1){
      std::cout << "total oxygen uptake " << uptake << " mg"  << std::endl;
      std::cout << "oxygen uptake in water " << uptake_water << " [mg/dm3 water]"  << std::endl;
      std::cout << "oxygen uptake in air   " << uptake_air << " [mg/dm3 air]"  << std::endl;
    }

    result[2] = uptake_water;
    result[3] = uptake_air/Mo;

    if (verbosityLevel)
      std::cout << "old " << x[2] << " new " << result[2]*dtsafe + x[2] << std::endl;

  }

  inline void growth_function(const RF & x0, const RF & x1,const RF & oxygen_water, RF &ft, RF &oxft) const
  {
    if (growthmodel == "monod")
      {
        ft = mumax * S(x1)/ (S(x1)+ks);
        oxft = oxmumax * S(x1)/ (S(x1)+oxks+1.e-40) * oxygen_water/(ko+oxygen_water+1.e-40);
      }
    else if (growthmodel == "tessier")
      {
        ft = mumax * (1-std::exp(-S(x1)/ks));
        oxft = oxmumax *  (1-std::exp(-S(x1)/oxks)) * (1-std::exp(-oxygen_water/ko));
      }
    else if (growthmodel == "contois")
      {
        ft = mumax * S(x1)/ (S(x1)+ks*x0+1.e-40);
        oxft = oxmumax * S(x1)/ (S(x1)+oxks*x0+1.e-40)*(oxygen_water/(ko*x0+oxygen_water+1.e-40));
      }
    else
      DUNE_THROW(Dune::Exception, "bacterial model " << growthmodel << " is not known.");

  }

  inline  RF S(RF s) const
  {
    return std::max(s-smin,0.0);
  }


  //! model derivative
  void f_x (const T& t, const V& x, Matrix& result) const
  {
    /* need only for implicit methods */
  }

  //! change verbosity level; 0 means completely quiet
  void setVerbosityLevel (int level)
  {
    verbosityLevel = level;
  }



  template<class Entity>
  inline RF getSl(const Entity & e) const
  {
    return sl;
  }

  void setExchangeType (std::string type_)
  {
    tp.setExchangeType(type_);
  }

private:
  TP& tp;
  const Dune::ParameterTree & param;

  RF sl,sg,phi;
  const std::string growthmodel;
  const RF oxmumax, oxks, ko, oxys, yo, mo, rd;
  const RF mumax, ks, ys, smin,xmax;
  const std::string exchangetype;
  const RF k_h, md, kappa;
  const RF Dw;
  const RF alpha,n,m;
  const RF sgmin;
  RF beta, agw;
  RF dtsafe;
  int verbosityLevel;
};


#endif
