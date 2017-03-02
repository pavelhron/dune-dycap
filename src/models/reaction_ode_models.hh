// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_REACTION_MODELS_HH
#define DUNE_DYCAP_REACTION_MODELS_HH

#include <vector>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

//==============================================================================
// reaction
//==============================================================================
//! Aerobic Growth Model for Ecoli
/*!
\tparam TP Reaction parameters
This class is used for solving aerobic growth model based on double monod kinetics
*/
template<typename TP>
class AerobicGrowthModel
{
  //! type for reaction parameters traits
  typedef typename TP::Traits Traits;
  //! type for range field
  typedef typename Traits::RangeFieldType RF;
  //! type for number type
  typedef RF N;
  //! type for time tipe
  typedef RF T;

  //! \enum number of component
  enum { Nr = 3};

public:
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
     AerobicGrowthModel (TP& tp_, int verbosityLevel_ = 0) :
    tp(tp_), verbosityLevel(verbosityLevel_)
  {
    //if (tp.getGridView().comm().rank()>0)
      verbosityLevel = 0;
  }


  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center
    Dune::FieldVector<typename Entity::ctype, Entity::dimension>
      x = e.geometry().center();
    // evaluate saturation
    mumax = tp.Mumax(e,x);
    Ks = tp.Ks(e,x);
    Ko = tp.Ko(e,x);
    Ys = tp.Ys(e,x);
    Yo = tp.Yo(e,x);
    rd = tp.Rd(e,x);
    exp = tp.Exp(e,x);

    if (verbosityLevel>0)
      {
        std::cout  << "\nmumax " << mumax
                   << "\nKs " << Ks
                  << "\nKo " << Ko
                  << "\nYs " << Ys
                  << "\nYo " << Yo
                  << "\nRd " << rd
                  << "\nexp " << exp
                  << std::endl;
      }
  }

  //! suggested time for the problem
  T suggest_dt (T dt, const V& x) const
  {
    static const T safety = 0.1;
    RF ft = mumax * std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko)*x[0];

    dt = std::min(dt, std::abs(x[0] / (ft-rd*x[0])) );
    dt = std::min(dt, std::abs(x[1] / (ft/Ys)) );
    //  dt = std::min(dt, std::abs(x[2] / (ft*Yo)) );
    if (dt<1e-3)
      std::cout << "time error " << std::endl;
    return safety * dt;
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
  }

  //! model evaluation \n
  /*!
    x[0] ... cell mass concentracion \n
    x[1] ... dissolved carbon concentration \n
    x[2] ... dissolved oxygen concentration \n
  */
  void f (const T& t, const V& x, V& result) const
  {

    RF ft = mumax * std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko);

    result[0] = (ft-rd)*x[0];
    result[1] = - (1./Ys*ft)*x[0];
    result[2] = - 0.0*Yo*ft;
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

private:
  TP& tp;
  RF mumax, Ks, Ko, Ys, Yo, rd, exp;
  int verbosityLevel;
};


//! Anaerobic Growth Model for Ecoli
/*!
\tparam TP Reaction parameters
This class is used for solving anaerobic growth model based on monod kinetics
*/
template<typename TP>
class AnaerobicGrowthModel
{
  //! type for reaction parameters traits
  typedef typename TP::Traits Traits;
  //! type for range field
  typedef typename Traits::RangeFieldType RF;
  //! type for number type
  typedef RF N;
  //! type for time tipe
  typedef RF T;

  //! \enum number of component
  enum { Nr = 2};

public:
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
  AnaerobicGrowthModel (TP& tp_, int verbosityLevel_ = 0 ) :
    tp(tp_), verbosityLevel(verbosityLevel_)
  {
    //if (tp.getGridView().comm().rank()>0)
      verbosityLevel = 0;
  }


  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center
    Dune::FieldVector<typename Entity::ctype, Entity::dimension>
      x = e.geometry().center();
    mumax = tp.Mumax(e,x);
    Ks = tp.Ks(e,x);
    Ys = tp.Ys(e,x);
    rd = tp.Rd(e,x);
    exp = tp.Exp(e,x);
    if (verbosityLevel>0)
      {
        std::cout  << "\nmumax " << mumax
                   << "\nKs " << Ks
                   << "\nYs " << Ys
                   << "\nRd " << rd
                   << "\nexp " << exp
                   << std::endl;
      }
  }

  //! suggested time for the problem
  T suggest_dt (T dt, const V& x) const
  {
    static const T safety = 0.1;
    RF ft = mumax * std::pow((x[1]/ (x[1]+Ks)),exp)*x[0];

    dt = std::min(dt, std::abs(x[0] / (ft-rd*x[0])) );
    dt = std::min(dt, std::abs(x[1] / (ft/Ys)) );
    if (dt<1e-3)
      std::cout << "time error " << std::endl;
    return safety * dt;
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
  }

  //! model evaluation \n
  /*!
    x[0] ... cell mass concentracion \n
    x[1] ... dissolved carbon concentration \n
  */
  void f (const T& t, const V& x, V& result) const
  {

    RF ft = mumax * std::pow(x[1]/ (x[1]+Ks),exp);

    result[0] = (ft-rd)*x[0];
    result[1] = - (1./Ys)*ft*x[0];
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

private:
  TP& tp;
  RF mumax, Ks, rd, Ys,exp;
  mutable RF cmax;
  int verbosityLevel;
};


/*! \brief Double Monod Model
  \tparam TP Reaction parametters\n

  \f$C_{0_2} \dots \quad\f$ dissolved oxygen concentration \n
  \f$C_{0_2,g} \dots \quad\f$ air oxygen concentration \n
  \f$C_{0_2}^* \dots \quad\f$ dissolved oxygen concentration at equilibrium    \n
  \f$k_l \, \alpha \dots \quad\f$  gas/liquid mass-transfer coefficient, oxygen mass transfer coefficient based on liquid volume and concetration driving force  \n
  \f$K_o\dots \quad\f$ oxygen half saturation constant  \n
  \f$K_s \dots \quad\f$  DOC half saturation constant \n
  \f$ C_{gluk} \dots \quad\f$  organic carbon (DOC) concentration  \n
  \f$X \dots \quad \f$  cell mass concentration \n
  \f$Y_{x,o} \dots \quad\f$ growth yield coefficient based on oxygen utilized  \n
  \f$Y_{x,g} \dots \quad\f$ growth yield coefficient based on DOC utilized  \n
  \f$\mu_{max} \dots \quad\f$ maximum specific growth rate  \n

  Equilibrium concentration \f$C_{0_2}^*\f$ for gases of low solubility is supplied by Henry's law in the form
  \f{equation*}{
  C_{0_2}^*= k_h \, C_{0_2,g} \,R \,T = k_H \, C_{0_2,g},
  \f}
  where \f$k_H \f$ is a Henry's law constant.


  The cell growth:
  \f{equation*}{
  \frac{d \, X}{d \, t} = \mu_{max}  \frac{C_{gluk}}{K_s + C_{gluk}} \frac{C_{0_2}}{K_s + C_{0_2}}X,
  \f}

  Organic carbon concentration change:
  \f{equation*}{
  \frac{d \, C_{gluk}}{d \, t} = -\frac{\mu_{max}}{Y_{x,g}}  \frac{C_{gluk}}{K_s + C_{gluk}} \frac{C_{0_2}}{K_s + C_{0_2}}X,
  \f}

  Dissolved oxygen concentration:
  \f{equation*}{
  \frac{d \, C_{0_2}}{d \, t} = k_l \left(C_{0_2}^* -C_{0_2}\right)  -\frac{\mu_{max}}{Y_{x,o}}  \frac{C_{gluk}}{K_s + C_{gluk}} \frac{C_{0_2}}{K_s + C_{0_2}}X.
  \f}

  Oxygen transfer from gas phase
  \f{equation*}{
  \frac{d \, C_{0_2,g}}{d \, t} = - k_l \left(C_{0_2}^* -C_{0_2}\right) \frac{s_l}{s_g}.
  \f}

*/
template<typename TP>
class MonodModel
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
  MonodModel (TP& tp_, int verbosityLevel_ = 0 ) :
    tp(tp_), verbosityLevel(verbosityLevel_), eps(1.e-3)
  {
    if (tp.getGridView().comm().rank()>0)
      verbosityLevel = 0;
    std::cout << "Monod Model vebosity level " << verbosityLevel << std::endl;
  }

  //! set parameters of each entity at the beginning of the computation
  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center
    Dune::FieldVector<typename Entity::ctype, Entity::dimension>
      x = e.geometry().center();
    // evaluate saturation
    sl = tp.sl(e,x);
    sg = 1.-sl;
    phi = tp.phi(e,x);
    kl = tp.k_l(e,x);
    k_H =tp.K_H(e,x);
    //aerobic model
    mumax = tp.Mumax(e,x);
    Ks = tp.Ks(e,x);
    Ko = tp.Ko(e,x);
    Ys = tp.Ys(e,x);
    Yo = tp.Yo(e,x);
    rd = tp.Rd(e,x);
    exp = tp.Exp(e,x);
    //anaerobic model
    mumax_an = tp.Mumax_an(e,x);
    Ks_an = tp.Ks_an(e,x);
    Ys_an = tp.Ys_an(e,x);
    rd_an = tp.Rd_an(e,x);
    exp_an = tp.Exp_an(e,x);
    switch_limit_min=tp.switch_limit_min(e,x);
    switch_limit_max=tp.switch_limit_max(e,x);
    exchangeswitchtype = tp.exchangeSwitchType();

    if (verbosityLevel>0)
      {
        std::cout << "sl " << sl
                  << "\nsg " << sg
                  << "\nphi " << phi
                  << "\nkl " << kl
                  << "\nk_H " << k_H
                  << "\nmumax " << mumax
                  << "\nmumax_an " << mumax_an
                  << "\nKs " << Ks
                  << "\nKs_an " << Ks_an
                  << "\nKo " << Ko
                  << "\nrd " << rd
                  << "\nrd_an " << rd_an
                  << "\nYs " << Ys
                  << "\nYs_an " << Ys_an
                  << "\nYo " << Yo
                  << "\nexponent " << exp
                  << "\nexponent_an " << exp_an
                  << "\nswitchtype " << exchangeswitchtype
                  << "\nswitchlimit_min " << switch_limit_min
                  << "\nswitchlimit_max " << switch_limit_max
                  << std::endl;
      }
  }

  //! suggested time for the problem
  T suggest_dt (T dt, const V& x) const
  {

    static const T safety = 0.2;
    const bool active_exchange = (sl<=(1.0-eps));


    RF ft = mumax * factor_function(x[2]) * std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko);
    RF ft_an = mumax_an * (1.-factor_function(x[2])) * std::pow(x[1]/ (x[1] + Ks_an),exp_an);

       //  std::cout << mumax << " " << factor_function(x[2]) << " " << std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko)  << " " <<(1.0/Yo*ft) << " " << x[0] << std::endl;

    dt = std::min(dt, 1./std::max(rd, rd_an) );
    dt = std::min(dt, std::abs(x[1]/( (1./Ys_an*ft_an + 1./Ys*ft )*x[0])) );
    dt = std::min(dt, std::abs(x[2] / ( (1.0/Yo*ft)*x[0] ) ));
    if (active_exchange)
      {
        //   dt = std::min(dt,
        //              std::abs(x[2] /  (std::min(kl*(k_H * x[3]-x[2]) - (1.0/Yo*ft)*x[0] , kl*(k_H * x[3]-x[2]))) ) );
        dt = std::min(dt,
                      std::abs(x[3] / (kl*(k_H*x[3] - x[2]) * sl / sg)) );
      }

    if (dt<1e-4){
      std::cout << active_exchange << " "  << x[2] << " " << x[3] << std::endl;
      std::cout << "time error " << dt << std::endl;
      DUNE_THROW(Dune::Exception,"timestep to small in microbiological growth");

      if (x[3] < 0.)
        std::cout << x[2] << " " <<  k_H * x[3] << " " << x[3] << std::endl;


    //    RF o2_max = 0.008009;
    //  RF o2_lim = 0.004;
    // RF o2_f = std::max(0., (o2_max-x[2])/(o2_max-o2_lim));
    //  RF o2_f = std::min(1., (x[2])/(o2_lim));

//  RF ft = mumax * factor_function(x[2]) * std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko);
//  RF ft_an = mumax_an * (1.-factor_function(x[2])) * std::pow(x[1]/ (x[1] + Ks_an),exp_an);
//  const bool active = sl<0.99;
//
//
//  dt = std::min(dt, std::abs(x[1]/(1./Ys*ft + 1./Ys_an*ft_an )));
//  dt = std::min(dt, std::abs(x[2]/(-kl*(std::max(k_H * x[3]-x[2],0.))+ (1.0/Yo*ft)*x[0])));
//
    }
    return safety * dt;
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
    x[3] ... oxygen in air concentration
  */
  void f (const T& t, const V& x, V& result) const
  {

    // RF o2_max = 0.008009;
    //RF o2_lim = 0.004;
    //  RF o2_f = std::min(1., (x[2])/(o2_lim));
    const bool active_exchange = (sl<=(1.0-eps));

    RF ft = mumax * factor_function(x[2]) * std::pow(x[1]/ (x[1] + Ks),exp) * x[2]/ (x[2] + Ko);
    RF ft_an = mumax_an * (1.-factor_function(x[2])) * std::pow(x[1]/ (x[1] + Ks_an),exp_an);

    result[0] = (ft+ft_an-std::max(rd, rd_an))*x[0];
    result[1] = - (1./Ys_an*ft_an + 1./Ys*ft )*x[0];
    result[3] = 0.;
    if (active_exchange)
      {
        result[2] = sg*kl*(k_H * x[3]-x[2]) - (1.0/Yo*ft)*x[0] ;
        result[3] = -sg*kl* sl/sg*(k_H * x[3] - x[2]);
      }
    else
      result[2] = - (1.0/Yo*ft)*x[0];

 }

  //! model derivative
  void f_x (const T& t, const V& x, Matrix& result) const
  {
    /* need only for implicit methods */
  }

  inline RF factor_function(RF o2) const
  {
    //RF sl_factor = std::pow(2.71,-std::pow(sl-0.5,2.0));
    // RF sl_factor = 1.;//std::abs(1-(sl-0.5));

    if (sg<switch_limit_min)
      return 0.;

    if (sg>switch_limit_max)
      return 1.;

    RF o2_f = (sg-switch_limit_min)/(switch_limit_max-switch_limit_min);

    if (exchangeswitchtype == "default")
      {
        if (o2 > switch_limit_min)
          return 1.;
        else return 0.;
      }

    if (exchangeswitchtype == "linear")
      return o2_f;
    if (exchangeswitchtype == "quadratic")
      return o2_f*o2_f;
    return 0.;

    /*   if (o2<switch_limit_min)
      return 0.;
    if (o2>switch_limit_max)
      return 1.*sl_factor;

    RF o2_f = (o2-switch_limit_min)/(switch_limit_max-switch_limit_min);

    if (exchangeswitchtype == "default")
      {
        if (o2 > switch_limit_min)
          return 1.*sl_factor;
        else return 0.;
      }

    if (exchangeswitchtype == "linear")
      return o2_f * sl_factor;
    if (exchangeswitchtype == "quadratic")
      return o2_f*o2_f * sl_factor;
    return 0.;
    */

    //std::min(1., (1.-sl)/(1.-Smax))*std::max(o2_f,0.);
  }

  //! change verbosity level; 0 means completely quiet
  void setVerbosityLevel (int level)
  {
    //if (s_ldgf.getGridView().comm().rank()>0)
    //   verbosityLevel = 0;
    //x else
    std::cout << "Monod Dopel Model vebostity is " << level << std::endl;
    verbosityLevel = level;
    //  std::cout << "new mmnd verbosity " << verbosityLevel << std::endl;
  }

private:
  TP& tp;
  RF sl, sg, phi, kl, k_H, mumax, Ks, Ko, Ys, Yo, rd, exp, mumax_an, Ks_an, Ys_an, rd_an, exp_an;
  RF switch_limit_min, switch_limit_max;
  std::string exchangeswitchtype;
  int verbosityLevel;
  const RF eps;
};

/*! \brief Dynamic Phase Exchange \n

  \f$C_{0_2} \dots \quad\f$ dissolved oxygen concentration \n
  \f$C_{0_2,g} \dots \quad\f$ air oxygen concentration \n
  \f$C_{0_2}^* \dots \quad\f$ dissolved oxygen concentration at equilibrium    \n
  \f$k_l \, \alpha \dots \quad\f$  gas/liquid mass-transfer coefficient, oxygen mass transfer coefficient based on liquid volume and concetration driving force  \n

  Equilibrium concentration \f$C_{0_2}^*\f$ for gases of low solubility is supplied by Henry's law in the form
  \f{equation*}{
  C_{0_2}^*= k_h \, C_{0_2,g} \,R \,T = k_H \, C_{0_2,g},
  \f}
  where \f$k_H \f$ is a Henry's law constant.


  Dissolved oxygen concentration:
  \f{equation*}{
  \frac{d \, C_{0_2}}{d \, t} = k_l \left(C_{0_2}^* -C_{0_2}\right).
  \f}

  Oxygen transfer from gas phase
  \f{equation*}{
  \frac{d \, C_{0_2,g}}{d \, t} = - k_l \left(C_{0_2}^* -C_{0_2}\right) \frac{s_l}{s_g}.
  \f}
*/

template<class TP>
class DynamicPhaseExchange
{
  typedef typename TP::Traits Traits;
  typedef typename Traits::RangeFieldType RF;
  typedef RF N;
  typedef RF T;


public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  typedef Dune::FieldVector<RF,2> V;
  typedef Dune::FieldMatrix<RF,2,2> Matrix;


  //! constructor stores parameter lambda
  DynamicPhaseExchange (TP& tp_) :
    tp(tp_) {}


  template<class Entity>
  void setup(const Entity & e)
  {
    // cell center
    Dune::FieldVector<typename Entity::ctype, Entity::dimension>
      x = e.geometry().center();
    // evaluate saturation
    sl = tp.sl(e,x);
    sg = 1.-sl;
    lambda = tp.k_l(e,x);
    K_H =tp.K_H(e,x);
  }

  T suggest_dt (T dt, const V& x) const
  {

    static const T safety = 0.1;
    dt = std::min(dt,
                  std::abs(x[0] / (lambda*(K_H*x[1] - x[0])) ) );
    dt = std::min(dt,
                  std::abs(x[1] / (lambda*(K_H*x[1] - x[0]) * sl[0] / sg[0]) ) );
    return safety * dt;
  }

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 2;
  }

  //! set initial state including time value
  void initialize (T& t0, V& x0) const
  {
    t0 = 0;
    x0[0] = 0;
    x0[1] = 0;
  }

  //! model evaluation
  /*! x[0] ... dissolved oxygen \n
    x[1] ... oxygen in the air */
  void f (const T& t, const V& x, V& result) const
  {
    result[0] = lambda*(K_H*x[1] - x[0]);
    result[1] = - result[0] * sl / sg;
  }

  //! model derivative
  void f_x (const T& t, const V& x, Matrix& result) const
  {
    // only needed for implicit methods
    result[0][0] = -lambda;
    result[0][1] = K_H*lambda;
    result[1][0] = lambda * sl / sg;
    result[1][1] = - K_H*lambda * sl / sg;
  }

private:
  TP& tp;
  RF sl;
  RF sg;
  RF K_H;
  RF lambda;
};




#endif
