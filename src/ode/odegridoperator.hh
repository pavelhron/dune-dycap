// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_ODEGRIDOPERATOR_HH
#define DUNE_DYCAP_ODEGRIDOPERATOR_HH

#include <vector>
#include<dune/common/parametertreeparser.hh>
#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include<dune/pdelab/common/elementmapper.hh>
#include<src/utilities/timer.hh>


//! Algebraic equation default class

/**
 *  this class is default class for solving of AE equations
 * and does not solve any AE
 */
class AEDefault
{
public:

  //! apply on an entity E and corresponding vector V
  template<class V, class E>
  void apply_local (V & xvector, E& e)
  {}

  //! apply solution of AE for the whole grid
  template<class V>
  void apply (std::vector<V*> & datavector)
  {}

  template<class V>
  void apply (V & v)
  {
    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(datavector);
  }

  template<class V>
  void apply (V& w, V & v)
  {
    if (&w==&v)
      DUNE_THROW(Dune::Exception, "ODEGridOperator::apply was called with 2 identical vectors.");

    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(datavector);
  }

};


/**
 * \brief ODEGridOperator, where the oxygen balance between two phases is
 *        solved in each timestep in ODE Solver
 * \tparam T          type to represent time values
 * \tparam GFS        grid function space
 * \tparam ODEMODEL   select the ODE model
 * \tparam ODESOLVER  select to ode solver which is applied for each cell
 * \tparam AEMODEL    select algebraic equation model
 * \tparam V          vector type to represent coefficients of solutions
 */
template<class T,class GV, class ODEMODEL, class ODESOLVER>
class EcoliGrowthGridOperator
{

  typedef typename ODEMODEL::number_type RF;
  typedef typename ODEMODEL::V Vector;
  enum {dim = GV::dimension};


public:
  EcoliGrowthGridOperator(const GV & gv_,
                          ODEMODEL & odemodel_,
                          ODESOLVER & odesolver_, int verbosity_=0)
    : gv(gv_), odemodel(odemodel_), odesolver(&odesolver_), verbosity(verbosity_)
  {

  }


  template<class V>
  void apply (T time, T dt, V & v)
  {
    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }

  template<class V>
  void apply (T time, T dt, V& w, V & v)
  {
    if (&w==&v)
      DUNE_THROW(Dune::Exception, "ODEGridOperator::apply was called with 2 identical vectors.");

    std::vector<V*> datavector;
    datavector.push_back(&w);
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }



  //! do one step;
  /*
   * \param[in]  time start of time step
   * \param[in]  dt suggested time step size
   * \param[out] vectors of data

   */
  template<class V>
  void apply (const T& time, const T& dt, std::vector<V*> & data)
  {
    assert(data.size()==ODEMODEL::Nr);


    // iterate over all cells
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    // map each cell to unique id
    Dune::PDELab::ElementMapper<GV> cell_mapper(gv);
    size_t counter;

    Dune::Dycap::DycapTimer watch;
    const RF TOL_old = odesolver->get_TOL();

    if (verbosity)
      std::cout << "reaction from " << time << " to " << time+dt << std::endl;
    // const RF dt_min;
    //odesolver->set_TOL(dt_min);

    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        if (it->partitionType()!=Dune::InteriorEntity)
          continue;

        watch.reset();
        counter = 0;

        // compute unique id
        typename GV::IndexSet::IndexType n = cell_mapper.map(*it);


        // setup vector with initial values
        Vector x(0);
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = (*data[i]).block(n);
            assert( i < x.size());
            x[i] = v;
          }

        if (verbosity)
          {
            for (auto it:x)
              std::cout << it << " ";
            std::cout << std::endl;
          }
        // setup problem for this cell
        odemodel.setup(*it);

        const bool adaptive = odesolver->adaptive();
        if (!adaptive)
          DUNE_THROW(Dune::Exception, "ODESolver must be with time adaptation!!");

        // set initial state
        odesolver->set_state(time, x);
        odesolver->set_dt(dt);

        while (odesolver->get_time()<time+dt-1.e-10){


          Vector y = odesolver->get_state();
          RF old_time = odesolver->get_time();
          const RF TOL =  odesolver->get_TOL();
          // advance model and save time and state
          odesolver->step();

          // last time step to T
          RF local_dt = odesolver->get_dt();

          if (!controlSolution(x) && TOL<1.e-10)
            {
              odesolver->set_TOL(TOL/static_cast<RF>(odesolver->get_order()));
              odesolver->set_state(old_time, y);
            }

          if ((odesolver->get_time()+local_dt) > time+dt)
            {
              odesolver->set_dt(time+dt-odesolver->get_time());
              if (odesolver->get_dt()<odesolver->get_dtmin())
                odesolver->set_dt(odesolver->get_dtmin()*10.);
            }

          if (verbosity)
            std::cout << "reaction timestep " << odesolver->get_dt() << " cell " << n << std::endl;

          counter++;
          if (odesolver->get_dt() < odesolver->get_dtmin())
            std::cout << "time " <<  odesolver->get_time() << " local_dt " << local_dt << " final time " << time+dt << " element " << n << std::endl;

        }

        odesolver->set_TOL(TOL_old);

        // write back the result
        x = odesolver->get_state();
        setSolution(x);
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = x[i];
            (*data[i]).block(n) = v;
          }

        if (verbosity>1)
          {
            Dune::FieldVector<RF, dim> center = it->geometry().center();
            std::cout << "element "<< n << " " << center[0] << " " << center[1] << " interationnr " <<   counter << " took " << std::setprecision(9) << std::fixed << watch.elapsed() << "s\n";
          }
      }
  }

  void setMethod(ODESOLVER& odesolver_)
  {
    odesolver = &odesolver_;
  }

  inline ODEMODEL& getODEModel()
  {
    return odemodel;
  }

private:

  inline const bool controlSolution(const Vector& x) const
  {
    for (unsigned int i=0; i<x.size(); i++)
      if (x[i]<-1.e-15)
        return false;
    return true;
  }

  inline const void setSolution(Vector& x)
  {
    for (unsigned int i=0; i<x.size(); i++)
      if (x[i]<0)
        x[i]=0;
  }

  const GV & gv;
  ODEMODEL & odemodel;
  ODESOLVER  *odesolver;
  int verbosity;
};



/**
 * \brief ODEGridOperator, where the oxygen balance between two phases is
 *        solved in each timestep in ODE Solver
 * \tparam T          type to represent time values
 * \tparam GFS        grid function space
 * \tparam ODEMODEL   select the ODE model
 * \tparam ODESOLVER  select to ode solver which is applied for each cell
 * \tparam AEMODEL    select algebraic equation model
 * \tparam V          vector type to represent coefficients of solutions
 */
template<class T,class GV, class ODEMODEL, class ODESOLVER, class AEMODEL = AEDefault>
class IterativeDAEGridOperator
{

  typedef typename ODEMODEL::number_type RF;
  typedef typename ODEMODEL::V Vector;
  enum {dim = GV::dimension};


public:
  IterativeDAEGridOperator(const GV & gv_,
                           ODEMODEL & odemodel_,
                           ODESOLVER & odesolver_)
    : gv(gv_), odemodel(odemodel_), odesolver(&odesolver_), aemodel(aeDefault()), verbosity(0)
  {
  }

  IterativeDAEGridOperator(const GV & gv_,
                           ODEMODEL & odemodel_,
                           ODESOLVER & odesolver_,
                           AEMODEL & aemodel_)
    : gv(gv_), odemodel(odemodel_), odesolver(&odesolver_), aemodel(aemodel_), verbosity(0)
  {
  }

  template<class V>
  void apply (T time, T dt, V & v)
  {
    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }

  template<class V>
  void apply (T time, T dt, V& w, V & v)
  {
    if (&w==&v)
      DUNE_THROW(Dune::Exception, "ODEGridOperator::apply was called with 2 identical vectors.");

    std::vector<V*> datavector;
    datavector.push_back(&w);
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }



  //! do one step;
  /*
   * \param[in]  time start of time step
   * \param[in]  dt suggested time step size
   * \param[out] vectors of data

   */
  template<class V>
  void apply (const T& time, const T& dt, std::vector<V*> & data)
  {
    // iterate over all cells
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    // map each cell to unique id
    Dune::PDELab::ElementMapper<GV> cell_mapper(gv);
    size_t counter;

    Dune::Dycap::DycapTimer watch;

    //std::cout << "reaction from " << time << " to " << time+dt << std::endl;


    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        watch.reset();
        counter = 0;

        // compute unique id
        typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
        // skip ghost and overlap
        if (it->partitionType()!=Dune::InteriorEntity)
          continue;


        // setup vector with initial values
        Vector x(0);
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = (*data[i]).block(n)[0];
            assert( i < x.size());
            x[i] = v;
          }

        // setup problem for this cell
        odemodel.setup(*it);



        // apply local equilibrium model
        aemodel.apply_local(x,*it);

        // set initial state
        odesolver->set_state(time, x);
        odesolver->set_dt(dt);

        // const RF dt_min;
        //odesolver->set_TOL(dt_min);

        /*
        // set dt
        bool last_dt = false;
        RF dt_local = dt;

        const bool adaptive = odesolver->adaptive();
        if (!adaptive)
        DUNE_THROW(Dune::Exception, "ODESolver must be with time adaptation!!");

        // compute ODE solution at T = time + dt
        while (odesolver->get_time() < time+dt-1.e-13)
        {
        // --- limit dt
        // don't walk further then the end time
        if (adaptive)
        dt_local = odesolver->get_dt();

        dt_local = std::min(dt_local,std::min(odemodel.suggest_dt(dt_local, x), time+dt - odesolver->get_time()));
        if (dt_local < odesolver->get_dtmin())
        dt_local = odesolver->get_dtmin()+1.e-40;

        // do dt
        Vector y = odesolver->get_state();
        RF old_time = odesolver->get_time();
        odesolver->step();
        x = odesolver->get_state();
        aemodel.apply_local(x,*it);

        if (!controlSolution(x))
        {
        dt_local = odesolver->get_dt();
        dt_local/=8.;
        odesolver->set_dt(dt_local);
        odesolver->set_state(old_time, y);
        // if (verbosity)
        std::cout << "solution control failure w dt "<< dt_local << std::endl;
        if (dt_local < odesolver->get_dtmin())
        {
        std::cout << "set solution to zero " << dt_local << std::endl;
        setSolution(x);
        }
        }

        if (counter>100000)
        std::cout << "odesolver time " << odesolver->get_time() << " timestep " << dt_local << " dt " << odesolver->get_dt() << std::endl;
        counter++;
        }
        */


        while (odesolver->get_time()<time+dt){
          // advance model and save time and state
          aemodel.apply_local(x,*it);
          odesolver->step();

          // last time step to T
          RF local_dt = odesolver->get_dt();

          if ((odesolver->get_time()+local_dt) > time+dt)
            {
              odesolver->set_dt(time+dt-odesolver->get_time());
              if (odesolver->get_dt()<odesolver->get_dtmin())
                odesolver->set_dt(odesolver->get_dtmin()*10.);
            }
          counter++;
          if (odesolver->get_dt() < odesolver->get_dtmin())
            std::cout << "time " <<  odesolver->get_time() << " local_dt " << local_dt << " final time " << time+dt << std::endl;

        }


        aemodel.apply_local(x,*it);

        // write back the result
        x = odesolver->get_state();
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = x[i];
            (*data[i]).block(n) = v;
          }

        if (verbosity)
          {
            Dune::FieldVector<RF, dim> center = it->geometry().center();
            std::cout << "element "<< n << " " << center[0] << " " << center[1] << " interationnr " <<   counter << " took " << std::setprecision(9) << std::fixed << watch.elapsed() << "s\n";
          }
      }
  }

  void setMethod(ODESOLVER& odesolver_)
  {
    odesolver = &odesolver_;
  }

  inline ODEMODEL& getODEModel()
  {
    return odemodel;
  }

private:

  inline const bool controlSolution(const Vector& x) const
  {
    for (unsigned int i=0; i<x.size(); i++)
      if (x[i]<0)
        return false;
    return true;
  }

  inline const void setSolution(Vector& x)
  {
    for (unsigned int i=0; i<x.size(); i++)
      if (x[i]<0)
        x[i]=0;
  }

  static AEMODEL & aeDefault()
  {
    static AEMODEL ae;
    return ae;
  }

  const GV & gv;
  ODEMODEL & odemodel;
  ODESOLVER  *odesolver;
  AEMODEL & aemodel;
  int verbosity;
};


/**
 * \tparam T          type to represent time values
 * \tparam GFS        grid function space
 * \tparam MODEL      select the ODE model
 * \tparam ODESOLVER  select to ode solver which is applied for each cell
 * \tparam V          vector type to represent coefficients of solutions
 */
template<class T,class GV, class MODEL, class ODESOLVER>
class ODEGridOperator
{

  typedef typename MODEL::number_type RF;
  typedef typename MODEL::V Vector;

  enum {dim = GV::dimension};

public:
  ODEGridOperator(GV & gv_,
                  MODEL & model_,
                  ODESOLVER & odesolver_)
    : gv(gv_), model(model_), odesolver(&odesolver_)
  {
  }



  template<class V>
  void apply (T time, T dt, V & v)
  {
    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }

  template<class V>
  void apply (T time, T dt, V& w, V & v)
  {
    if (&w==&v)
      DUNE_THROW(Dune::Exception, "ODEGridOperator::apply was called with 2 identical vectors.");

    std::vector<V*> datavector;
    datavector.push_back(&v);
    apply(time,dt,datavector);
  }


  //! do one step;
  /*
   * \param[in]  time start of time step
   * \param[in]  dt suggested time step size
   * \param[out] data

   */
  template<class V>
  void apply (T time, T dt, std::vector<V*> & data)
  {
    // iterate over all cells
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    // map each cell to unique id
    Dune::PDELab::ElementMapper<GV> cell_mapper(gv);

    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        if (it->partitionType()!=Dune::InteriorEntity)
          continue;

        // compute unique id
        typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
        // skip ghost and overlap
        //  if (it->partitionType()!=Dune::InteriorEntity)
        //  continue;

        // setup vector with initial values
        Vector x;
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = (*data[i]).block(n);
            assert( i < x.size() );
            x[i] = v;
          }


        // setup problem for this cell
        model.setup(*it);

        // set initial state
        odesolver->set_state(time, x);

        // set dt
        odesolver->set_dt(model.suggest_dt(dt, x));
        bool last_dt = false;

        // compute ODE solution at T = time + dt
        while (odesolver->get_time() < time+dt-1e-8)
          {
            odesolver->set_dt(model.suggest_dt(dt, x));
            // --- limit dt
            // don't walk further then the end time
            odesolver->set_dt(std::min(odesolver->get_dt(), time+dt - odesolver->get_time()));
            // avoid mini dts
            if (odesolver->get_time() + 2*odesolver->get_dt() > time+dt+1e-12 && !last_dt)
              {
                odesolver->set_dt(0.5*(time+dt - odesolver->get_time()));
                last_dt = true;
                //  std::cout << "lastdt" << std::endl;
              }
            // do dt
            odesolver->step();
            x = odesolver->get_state();
            T tt = odesolver->get_time();
            odesolver->set_state(tt, x);
            //  odesolver->get_time();
            // std::cout <<  " " << tt << std::endl;
          }

        // write back the result
        x = odesolver->get_state();
        for (unsigned int i=0; i<data.size(); i++)
          {
            RF v = x[i];
            (*data[i]).block(n) = v;
          }
      }
    // gfs.gridView().comm().barrier();

  }
private:
  GV & gv;
  MODEL & model;
  ODESOLVER* odesolver;
};


/*! \brief Phase Exchange Equilibrium \n

  \f$C_{0_2} \dots \quad\f$ dissolved oxygen concentration \n
  \f$C_{0_2,g} \dots \quad\f$ air oxygen concentration \n
  \f$C_{0_2}^* \dots \quad\f$ dissolved oxygen concentration at equilibrium    \n

  Equilibrium concentration \f$C_{0_2}^*\f$ for gases of low solubility is supplied by Henry's law in the form
  \f{equation*}{
  C_{0_2}^*= k_h \, C_{0_2,g} \,R \,T = k_H \, C_{0_2,g},
  \f}
  where \f$k_H \f$ is a Henry's law constant.

  If the phase exchange is sufficiently fast, we can use phase exchange equilibrium model in the form
  \f{equation*}{
  \begin{pmatrix}
  s_l & s_g \\
  1   & -k_H
  \end{pmatrix}
  \begin{pmatrix}
  C_{0_2}^*  \\
  C_{0_2,g}^*
  \end{pmatrix}
  =
  \begin{pmatrix}
  s_l & s_g \\
  0   & 0
  \end{pmatrix}
  \begin{pmatrix}
  C_{0_2}  \\
  C_{0_2,g}
  \end{pmatrix},
  \f}
  where \f$C_{0_2}^*\f$ and \f$C_{0_2,g}^*\f$ are the new concentrations at the equilibrium.
*/



template<class RF, class GV, class SL, class SG>
class PhaseExchangeEquilibriumGridOperator
{
  typedef RF N;
  typedef RF T;

  typedef typename GV::Grid::ctype DF;
  enum {dim = GV::dimension};

public:

  //! constructor stores parameter lambda
  PhaseExchangeEquilibriumGridOperator (Dune::ParameterTree p, GV & gv_, const SL & s_l, const SG & s_g) :
    K_H(p.get<RF>("K_H")), gv(gv_),
    s_ldgf(s_l), s_gdgf(s_g)
  {
    if (gv.comm().rank()==0)
      std::cout << "PhaseExchangeEquilibriumModel::Warning:\n " <<
        "apply method can be called only with 2 components \n" <<
        " O2water, O2 gas" << std::endl;
  }

  template<class V, class E>
  void apply_local (V & xvector, E& e)
  {


    std::vector<RF*> datavector;
    datavector.push_back(&xvector[0]);
    datavector.push_back(&xvector[1]);

    // setup vector with initial values
    for (unsigned int i=0; i<datavector.size(); i++)
      {
        RF v = (*datavector[i]);
        assert( i < x.dim() );
        x[i] = v;
      }

    // cell center
    Dune::FieldVector<RF, dim> center = e.geometry().center();
    typename SL::Traits::RangeType sl;
    typename SG::Traits::RangeType sg;
    // evaluate saturation
    s_ldgf.evaluate(e, center, sl);
    s_gdgf.evaluate(e, center, sg);

    if (sg<1.e-5)
      return;
    /*
      /           \  /     \   /           \  /         \
      | S_l   S_g |  | C_l |   | S_l   S_g |  | ~C_l |
      |           |  |     | = |           |  |      |
      |   1   -KH |  | C_g |   |   0     0 |  | ~C_g |
      \           /  \     /   \           /  \      /
    */
    A[0][0] =sl;
    A[0][1] =sg;
    A[1][0] = 0.0;
    A[1][1] = 0.0;
    b = 0.0;
    A.umv(x,b);
    A[1][0] = 1.0;
    A[1][1] = -K_H;
    A.solve(x,b);

    // write back the result
    // x = odesolver.get_state();
    for (unsigned int i=0; i<datavector.size(); i++)
      {
        RF v = x[i];
        (*datavector[i]) = v;
      }

  }

  template<class V>
  void apply (std::vector<V*> & datavector)
  {
    // iterate over all cells
    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

    Dune::PDELab::ElementMapper<GV> cell_mapper(gv);
    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {

        if (it->partitionType()!=Dune::InteriorEntity)
          continue;

        // compute unique id
        typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
        // setup vector with initial values
        for (unsigned int i=0; i<datavector.size(); i++)
          {
            RF v = (*datavector[i]).block(n);
            assert( i < x.dim() );
            x[i] = v;
          }

        // solve local problem
        {
          // cell center
          Dune::FieldVector<DF, dim> center = it->geometry().center();
          // evaluate saturation
          s_ldgf.evaluate(*it, center, sl);
          s_gdgf.evaluate(*it, center, sg);
          /*
            /           \  /     \   /           \  /      \
            | S_l   S_g |  | C_l |   | S_l   S_g |  | ~C_l |
            |           |  |     | = |           |  |      |
            |   1   -KH |  | C_g |   |   0     0 |  | ~C_g |
            \           /  \     /   \           /  \      /
          */
          A[0][0] = sl;
          A[0][1] = sg;
          A[1][0] = 0.0;
          A[1][1] = 0.0;
          b = 0.0;
          A.umv(x,b);
          A[1][0] = 1.0;
          A[1][1] = -K_H;
          A.solve(x,b);
        }

        // write back the result
        // x = odesolver.get_state();
        for (unsigned int i=0; i<datavector.size(); i++)
          {
            RF v = x[i];
            (*datavector[i]).block(n) = v;
          }
      }
  }

private:
  RF K_H;
  GV & gv;
  const SL & s_ldgf;
  const SG & s_gdgf;
  typename SL::Traits::RangeType sl;
  typename SG::Traits::RangeType sg;
  Dune::FieldMatrix<RF,2,2> A;
  Dune::FieldVector<RF,2> b;
  Dune::FieldVector<RF,2> x;
};


template<class RF, class SL, class SG>
class PhaseExchangeEquilibriumModel
{
  typedef RF N;
  typedef RF T;
  enum {dim = SL::Traits::GridViewType::dimension};
public:

  //! constructor stores parameter lambda
  PhaseExchangeEquilibriumModel (Dune::ParameterTree p, const SL & s_l, const SG & s_g, const RF exchange_lower_boundary_ = 0) :
    K_H(p.get<RF>("K_H")), s_ldgf(s_l), s_gdgf(s_g), exchange_lower_boundary(exchange_lower_boundary_)
  {
    if (s_ldgf.getGridView().comm().rank()==0)
      std::cout << "PhaseExchangeEquilibriumModel::Warning:\n " <<
        "apply method can be called only with 4 components \n" <<
        " (DOC, O2water, O2 gas, Ecoli cells" << std::endl;
  }

  template<class V, class E>
  void apply_local (V & xvector, E& e)
  {


    std::vector<RF*> datavector;
    datavector.push_back(&xvector[2]);
    datavector.push_back(&xvector[3]);

    // setup vector with initial values
    for (unsigned int i=0; i<datavector.size(); i++)
      {
        RF v = (*datavector[i]);
        assert( i < x.dim() );
        x[i] = v;
      }

    // cell center
    Dune::FieldVector<RF, dim> center = e.geometry().center();
    typename SL::Traits::RangeType sl;
    typename SG::Traits::RangeType sg;
    // evaluate saturation
    s_ldgf.evaluate(e, center, sl);
    s_gdgf.evaluate(e, center, sg);

    if (sg<exchange_lower_boundary)
      return;
    /*
      /           \  /     \   /           \  /         \
      | S_l   S_g |  | C_l |   | S_l   S_g |  | ~C_l |
      |           |  |     | = |           |  |      |
      |   1   -KH |  | C_g |   |   0     0 |  | ~C_g |
      \           /  \     /   \           /  \      /
    */
    A[0][0] =sl;
    A[0][1] =sg;
    A[1][0] = 0.0;
    A[1][1] = 0.0;
    b = 0.0;
    A.umv(x,b);
    A[1][0] = 1.0;
    A[1][1] = -K_H;
    A.solve(x,b);

    // write back the result
    for (unsigned int i=0; i<datavector.size(); i++)
      {
        RF v = x[i];
        (*datavector[i]) = v;
      }

  }

  template<class V>
  void apply (std::vector<V*> & datavector)
  {

  }

private:
  RF K_H;
  const SL & s_ldgf;
  const SG & s_gdgf;
  typename SL::Traits::RangeType sl;
  typename SG::Traits::RangeType sg;
  const RF exchange_lower_boundary;
  Dune::FieldMatrix<RF,2,2> A;
  Dune::FieldVector<RF,2> b;
  Dune::FieldVector<RF,2> x;
};


#endif
