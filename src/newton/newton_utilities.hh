// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_NEWTON_UTILS
#define DUNE_DYCAP_NEWTON_UTILS

#include<dune/common/parametertree.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/grid.hh>
#include <dune/geometry/referenceelements.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/backend/istl.hh>

namespace Dune {
  namespace PDELab {

    //! interpret a parameter tree as a set of options for the newton solver
    /**

       example configuration:

       \code
       [NewtonParameters]

       ReassembleThreshold = 0.1
       LineSearchMaxIterations = 10
       MaxIterations = 7
       AbsoluteLimit = 1e-6
       Reduction = 1e-4
       LinearReduction = 1e-3
       LineSearchDamping  = 0.9
       Verbosity = 2
       \endcode
    */
    struct NewtonParameters : public Dune::ParameterTree
    {
    public:
      NewtonParameters(const Dune::ParameterTree p) :
        Dune::ParameterTree(p)
      {
      }
      NewtonParameters & operator = (const Dune::ParameterTree & p)
      {
        static_cast<ParameterTree &>(*this) = p;
        return *this;
      }

      //! apply the parameter set to an instance of Dune::PDELab::Newton
      template<typename N>
      void set(N & newton) const
      {
        if (this->hasKey("ReassembleThreshold"))
          newton.setReassembleThreshold(
                                        this->get<double>("ReassembleThreshold"));
        if (this->hasKey("LineSearchMaxIterations"))
          newton.setLineSearchMaxIterations(
                                            this->get<int>("LineSearchMaxIterations"));
        if (this->hasKey("MaxIterations"))
          newton.setMaxIterations(
                                  this->get<int>("MaxIterations"));

        if (this->hasKey("ForceIteration"))
          newton.setForceIteration(
                                   this->get<bool>("ForceIteration"));

        if (this->hasKey("LineSearchStrategy"))
          newton.setLineSearchStrategy(
                                       this->get<std::string>("LineSearchStrategy"));

        if (this->hasKey("AbsoluteLimit"))
          newton.setAbsoluteLimit(
                                  this->get<double>("AbsoluteLimit"));
        if (this->hasKey("Reduction"))
          newton.setReduction(
                              this->get<double>("Reduction"));
        if (this->hasKey("LinearReduction"))
          newton.setMinLinearReduction(
                                       this->get<double>("LinearReduction"));
        if (this->hasKey("LineSearchDamping"))
          newton.setLineSearchDampingFactor(
                                            this->get<double>("LineSearchDamping"));
        if (this->hasKey("Verbosity"))
          newton.setVerbosityLevel(
                                   this->get<int>("Verbosity"));
      }
    };

  }
}

template<typename V>
class NewtonHelperBase
{
public:
  virtual void set_newton_iteration(int iteration_){}
  virtual void ls_iteration(int ls, V &v){}
  virtual void residual_change(V& r){}
};


template<class Writer, class T, class Residuum, class V, class ResiduumStationary = Residuum, class VS = V>
class NewtonHelper : public NewtonHelperBase<V>
{
public:

  typedef typename T::Traits::RangeFieldType RF;

  NewtonHelper(Writer & writer_, T & t_, Residuum & residuum_, V & resvector_)
    : writer(writer_), t(t_), residuum(residuum_), resvector(resvector_), residuum_stationary(residuum_), resvector_stationary(resvector_), mass_computed(false), mass_liquid(0.), mass_gas(0.), iteration(0)
  {
  }

  NewtonHelper(Writer & writer_, T & t_, Residuum & residuum_, V & resvector_, ResiduumStationary & residuum_stationary_, VS & resvector_stationary_)
    : writer(writer_), t(t_), residuum(residuum_), resvector(resvector_), residuum_stationary(residuum_stationary_), resvector_stationary(resvector_stationary_),  mass_computed(false), mass_liquid(0.), mass_gas(0.), iteration(0)
  {
  }

  //! set new newton iteration
  void set_newton_iteration(int iteration_)
  {
    iteration = -iteration_;
    if (iteration_>=0)
      writer.set_iteration(iteration_);
    mass_computed = false;
  }

  inline int get_newton_iteration()
  {
    return iteration;
  }

  //! ls linesearch iteration
  //! v solution vector
  template<class VV>
  void ls_iteration(int ls, VV &v)
  {

    resvector = residuum.defect(v);
    resvector_stationary = residuum_stationary.defect(v);
    writer.write(ls);
  }


  void residual_change(V& r)
  {
    if (!t.compute_correction()) return;

    if (!mass_computed)
      compute_mass();

    typedef typename T::Traits::GridViewType::Traits::template Codim<0>::Iterator ElementIterator;
    typedef typename T::Traits::GridViewType::Traits::template Codim<0>::Entity Element;
    typedef typename T::Traits::GridViewType GV;


    //! \brief export type for range field
    // typedef typename  T::Traits::RangeFieldType RF;
    enum { dim =  T::Traits::dimDomain };

    typedef typename Dune::PDELab::ElementGeometry<Element> EG;

    Dune::PDELab::ElementMapper<GV> cell_mapper(t.getGridView());
    // loop once over the grid
    for (ElementIterator it = t.getGridView().template begin<0>();
         it!=t.getGridView().template end<0>(); ++it)
      {
        typename GV::IndexSet::IndexType id = cell_mapper.map(*it);
        EG eg(*it);
        // cell geometry
        //const Dune::FieldVector<RF,dim>&
        //  cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

        typename V::ElementType cont;
        //cont[0]= mass_liquid;
        // cont[1]=mass_gas;
        //  back[0] = mass_liquid;
        // back[1] = mass_gas;

        r.block(0) *= mass_liquid;//t.get_mass_l(eg.entity(),cell_center_local);
        r.block(1) *= mass_gas;//t.get_mass_g(eg.entity(),cell_center_local);
        /*  r.block(2) *= mass_liquid;
        r.block(3) *= mass_gas;
        */
      }

  }


  RF getMassL()
  {
    return mass_liquid;
  }

  RF getMassG()
  {
    return mass_gas;
  }

  void compute_mass()
  {
    RF mass_liquid_local = 0.0;
    RF mass_gas_local = 0.0;

    typedef typename T::Traits::GridViewType::Traits::template Codim<0>::Iterator ElementIterator;
    typedef typename T::Traits::GridViewType::Traits::template Codim<0>::Entity Element;

    //! \brief export type for range field
    typedef typename  T::Traits::RangeFieldType RF;
    enum { dim =  T::Traits::dimDomain };

    typedef typename Dune::PDELab::ElementGeometry<Element> EG;

    // loop once over the grid
    for (ElementIterator it = t.getGridView().template begin<0>();
         it!=t.getGridView().template end<0>(); ++it)
      {
        EG eg(*it);
        // cell geometry
        const Dune::FieldVector<RF,dim>&
          cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

        // skip ghost and overlap
        if (it->partitionType()!=Dune::InteriorEntity)
          continue;
        mass_liquid_local+=t.get_mass_l(eg.entity(),cell_center_local);
        mass_gas_local+=t.get_mass_g(eg.entity(),cell_center_local);
      }
    mass_liquid = t.getGridView().comm().sum(mass_liquid_local);
    mass_gas = t.getGridView().comm().sum(mass_gas_local);

    // if (t.getGridView().comm().rank() == 0)
    std::cout << "mass for reduction, liquid: " << mass_liquid << " gas: " << mass_gas << std::endl;

    mass_computed = true;
  }

private:

  Writer & writer;
  T & t;
  Residuum & residuum;
  V& resvector;
  ResiduumStationary & residuum_stationary;
  VS& resvector_stationary;
  bool mass_computed;
  RF mass_liquid;
  RF mass_gas;
  int iteration;

};

/**
 * \brief mass correction which can be used for twophase problem
 *
 * \tparam T     twophase traits to represent twophase problem
 * \tparam PL    liquid pressure dgf
 * \tparam PG    gass pressure dgf
 */
template<typename  T, typename PL, typename PG>
class NewtonMassCorrection

{
  const T& t;
  const PL& pl;
  const PG& pg;
  const bool correction;

public:

  typedef typename T::Traits Traits;
  typedef typename Traits::RangeFieldType RF;

  NewtonMassCorrection (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_), correction((t.getParam().sub("Setup")).get("newtoncorrection",false))
  {}

  //! get liquid mass on the element
  typename Traits::RangeFieldType
  get_mass_l(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (!correction)
      return 1.;
    typename PL::Traits::RangeType pl_value, pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    RF cell_volume = e.geometry().volume();
    return t.phi(e,x)*t.rho_l(e,x,pg_value)*(t.s_l(pg_value-pl_value))*cell_volume;
  }

  //! get gas mass on the element
  typename Traits::RangeFieldType
  get_mass_g(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    if (!correction)
      return 1.;
    typename PL::Traits::RangeType pl_value, pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    RF cell_volume = e.geometry().volume();
    return t.phi(e,x)*t.rho_g(e,x,pg_value)*(1.-t.s_l(pg_value-pl_value))*cell_volume;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }

  inline const bool compute_correction()
  {
    return correction;
  }

};

#endif // DUNE_DYCAP_NEWTON_UTILS
