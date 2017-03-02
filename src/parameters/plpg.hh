// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PLPG_HH
#define DUNE_DYCAP_PLPG_HH

#include <dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

template<typename GV, typename RF>
class TwoPhaseParameter;

// initial conditions for phase pressures
template<typename GV, typename RF>
class P_l
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_l<GV,RF> >
{
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_l<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  P_l (const GV& gv, const TwoPhaseParameter<GV,RF>& tp_) : BaseT(gv), tp(tp_) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = tp.g_l_global(x,0);
    //    std::cout << x << " " << tp.g_l_global(x,0) << std::endl;
  }
};


template<typename GV, typename RF>
class P_g
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_g<GV,RF> >
{
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_g<GV,RF> > BaseT;

  P_g (const GV& gv_, const TwoPhaseParameter<GV,RF>& tp_) :
    BaseT(gv_), tp(tp_)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    RF s_l = tp.s_l(tp.g_g_global(x,0)-tp.g_l_global(x,0));
    RF p_l =  tp.g_l_global(x,0);

    y = tp.pc(s_l) + p_l;
    //y = tp.g_g_global(x,0);
  }
};



template<typename  T, typename PL, typename PG>
class P_c
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                   typename PL::Traits::RangeFieldType,
                                   PL::Traits::dimRange,
                                   typename PL::Traits::RangeType>,
  P_c<T,PL,PG> >
{
  const T& t;
  const PL& pl;
  const PG& pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,
                                           PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,P_c<T,PL,PG> > BaseT;

  P_c (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = pg_value - pl_value;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

template<typename  T, typename PL, typename PG>
class S_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_l<T,PL,PG> >
{
  const T& t;
  const PL& pl;
  const PG& pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_l<T,PL,PG> > BaseT;

  S_l (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = t.s_l(e,x,pg_value-pl_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

template<typename  T, typename PL, typename PG>
class S_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_g<T,PL,PG> >
{
  const T& t;
  const PL& pl;
  const PG& pg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_g<T,PL,PG> > BaseT;

  S_g (const T& t_, const PL& pl_, const PG& pg_) : t(t_), pl(pl_), pg(pg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType pl_value,pg_value;
    pl.evaluate(e,x,pl_value);
    pg.evaluate(e,x,pg_value);
    y = 1.0-t.s_l(e,x,pg_value-pl_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

#endif
