
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PLPC_HH
#define DUNE_DYCAP_PLPC_HH

#include <dune/pdelab/common/function.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

// initial conditions for liquid pressure
template<typename GV, typename RF, typename TP>
class P_l
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_l<GV,RF,TP> >
{
  const TP& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_l<GV,RF,TP> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  P_l (const GV& gv, const TP& tp_) : BaseT(gv), tp(tp_) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = tp.g_l_global(x,0);
  }
};

// initial conditions for capillary pressure
template<typename GV, typename RF, typename TP>
class P_c
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  P_c<GV,RF,TP> >
{
  const TP& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,P_c<GV,RF,TP> > BaseT;

  P_c (const GV& gv_, const TP& tp_) :
    BaseT(gv_), tp(tp_)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    RF s_l = tp.s_l(tp.g_g_global(x,0)-tp.g_l_global(x,0));
    y = std::max(tp.pc(s_l),0.);
    //y=tp.g_g_global(x,0)-tp.g_l_global(x,0);
  }
};


/********************************************************************/
/************** saturation and gas pressure output ******************/
/********************************************************************/

template<typename  T, typename PL, typename PC>
class P_g
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename PC::Traits::GridViewType,
                                   typename PC::Traits::RangeFieldType,
                                   PC::Traits::dimRange,
                                   typename PC::Traits::RangeType>,
  P_g<T,PL,PC> >
{
  const T& t;
  const PL& pl;
  const PC& pc;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PC::Traits::GridViewType,
                                           typename PC::Traits::RangeFieldType,
                                           PC::Traits::dimRange,
                                           typename PC::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,P_g<T,PL,PC> > BaseT;

  P_g (const T& t_, const PL& pl_, const PC& pc_) : t(t_), pl(pl_), pc(pc_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PC::Traits::RangeType pl_value,pc_value;
    pl.evaluate(e,x,pl_value);
    pc.evaluate(e,x,pc_value);
    y = pc_value + pl_value;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};



template<typename  T, typename PL, typename PC>
class S_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_l<T,PL,PC> >
{
  const T& t;
  const PL& pl;
  const PC& pc;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_l<T,PL,PC> > BaseT;

  S_l (const T& t_, const PL& pl_, const PC& pc_) : t(t_), pl(pl_), pc(pc_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType pc_value;
    pc.evaluate(e,x,pc_value);
    y = t.s_l(e,x,pc_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

template<typename  T, typename PL, typename PC>
class S_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          S_g<T,PL,PC> >
{
  const T& t;
  const PL& pl;
  const PC& pc;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_g<T,PL,PC> > BaseT;

  S_g (const T& t_, const PL& pl_, const PC& pc_) : t(t_), pl(pl_), pc(pc_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType pc_value;
    pc.evaluate(e,x,pc_value);
    y = 1.0-t.s_l(e,x,pc_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};

#endif
