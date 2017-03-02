// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PLSG_HH
#define DUNE_DYCAP_PLSG_HH

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
  }
};


template<typename GV, typename RF>
class S_g
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                  S_g<GV,RF> >
{
  const TwoPhaseParameter<GV,RF>& tp;
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,S_g<GV,RF> > BaseT;

  S_g (const GV& gv_, const TwoPhaseParameter<GV,RF>& tp_) :
    BaseT(gv_), tp(tp_)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    RF p_g = tp.g_g_global(x,0);
    RF p_l = tp.g_l_global(x,0);
    y = 1. - tp.s_l(p_g-p_l);
  }
};



template<typename  T, typename SG>
class S_l
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename SG::Traits::GridViewType,
                                                                           typename SG::Traits::RangeFieldType,SG::Traits::dimRange,
                                                                           typename SG::Traits::RangeType>,
                                          S_l<T,SG> >
{
  const T& t;
  const SG& sg;


public:
  typedef Dune::PDELab::GridFunctionTraits<typename SG::Traits::GridViewType,
                                           typename SG::Traits::RangeFieldType,SG::Traits::dimRange,
                                           typename SG::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,S_l<T,SG> > BaseT;

  S_l (const T& t_, const SG& sg_) : t(t_), sg(sg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename SG::Traits::RangeType sg_value;
    sg.evaluate(e,x,sg_value);

    y = 1. - sg_value;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return sg.getGridView();
  }
};


template<typename  T, typename SG>
class Zero
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename SG::Traits::GridViewType,
                                                                           typename SG::Traits::RangeFieldType,SG::Traits::dimRange,
                                                                           typename SG::Traits::RangeType>,
                                          Zero<T,SG> >
{
  const T& t;
  const SG& sg;


public:
  typedef Dune::PDELab::GridFunctionTraits<typename SG::Traits::GridViewType,
                                           typename SG::Traits::RangeFieldType,SG::Traits::dimRange,
                                           typename SG::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,Zero<T,SG> > BaseT;

  Zero (const T& t_, const SG& sg_) : t(t_), sg(sg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y = 1.;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return sg.getGridView();
  }
};


template<typename  T, typename PL, typename SG>
class P_c
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          P_c<T,PL,SG> >
{
  const T& t;
  const PL& pl;
  const SG& sg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,P_c<T,PL,SG> > BaseT;

  P_c (const T& t_, const PL& pl_, const SG& sg_) : t(t_), pl(pl_), sg(sg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType sg_value;

    sg.evaluate(e,x,sg_value);

    y = t.pc(e,x,1.-sg_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};


template<typename  T, typename PL, typename SG>
class P_g
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                                                           typename PL::Traits::RangeType>,
                                          P_g<T,PL,SG> >
{
  const T& t;
  const PL& pl;
  const SG& sg;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                           typename PL::Traits::RangeFieldType,PL::Traits::dimRange,
                                           typename PL::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,P_g<T,PL,SG> > BaseT;

  P_g (const T& t_, const PL& pl_, const SG& sg_) : t(t_), pl(pl_), sg(sg_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename PL::Traits::RangeType sg_value, pl_value;

    sg.evaluate(e,x,sg_value);
    pl.evaluate(e,x,pl_value);

    y = pl_value+t.pc(e,x,1.-sg_value);
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return pl.getGridView();
  }
};


template<typename V>
class VelocityX
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename V::Traits::GridViewType,
                                                                           typename V::Traits::RangeFieldType,V::Traits::dimRange,
                                                                           typename V::Traits::RangeType>,
                                          VelocityX<V> >
{
  const V& v;


public:
  typedef Dune::PDELab::GridFunctionTraits<typename V::Traits::GridViewType,
                                           typename V::Traits::RangeFieldType,V::Traits::dimRange,
                                           typename V::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,VelocityX<V> > BaseT;

  VelocityX (const V& v_) : v(v_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename Traits::RangeType velo;
    v.evaluate(e,x,velo);
    y = velo;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return v.getGridView();
  }
};

#endif
