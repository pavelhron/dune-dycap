// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef __DUNE_DYCAP_GFUTILITIES_HH__
#define __DUNE_DYCAP_GFUTILITIES_HH__

#include <dune/common/debugstream.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/common/fvector.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>

using namespace Dune::PDELab;

template<typename GF>
void maxGridFunction(const GF& gf,
                     typename GF::Traits::RangeType& sum,
                     unsigned qorder = 1) {
  typedef typename GF::Traits::GridViewType GV;
  typedef typename GV::template Codim<0>::
    template Partition<Dune::Interior_Partition>::Iterator EIterator;
  typedef typename GV::template Codim<0>::Geometry Geometry;
  typedef typename GF::Traits::RangeType Range;
  typedef typename GF::Traits::DomainFieldType DF;
  static const int dimD = GF::Traits::dimDomain;
  typedef Dune::QuadratureRule<DF,dimD> QR;
  typedef Dune::QuadratureRules<DF,dimD> QRs;
  typedef typename QR::const_iterator QIterator;

  sum = 0;
  Range val;
  const EIterator eend = gf.getGridView().template end<0,
                                                       Dune::Interior_Partition>();
  for(EIterator eit = gf.getGridView().template begin<0,
                                                      Dune::Interior_Partition>(); eit != eend; ++eit) {
    const Geometry& geo = eit->geometry();
    Dune::GeometryType gt = geo.type();
    const QR& rule = QRs::rule(gt,qorder);
    const QIterator qend = rule.end();

    for (QIterator qit=rule.begin(); qit != qend; ++qit)
      {
        // evaluate the given grid functions at integration point
        gf.evaluate(*eit,qit->position(),val);
        if (std::abs(val)>sum)
            sum = std::abs(val);
            }
      }
  }

  template<typename GV, typename RF>
    class ConstantDiscreteGridFunction
      : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>,
                                                      ConstantDiscreteGridFunction<GV,RF> >
  {
  public:
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ConstantDiscreteGridFunction<GV,RF> > BaseT;
    enum {dim=Traits::DomainType::dimension};

    ConstantDiscreteGridFunction (const GV& gv,  typename Traits::RangeFieldType v_) : BaseT(gv), v(v_)
    {}

    inline void evaluateGlobal (const typename Traits::DomainType& x,
                                typename Traits::RangeType& y) const
    {
      y = v;
    }
  private:
    typename Traits::RangeFieldType v;
  };


// constant analytic grid function space class
template<typename GV, typename RF>
class ConstADGF
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1>, ConstADGF<GV,RF> >
{

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,1> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,ConstADGF<GV,RF> > BaseT;
  enum {dim=Traits::DomainType::dimension};

  ConstADGF (const GV& gv) : BaseT(gv) {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
    y = 1.0;
  }
};


//! turn a std::vector of boundary data into a BoundaryGridFunction
template<typename GV, typename RF, int dimR, typename Mapper,typename Data>
class P0BoundaryGridFunction :
  public Dune::PDELab::BoundaryGridFunctionBase<
  Dune::PDELab::BoundaryGridFunctionTraits<GV, RF, dimR, Dune::FieldVector<RF, dimR>
                                           >,  P0BoundaryGridFunction<GV, RF, dimR, Mapper, Data> >
{
  typedef Dune::FieldVector<RF, dimR> Range;

public:
  typedef Dune::PDELab::BoundaryGridFunctionTraits<GV, RF, dimR, Range>
  Traits;

private:
  const GV &gv;
  const Mapper &mapper;
  const Data &data;

public:
  //! constructor
  P0BoundaryGridFunction(const GV &gv_,
                         const Mapper &mapper_,
                         const Data &data_) :
    gv(gv_), mapper(mapper_), data(data_)
  { }

  //! evaluate
  void evaluate(const Dune::PDELab::IntersectionGeometry
                <typename GV::Intersection> &ig,
                const typename Traits::DomainType &x,
                Range &y) const
  {
    y = data[mapper.map(*(ig.inside()), ig.indexInInside(), dimR)];
  }

  //! get gridview
  const GV &getGridView() const { return gv; }
};

//! turn a std::vector of boundary data into a BoundaryGridFunction
template<typename GV, typename RF, int dimR, typename Mapper,typename Data>
class P0GridFunction :
  public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<GV, RF, dimR, Dune::FieldVector<RF, dimR>
                                           >,  P0GridFunction<GV, RF, dimR, Mapper, Data> >
{
  typedef Dune::FieldVector<RF, dimR> Range;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV, RF, dimR, Range> Traits;

private:
  const GV &gv;
  const Mapper &mapper;
  const Data &data;

public:
  //! constructor
  P0GridFunction(const GV &gv_,
                 const Mapper &mapper_,
                 const Data &data_) :
    gv(gv_), mapper(mapper_), data(data_)
  { }

  void evaluate (const typename Traits::ElementType& e,
                 const typename Traits::DomainType& x,
                 typename Traits::RangeType& y) const
  {
    y = data[mapper.map(e)];
  }

  //! get gridview
  const GV &getGridView() const { return gv; }
};


// The difference between two grid functions
 template<typename GV, typename A, typename B>
   class MinusGridFunction
 {
 private:
   const A & a; // discrete grid function
   const B & b; // discrete grid function
   typedef typename A::Traits::DomainType DT;
   typedef typename A::Traits::RangeType RT;
   typedef typename GV::template Codim<0>::Entity Entity;

 public:

   typedef typename A::Traits Traits;

   MinusGridFunction(const A & a_, const B & b_)
     : a(a_), b(b_)
   {}

   void evaluate(const Entity & e, const DT & lx, RT & y) const
   {
     RT av; a.evaluate(e,lx,av);
     RT bv; b.evaluate(e,lx,bv);
     y = av - bv;
   }
 };




 template<typename GV, typename A, typename B>
   class PlusGridFunction
 {
 private:
   const A & a; // discrete grid function
   const B & b; // discrete grid function
   typedef typename A::Traits::DomainType DT;
   typedef typename A::Traits::RangeType RT;
   typedef typename GV::template Codim<0>::Entity Entity;

 public:

   typedef typename A::Traits Traits;

   PlusGridFunction(const A & a_, const B & b_)
     : a(a_), b(b_)
   {}

   void evaluate(const Entity & e, const DT & lx, RT & y) const
   {
     RT av; a.evaluate(e,lx,av);
     RT bv; b.evaluate(e,lx,bv);
     y = av + bv;
   }

   inline const typename Traits::GridViewType& getGridView () const
   {
     return a.getGridView();
   }
 };


// The difference between two grid functions
 template<typename GV, typename A, typename B>
   class MinGridFunction
 {
 private:
   const A & a; // discrete grid function
   const B & b; // discrete grid function
   typedef typename A::Traits::DomainType DT;
   typedef typename A::Traits::RangeType RT;
   typedef typename GV::template Codim<0>::Entity Entity;

 public:

   typedef typename A::Traits Traits;

   MinGridFunction(const A & a_, const B & b_)
     : a(a_), b(b_)
   {}

   void evaluate(const Entity & e, const DT & lx, RT & y) const
   {
     RT av; a.evaluate(e,lx,av);
     RT bv; b.evaluate(e,lx,bv);
     y = std::min(av,bv);
   }
 };

 // The product of two grid functions
 template<typename GV,typename A, typename B = ConstantDiscreteGridFunction<GV, typename A::Traits::RangeFieldType> >
   class ProductGridFunction
 {
 private:
   const A & a; // discrete grid function
   const B & b; // discrete grid function
   typedef typename A::Traits::DomainType DT;
   typedef typename A::Traits::RangeType RT;
   const RT weight;
   typedef typename GV::template Codim<0>::Entity Entity;

 public:

   typedef typename A::Traits Traits;

   ProductGridFunction(const A & a_, const B & b_, const typename A::Traits::RangeType weight_=1.0)
     : a(a_), b(b_), weight(weight_)
   {}

   ProductGridFunction(const A & a_, const typename A::Traits::RangeType weight_=1.0)
     : a(a_), b(B(a.getGridView(),1.0)), weight(weight_)
   {
   }

   void evaluate(const Entity & e, const DT & lx, RT & y) const
   {
     RT av; a.evaluate(e,lx,av);
     RT bv; b.evaluate(e,lx,bv);
     y = av * bv * weight;
   }

   inline const typename Traits::GridViewType& getGridView () const
   {
     return a.getGridView();
   }

 };



 /*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

   \tparam T1  a grid function type
   \tparam T2  a grid function type
 */
 template<typename T1, typename T2>
   class DifferenceAdapter
     : public Dune::PDELab::GridFunctionBase<
       Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                        typename T1::Traits::RangeFieldType,
                                        1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
       ,DifferenceAdapter<T1,T2> >
 {
 public:
   typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                            typename T1::Traits::RangeFieldType,
                                            1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

   //! constructor
   DifferenceAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

   //! \copydoc GridFunctionBase::evaluate()
   inline void evaluate (const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x,
                         typename Traits::RangeType& y) const
   {
     typename T1::Traits::RangeType y1;
     t1.evaluate(e,x,y1);
     typename T2::Traits::RangeType y2;
     t2.evaluate(e,x,y2);
     y1 -= y2;
     y = y1;
   }

   inline const typename Traits::GridViewType& getGridView () const
   {
     return t1.getGridView();
   }

 private:
   const T1& t1;
   const T2& t2;
 };


 /*! \brief Adapter returning f1(x)-f2(x) for two given grid functions

   \tparam T1  a grid function type
   \tparam T2  a grid function type
 */
 template<typename T1, typename T2>
   class DifferenceAbsoluteAdapter
     : public Dune::PDELab::GridFunctionBase<
       Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                        typename T1::Traits::RangeFieldType,
                                        1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
       ,DifferenceAbsoluteAdapter<T1,T2> >
 {
 public:
   typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                            typename T1::Traits::RangeFieldType,
                                            1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

   //! constructor
   DifferenceAbsoluteAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

   //! \copydoc GridFunctionBase::evaluate()
   inline void evaluate (const typename Traits::ElementType& e,
                         const typename Traits::DomainType& x,
                         typename Traits::RangeType& y) const
   {
     typename T1::Traits::RangeType y1;
     t1.evaluate(e,x,y1);
     typename T2::Traits::RangeType y2;
     t2.evaluate(e,x,y2);
     y1 -= y2;
     y = y1.one_norm();
   }

   inline const typename Traits::GridViewType& getGridView () const
   {
     return t1.getGridView();
   }

 private:
   const T1& t1;
   const T2& t2;
 };


template<typename T1, typename Flux>
class FluxAdapter
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,FluxAdapter<T1,Flux> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor
  FluxAdapter (const T1& t1_, const Flux& flux_) : t1(t1_), flux(flux_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename Traits::RangeType y1;
    //t1.evaluate(e,x,y1);

    //typename Traits::RangeType y2;
    flux.evaluate(e,x,y1,y[0]);
    //std::cout << y[0] << std::endl;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }

private:
  const T1& t1;
  const Flux& flux;
};

/*! \brief Adapter returning ||f1(x)-f2(x)||^2 for two given grid functions

  \tparam T1  a grid function type
  \tparam T2  a grid function type
*/
template<typename T1, typename T2>
class DifferenceSquaredAdapter
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                   typename T1::Traits::RangeFieldType,
                                   1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> >
  ,DifferenceSquaredAdapter<T1,T2> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename T1::Traits::GridViewType,
                                           typename T1::Traits::RangeFieldType,
                                           1,Dune::FieldVector<typename T1::Traits::RangeFieldType,1> > Traits;

  //! constructor
  DifferenceSquaredAdapter (const T1& t1_, const T2& t2_) : t1(t1_), t2(t2_) {}

  //! \copydoc GridFunctionBase::evaluate()
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename T1::Traits::RangeType y1;
    t1.evaluate(e,x,y1);
    typename T2::Traits::RangeType y2;
    t2.evaluate(e,x,y2);
    y1 -= y2;
    y = y1.two_norm2();
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return t1.getGridView();
  }

private:
  const T1& t1;
  const T2& t2;
};





/**
   \brief This grid function computes the solution which is defined
   on coarse grid for elements in the finer grid. Can be useful if you want
   to interpolate solution from coarser grid to finer grid.

   The functions are expected to be defined on level grid views and
   the reference function is expected to be defined on a higher
   level. The evaluate function should only be called for elements
   from the grid view of the reference grid function.

   \tparam RF Reference grid function (on a fine level)
   \tparam TF Test grid function (on a coarser level than RF)
*/
template<typename RF, typename TF>
class MultiLevelDifferenceAdapter
  : public Dune::TypeTree::LeafNode, Dune::PDELab::GridFunctionInterface<
  Dune::PDELab::GridFunctionTraits<
    typename RF::Traits::GridViewType,
    typename RF::Traits::RangeFieldType,
    1,
    Dune::FieldVector<typename RF::Traits::RangeFieldType,2>
    >,
  MultiLevelDifferenceAdapter<RF,TF>
  >
{
public:

  typedef Dune::PDELab::GridFunctionTraits<
  typename RF::Traits::GridViewType,
  typename RF::Traits::RangeFieldType,
  2,
  Dune::FieldVector<typename RF::Traits::RangeFieldType,2>
  > Traits;

  MultiLevelDifferenceAdapter (const RF & rf_, const TF & tf_)
    : rf(rf_), tf(tf_), tf_level(tf.getGridView().template begin<0>()->level()),
      rf_level(rf.getGridView().template begin<0>()->level())
  {}

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    assert(e.level() == rf_level);
    typename RF::Traits::RangeType ry;
    rf.evaluate(e,x,ry);

    typedef typename Traits::ElementType::EntityPointer EP;
    EP ep(e);
    while(ep->level() != tf_level){
      assert(ep->hasFather());
      ep = ep->father();
    }

    typename TF::Traits::RangeType ty;
    typename Traits::DomainType tx = ep->geometry().local(e.geometry().global(x));
    tf.evaluate(*ep,tx,ty);

    ry-=ty;
    y = ry;
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView () const
  {
    return rf.getGridView();
  }

private:
  const RF & rf;
  const TF & tf;
  const int tf_level;
  const int rf_level;
};


/**
   \brief This grid function computes the squared difference of two given
   grid functions.

   The functions are expected to be defined on level grid views and
   the reference function is expected to be defined on a higher
   level. The evaluate function should only be called for elements
   from the grid view of the reference grid function.

   \tparam RF Reference grid function (on a fine level)
   \tparam TF Test grid function (on a coarser level than RF)
*/
template<typename RF, typename TF>
class MultiLevelDifferenceSquaredAdapter
  : public Dune::TypeTree::LeafNode, Dune::PDELab::GridFunctionInterface<
  Dune::PDELab::GridFunctionTraits<
    typename RF::Traits::GridViewType,
    typename RF::Traits::RangeFieldType,
    1,
    Dune::FieldVector<typename RF::Traits::RangeFieldType,1>
    >,
  MultiLevelDifferenceSquaredAdapter<RF,TF>
  >
{
public:

  typedef Dune::PDELab::GridFunctionTraits<
  typename RF::Traits::GridViewType,
  typename RF::Traits::RangeFieldType,
  1,
  Dune::FieldVector<typename RF::Traits::RangeFieldType,1>
  > Traits;

  MultiLevelDifferenceSquaredAdapter (const RF & rf_, const TF & tf_)
    : rf(rf_), tf(tf_), tf_level(tf.getGridView().template begin<0>()->level()),
      rf_level(rf.getGridView().template begin<0>()->level())
  {}

  // Evaluate
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    assert(e.level() == rf_level);
    typename RF::Traits::RangeType ry;
    rf.evaluate(e,x,ry);

    typedef typename Traits::ElementType::EntityPointer EP;
    EP ep(e);
    while(ep->level() != tf_level){
      assert(ep->hasFather());
      ep = ep->father();
    }
    typename TF::Traits::RangeType ty;
    typename Traits::DomainType tx = ep->geometry().local(e.geometry().global(x));
    tf.evaluate(*ep,tx,ty);

    ry -= ty;
    y = ry * ry;
  }

  //! get a reference to the GridView
  inline const typename Traits::GridViewType& getGridView () const
  {
    return rf.getGridView();
  }

private:
  const RF & rf;
  const TF & tf;
  const int tf_level;
  const int rf_level;
};

/**
   \brief This grid function computes the solution which is defined
   on coarse grid for elements in the finer grid. Can be useful if you want
   to interpolate solution from coarser grid to finer grid.

   The functions are expected to be defined on level grid views and
   the reference function is expected to be defined on a higher
   level. The evaluate function should only be called for elements
   from the grid view of the reference grid function.

   \tparam RF Reference grid function (on a fine level)
   \tparam TF Test grid function (on a coarser level than RF)
*/
template<typename LGV, typename TF>
class MultiLevelAdapter
  : public Dune::TypeTree::LeafNode, Dune::PDELab::GridFunctionInterface<
  Dune::PDELab::GridFunctionTraits<
    LGV,
    typename TF::Traits::RangeFieldType,
    1,
    Dune::FieldVector<typename TF::Traits::RangeFieldType,1>
    >,
  MultiLevelAdapter<LGV,TF>
  >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<
  LGV,
  typename TF::Traits::RangeFieldType,
  1,
  Dune::FieldVector<typename TF::Traits::RangeFieldType,1>
  > Traits;


MultiLevelAdapter (const LGV & gv_, const TF & tf_)
  : gv(gv_), tf(tf_), tf_level(tf.getGridView().template begin<0>()->level()), rf_level(gv.template begin<0>()->level())

{}

// Evaluate
inline void evaluate (const typename Traits::ElementType& e,
                      const typename Traits::DomainType& x,
                      typename Traits::RangeType& y) const
{
  assert(e.level() == rf_level);

  typedef typename Traits::ElementType::EntityPointer EP;
  EP ep(e);
  while(ep->level() != tf_level){
    assert(ep->hasFather());
    ep = ep->father();
  }

  typename TF::Traits::RangeType ty;
  typename Traits::DomainType tx = ep->geometry().local(e.geometry().global(x));
  tf.evaluate(*ep,tx,ty);

  y = ty;
}

//! get a reference to the GridView
inline const typename Traits::GridViewType& getGridView () const
{
  return gv;
}


private:
const LGV & gv;
const TF & tf;
const int tf_level;
const int rf_level;
};


  /** \brief DiscreteGridFunction for vector-valued functions
     *
     * convert a power function space of scalar function spaces into a
     * vector-valued grid function this is just an intermediate
     * solution to provide VTK output
     *
     * \tparam T Type of PowerGridFunctionSpace
     * \tparam X Type of coefficients vector
     * \tparam dimR Force a different number of components for the resulting
     *              GridFunction than the PowerGridFunctionSpace.
     */
    template<typename T, typename X, std::size_t dimR = T::CHILDREN>
    class DiscreteGridFunctionFromPGF
      : public GridFunctionInterface<
          GridFunctionTraits<
            typename T::Traits::GridViewType,
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            1,
            Dune::FieldVector<
              typename T::template Child<0>::Type::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType,
              1
              >
            >,
          DiscreteGridFunctionFromPGF<T,X>
          >,
        public Dune::TypeTree::LeafNode
    {
      typedef T GFS;

      typedef GridFunctionInterface<
        GridFunctionTraits<
          typename T::Traits::GridViewType,
          typename T::template Child<0>::Type::Traits::FiniteElementType
                   ::Traits::LocalBasisType::Traits::RangeFieldType,
          1,
          Dune::FieldVector<
            typename T::template Child<0>::Type::Traits::FiniteElementType
                     ::Traits::LocalBasisType::Traits::RangeFieldType,
            1
            >
          >,
        DiscreteGridFunctionFromPGF<T,X,dimR>
        > BaseT;

    public:
      typedef typename BaseT::Traits Traits;
      typedef typename T::template Child<0>::Type ChildType;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeFieldType RF;
      typedef typename ChildType::Traits::FiniteElementType
                       ::Traits::LocalBasisType::Traits::RangeType RT;

      //! construct
      /**
       * \param gfs   GridFunctionSpace.
       * \param x_    Coefficient vector.
       * \param start Number of first child of gfs to use.
       */
      DiscreteGridFunctionFromPGF(const GFS& gfs, const X& x_,
                                  std::size_t k_)
      : pgfs(stackobject_to_shared_ptr(gfs))
      , k(k_)
      , lfs(gfs)
      , lfs_cache(lfs)
      , x_view(x_)
      , xl(gfs.maxLocalSize())
      , yb(gfs.maxLocalSize())
      {
        for(std::size_t i = 0; i < dimR; ++i)
          remap[i] = i;
      }


      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        lfs.bind(e);
        lfs_cache.update();
         x_view.bind(lfs_cache);
        x_view.read(xl);
        x_view.unbind();
        lfs.child(remap[k]).finiteElement().localBasis().
          evaluateFunction(x,yb);
        y = 0.0;
        for (unsigned int i=0; i<yb.size(); i++)
          y += xl[lfs.child(remap[k]).localIndex(i)]*yb[i];

      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return pgfs->gridView();
      }

    private:
      typedef LocalFunctionSpace<GFS> LFS;
      typedef LFSIndexCache<LFS> LFSCache;
      typedef typename X::template ConstLocalView<LFSCache> XView;

      Dune::shared_ptr<GFS const> pgfs;
      const size_t k;
      std::size_t remap[dimR];
      mutable LFS lfs;
      mutable LFSCache lfs_cache;
      mutable XView x_view;
      mutable std::vector<RF> xl;
      mutable std::vector<RT> yb;
      Dune::shared_ptr<const X> px; // FIXME: dummy pointer to make sure we take ownership of X
    };

#endif // __DUNE_DYCAP_GFUTILITIES_HH__
