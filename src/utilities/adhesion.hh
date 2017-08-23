// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifndef DUNE_DYCAP_ADHESION_HH
#define DUNE_DYCAP_ADHESION_HH

#include <vector>
#include <dune/common/parametertreeparser.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

/**
 * \todo{Write the general interface for bacterial growth}
 */

//! class to compute and store the bacterial adhesion
/**
 * \brief class to compute and store the bacterial adhesion
 *
 * \tparam Traits     traits to represent the data types
 * \tparam GFS        grid function space, bacterium
 * \tparam VEL        velocity discrete grid function, should represent pore velocity
 */
template <typename GFS, typename VEL>
class Adhesion
{
public:
  typedef typename GFS::Traits::GridViewType GridView;
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GridView,typename GridView::ctype,1> Traits;

  enum { dim = GridView::dimension};
  //! \brief export type for range field
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef RangeFieldType RF;

  //! \brief export type for doman field
  typedef typename Traits::DomainFieldType DomainFieldType;

  //! \brief export type for range type
  typedef typename Traits::RangeType RangeType;

  //! \brief the grid view
  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;

  //! \brief export type of the vector for the bacterial concentration
  typedef typename Dune::PDELab::Backend::Vector<GFS,RangeFieldType>::Type CV;

public:
  Adhesion(const GFS & gfs_, VEL & vel_, const Dune::ParameterTree param) :
    gfs(gfs_),
    lfs(gfs),
    vel(vel_),
    adhesion(gfs,0.),
    ma(param.sub("Adhesion").get<RF>("ma")),
    mb(param.sub("Adhesion").get<RF>("mb")),
    mc(param.sub("Adhesion").get<RF>("mc"))
  {
  }

  //! Update the adhesion from the old and new concentration (before and after growth)
  /*!
    \param cold concentration at the old time.
    \param cnew concentration at the new time.
    \return The test results
  */
  void update(CV & cold, CV & cnew)
  {
    // for each cell
    typedef typename GridView::template Codim <0>::Iterator ElementIterator;
    ElementIterator it = gfs.gridView().template begin<0>();
    ElementIterator endit = gfs.gridView().template end<0>();
    for (; it != endit; ++it)
      {

        lfs.bind(*it);
        unsigned int n = lfs.globalIndex(0);

        // compute cell value
        typename VEL::Traits::RangeType velo;

        // is not good, should bee
        vel.evaluate(*it, it->geometry().local(it->geometry().center()), velo);
        RangeFieldType adh = (cnew[n]-cold[n]);
        if (adh<0)
          adhesion[n]+=adh;
        else
          adhesion[n]+=evaluateAdhesion(velo)*adh;
      }
  }

  //! returns adhesion in %, see experiments and corresponding fit
  RangeFieldType evaluateAdhesion(typename VEL::Traits::RangeType & v)
  {
    RangeFieldType vabs = v.two_norm();
    RF f = ma/(vabs + mb) + mc;
    return std::min(f/100.,0.75);
  }

  //! returns adhesion in %, see experiments and corresponding fit
  RangeFieldType evaluateAdhesion(const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x)
  {
    // compute cell value
    typename VEL::Traits::RangeType velo;
    vel.evaluate(e,x, velo);
    return evaluateAdhesion(velo);
  }


  //! get the vector with adhesive bacterium
  CV & getAdhesion()
  {
    return adhesion;
  }

    inline const typename Traits::GridViewType& getGridView ()
  {
    return gfs.getGridView();
  }


private:
  // functionspace related variables
  const GFS& gfs;   //!< grid function space
  LFS lfs;          //!< local function space
  VEL & vel;        //!< pore velocity!
  CV adhesion;      //!< vector with adhesion
  const RF ma;
  const RF mb;
  const RF mc;
 };


//! adhesion visualization
template<typename ADH>
class AdhesionVisualize
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename ADH::Traits::GridViewType,
                                   typename ADH::Traits::RangeFieldType,
                                   1,
                                   typename ADH::Traits::RangeType>, AdhesionVisualize<ADH> >
{
  ADH& adh;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename ADH::Traits::GridViewType,
                                           typename ADH::Traits::RangeFieldType,
                                           1,
                                           typename ADH::Traits::RangeType> Traits;

  AdhesionVisualize (ADH& adh_) : adh(adh_)  {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y =adh.evaluateAdhesion(e,x);
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return adh.getGridView();
  }
};


#endif
