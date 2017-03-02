// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_INTERPOLATE_HH
#define DUNE_DYCAP_INTERPOLATE_HH

#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

/**
   \brief interpolate local vectors to global block vector

   * \tparam XF block vector with solution for every GFS
   * \tparam GFS GridFunctionSpace for concentration
   * \tparam X1 vector with solution for each GFS
   */
template<typename XF, typename GFS, typename X1, typename X2, typename X3, typename X4>
void interpolateSolutionToGFS (const XF& xf, GFS& gfs, X1& x1, X2& x2, X3& x3, X4& x4)
{
  #warning "interpolateSolutionToGFS is deprecated, use copy_dofs_parent_to_child!"
  /*
  Dune::PDELab::copy_dofs_parent_to_child(xf,x1,gfs,0);
  Dune::PDELab::copy_dofs_parent_to_child(xf,x2,gfs,1);
  Dune::PDELab::copy_dofs_parent_to_child(xf,x3,gfs,2);
  Dune::PDELab::copy_dofs_parent_to_child(xf,x4,gfs,3);
  */

  /*
  // get some types
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  //typedef typename GV::Traits::template Codim<0>::Entity Element;

  // make local function space
  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;

  LFS lfs(gfs);

  // map each cell to unique id
  Dune::PDELab::MultiGeomUniqueIDMapper<GV> cell_mapper(gfs.gridView());

  // loop once over the grid
  for (ElementIterator it = gfs.gridView().template begin<0>();
       it!=gfs.gridView().template end<0>(); ++it)
    {
      // bind local function space to element
      lfs.bind(*it);

      // compute unique id
      const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

      // copy values from solution for each component to block vector containing solution for every components
      x1[ids]=xf[ids][0];
      x2[ids]=xf[ids][1];
      x3[ids]=xf[ids][2];
      x4[ids]=xf[ids][3];
    }
  */
}

/**
   \brief interpolate global block vector to local vectors

   * \tparam XF block vector with solution for every GFS
   * \tparam GFS GridFunctionSpace for concentration
   * \tparam X1 vector with solution for each GFS
   */
template<typename XF, typename GFS, typename X1, typename X2, typename X3, typename X4>
void interpolateSolutionToPGFS (XF& xf, const GFS& gfs, const X1& x1, const X2& x2, const X3& x3, const X4& x4)
{

   #warning "interpolateSolutionToPGFS is deprecated, use copy_dofs_child_to_parent!"
  /*  Dune::PDELab::copy_dofs_child_to_parent(x1,xf,gfs,0);
  Dune::PDELab::copy_dofs_child_to_parent(x2,xf,gfs,1);
  Dune::PDELab::copy_dofs_child_to_parent(x3,xf,gfs,2);
  Dune::PDELab::copy_dofs_child_to_parent(x4,xf,gfs,3);
  */
   /*
  // get some types
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  //typedef typename GV::Traits::template Codim<0>::Entity Element;

  // make local function space
  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;

  LFS lfs(gfs);

  // map each cell to unique id
  Dune::PDELab::MultiGeomUniqueIDMapper<GV> cell_mapper(gfs.gridView());

  // loop once over the grid
  for (ElementIterator it = gfs.gridView().template begin<0>();
       it!=gfs.gridView().template end<0>(); ++it)
    {
      // bind local function space to element
      lfs.bind(*it);

      // compute unique id
      const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

      // copy values from solution for each component to block vector containing solution for every components
      xf[ids][0]=x1[ids];
      xf[ids][1]=x2[ids];
      xf[ids][2]=x3[ids];
      xf[ids][3]=x4[ids];
    }
  */
}
#endif
