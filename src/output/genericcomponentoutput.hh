// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_GENERICCOMPONENTOUTPUT_HH
#define DUNE_DYCAP_GENERICCOMPONENTOUTPUT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <memory>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/typetree/typetree.hh>

namespace VTKGridFunctionImp{

  /** \brief TMP to fill a pvdwriter with all DGF representing all components
      using the corresponding sub grid view.
  */
  template<typename GFS, typename U, typename PVD, int k>
  struct FillVTKFunctions
  {
    static void fill
    (const std::string label, const GFS & gfs,
     const U & u, PVD & pvd)
    {
      typedef Dune::PDELab::GridFunctionSubSpace<GFS,Dune::TypeTree::TreePath<GFS::CHILDREN-k> > SubGFS;
      std::shared_ptr<SubGFS> subgfs(new SubGFS(gfs));

      typedef Dune::PDELab::DiscreteGridFunction<SubGFS,U> ClDGF;
      std::shared_ptr<ClDGF> cldgf(new ClDGF(subgfs,Dune::stackobject_to_shared_ptr(u)));

      char basename[255];
      sprintf(basename,"component_%u",GFS::CHILDREN-k);
      typedef Dune::PDELab::VTKGridFunctionAdapter<ClDGF> VTK;
      pvd.addCellData(std::shared_ptr<VTK>(new VTK(cldgf,basename)));

      FillVTKFunctions<GFS,U,PVD,k-1>::fill(label,gfs,u,pvd);
    }
  };

  //! default imlementation
  template<typename GFS, typename U, typename PVD>
  struct FillVTKFunctions<GFS, U, PVD, 0>
  {
    static void fill
    (const std::string label, const GFS & gfs,
     const U & u, PVD & pvd)
    {
    }
  };

} // end namespace

#endif
