// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TWOPHASEVELOCITY_BOILERPLATE_HH
#define DUNE_DYCAP_TWOPHASEVELOCITY_BOILERPLATE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>

#include<src/utilities/rt0qfem.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>

#include<src/models/twophasevelocity.hh>
#include<src/output/output.hh>

namespace Dune {
  namespace Dycap{

    template<typename GV, typename RF, typename DF, typename TP, typename P_lDGF, typename P_cDGF, bool vflux = false>
    class TwoPhaseVelocity
    {
    public:
      enum {dim = GV::dimension};
      typedef VelocityLiquid<TP,P_lDGF,P_cDGF> VliquidDGF;
      typedef VelocityGas<TP,P_lDGF,P_cDGF> VgasDGF;

      typedef Dune::PDELab::PiolaBackwardAdapter<VliquidDGF> VliquidBACK;
      typedef Dune::PDELab::PiolaBackwardAdapter<VgasDGF> VgasBACK;
      typedef Dune::PDELab::RT0QLocalFiniteElementMap<GV,DF,RF,dim> RT0FEM;

      typedef Dune::PDELab::GridFunctionSpace<GV,RT0FEM,Dune::PDELab::P0ParallelConstraints,Dune::PDELab::istl::VectorBackend<> > RT0GFS;

      typedef typename Dune::PDELab::Backend::Vector<RT0GFS,RF> RT0VEC;
      typedef Dune::PDELab::DiscreteGridFunctionPiola<RT0GFS,RT0VEC> RT0DGF;

      TwoPhaseVelocity(const GV &gv_, TP & tp, P_lDGF & p_ldgf, P_cDGF & p_cdgf)
        : gv(gv_),
          vliquiddgf(tp,p_ldgf,p_cdgf,vflux),
          vgasdgf(tp,p_ldgf,p_cdgf,vflux),
          vliquiddgf_back(vliquiddgf),
          vgasdgf_back(vgasdgf),
          rt0fem(gv),
          rt0gfs(gv,rt0fem),
          rt0vec_l(rt0gfs,0.0),
          rt0vec_zero(rt0gfs,0.0),
          rt0vec_g(rt0gfs,0.0),
          rt0vec_l_old(rt0gfs,0.0),
          rt0vec_g_old(rt0gfs,0.0),
          rt0_l(rt0gfs,rt0vec_l),
          rt0_zero(rt0gfs,rt0vec_zero),
          rt0_g(rt0gfs,rt0vec_g),
          rt0_l_old(rt0gfs,rt0vec_l_old),
          rt0_g_old(rt0gfs,rt0vec_g_old),
          verbosity(1)
      {
        vliquiddgf.set_time(0.);
        vgasdgf.set_time(0.);

        if (gv.comm().rank()>0)
          verbosity = 0;
        interpolate(0.);
      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter, int vtklevel = 1)
      {
        typedef Dune::PDELab::VTKGridFunctionAdapter<RT0DGF> VTK;
        pvdwriter.addVertexData(std::shared_ptr<VTK>(new VTK(rt0_l,"liquid velocity")));
        pvdwriter.addVertexData(std::shared_ptr<VTK>(new VTK(rt0_g,"gas velocity")));
        if (vtklevel>1)
          {
            pvdwriter.addVertexData(std::shared_ptr<VTK>(new VTK(rt0_l_old,"liquid velocity old")));
            pvdwriter.addVertexData(std::shared_ptr<VTK>(new VTK(rt0_g_old,"gas velocity old")));
          }
        if (vtklevel>2)
          {
            typedef Dune::PDELab::VTKGridFunctionAdapter<VgasDGF> VTK2;
            typedef Dune::PDELab::VTKGridFunctionAdapter<VliquidDGF> VTK3;

            pvdwriter.addVertexData(std::shared_ptr<VTK2>(new VTK2(vgasdgf,"gas")));
            pvdwriter.addVertexData(std::shared_ptr<VTK3>(new VTK3(vliquiddgf,"liquid")));
          }
      }

      // interpolate new velocities to RT0 space
      void interpolate(RF t)
      {
        if (verbosity)
          std::cout << "compute velocities at time "<< t << " \n";
        watch.reset();
        rt0vec_l_old = rt0vec_l;
        rt0vec_g_old = rt0vec_g;


        vliquiddgf.set_time(t);
        vgasdgf.set_time(t);


      //  Dune::PDELab::CopyDataHandle<RT0GFS,RT0VEC> copydh(rt0gfs,rt0vec_l);
      //  if (rt0gfs.gridView().comm().size()>1) // collect correction
      //    rt0gfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
        Dune::PDELab::interpolate(vgasdgf_back,rt0gfs,rt0vec_g);
        Dune::PDELab::interpolate(vliquiddgf_back,rt0gfs,rt0vec_l);

        if (verbosity)
          std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
      }

      void zeroVelocity(bool both)
      {
        rt0vec_l = 0.0;
        rt0vec_l_old = 0.0;
        if (both)
          {
            rt0vec_g = 0.0;
            rt0vec_g_old = 0.0;
          }
      }


      RT0DGF & getVl(){return rt0_l;}
      const RT0DGF & getVl() const {return rt0_l;}

      RT0DGF & getVg(){return rt0_g;}
      const RT0DGF & getVg() const {return rt0_g;}

      RT0DGF & getVlOld(){return rt0_l_old;}
      const RT0DGF & getVlOld() const {return rt0_l_old;}

      RT0DGF & getVgOld(){return rt0_g_old;}
      const RT0DGF & getVgOld() const {return rt0_g_old;}

      RT0DGF & getVzero(){return rt0_zero;}
      const RT0DGF & getVzero() const {return rt0_zero;}

    private:
      const GV& gv;
      Dune::Timer watch;

      VliquidDGF vliquiddgf;
      VgasDGF vgasdgf;

      VliquidBACK vliquiddgf_back;
      VgasBACK vgasdgf_back;

      RT0FEM rt0fem;
      RT0GFS rt0gfs;

      RT0VEC rt0vec_l;
      RT0VEC rt0vec_zero;
      RT0VEC rt0vec_g;

      RT0VEC rt0vec_l_old;
      RT0VEC rt0vec_g_old;

      RT0DGF rt0_l;
      RT0DGF rt0_zero;
      RT0DGF rt0_g;

      RT0DGF rt0_l_old;
      RT0DGF rt0_g_old;

      int verbosity;
    };


  } // end namespace Dycap
} // end namespace Dune

#endif
