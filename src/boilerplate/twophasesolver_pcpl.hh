// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TWOPHASESOLVER_PCPL_BOILERPLATE_HH
#define DUNE_DYCAP_TWOPHASESOLVER_PCPL_BOILERPLATE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>

//#include<dune/pdelab/"../oldpdelab/pdelab.hh"
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/backend/istl.hh>

#include<src/newton/newton.hh>
#include<src/newton/newton_utilities.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<src/utilities/utilities.hh>

#include<src/models/pcpl_twophaseccfv_incompressible.hh>
#include<src/parameters/plpc.hh>
#include<src/output/output.hh>


#define PLPC TRUE

#if PLPC
#include<src/models/pcpl_twophaseccfv_incompressible.hh>
#include<src/parameters/plpc.hh>

#elif PLPG
#include<src/models/twophaseccfv.hh>
#include<src/parameters/plpg.hh>

#elif PLSG
#include<src/models/plsg_twophaseccfv.hh>
#include<src/parameters/plsg.hh>

#elif PLPCFLUX
#include<src/models/pcpl_twophaseccfv_totalflux.hh>
#include<src/parameters/plpc.hh>
#endif

#include<src/parameters/twophaseparameters.hh>

#include<src/output/output.hh>
#include"twophasevelocity.hh"



namespace Dune {
  namespace Dycap{

    template<typename GV, typename TP>
    class TwoPhaseSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE0;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;
      typedef Dune::PDELab::istl::VectorBackend
      <Dune::PDELab::istl::Blocking::fixed,2> VBE;
      typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,VBE,
                                                   Dune::PDELab::EntityBlockedOrderingTag> TPGFS;

      // subspaces for visualization
      typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0>> P_lSUB;
      typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1>> P_cSUB;

      // grid operator types
      typedef typename TPGFS::template ConstraintsContainer<RF>::Type C;
      typedef Dune::PDELab::TwoPhaseTwoPointFluxOperator<TP> LOP;
      typedef Dune::PDELab::TwoPhaseOnePointTemporalOperator<TP> MLOP;
      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
      typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
      typedef typename IGO::Traits::Domain V;

      // discrete grid functions
      typedef Dune::PDELab::DiscreteGridFunction<P_lSUB,V> P_lDGF;
      typedef Dune::PDELab::DiscreteGridFunction<P_cSUB,V> P_cDGF;
      typedef S_l<TP,P_lDGF,P_cDGF> S_lDGF;
      typedef S_g<TP,P_lDGF,P_cDGF> S_gDGF;
      typedef P_g<TP,P_lDGF,P_cDGF> P_gDGF;

      //inital grid functions
      typedef P_l<GV,RF,TP> P_lType;
      typedef P_c<GV,RF, TP> P_cType;
      typedef Dune::PDELab::CompositeGridFunction<P_lType,P_cType> PType;

      typedef Dune::PDELab::ImplicitEulerParameter<RF> TSMETHOD;
      typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
      //typedef ISTLBackend_OVLP_BCGS_SuperLU<TPGFS,C> LS;
      // typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<TPGFS,C> LS;
      // typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<TPGFS,C> LS;
      //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILU0<TPGFS,C> LS;

      typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
      typedef Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,V,V> OSM;


      typedef TwoPhaseVelocity<GV,RF,DF,TP,P_lDGF,P_cDGF> TPV;

      typedef typename TPV::RT0DGF VELDGF;
      typedef PermeabilityField<TP> PF;


      /** \brief Constructor
          \param [in] gv_ The grid view
          \param [in] fem_ The finite element map
          \param [in] level_ The level of the grif
      */
      TwoPhaseSolver(const GV &gv_, TP& tp_,  const Dune::ParameterTree & param_)
        : param(param_),
          gv(gv_), tp(tp_),
          verbosity(1),
          graphics(param.get<bool>("writeVTKTwoPhase", true)),
          vtklevel(param.get<int>("vtklevel", 1)),
          fem(getCube()),
          con(),
          gfs(gv,fem,con), tpgfs(gfs), p_lsub(tpgfs), p_csub(tpgfs),
          lop(tp), mlop(tp), mbe(5),
          cg(getContainer()),
          go0(tpgfs,cg,tpgfs,cg,lop,mbe),
          go1(tpgfs,cg,tpgfs,cg,mlop,mbe),
          igo(go0,go1),
          pold(tpgfs, 0.0),
          pnew(tpgfs, 0.0),
          p_ldgf(p_lsub,pnew),  p_cdgf(p_csub,pnew), p_gdgf(tp,p_ldgf,p_cdgf),
          pold_ldgf(p_lsub,pold),  pold_cdgf(p_csub,pold), pold_gdgf(tp,pold_ldgf,pold_cdgf),
          s_ldgf(tp,p_ldgf,p_cdgf),
          sold_ldgf(tp,pold_ldgf,pold_cdgf),
          s_gdgf(tp,p_ldgf,p_cdgf),
          sold_gdgf(tp,pold_ldgf,pold_cdgf),
          pf(tp),
          tsmethod(),
          //ls(tpgfs,cg,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0)), // superlu
          ls(tpgfs,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0),false), // amg
          newton(igo,ls),
          osm(tsmethod,igo,newton),
          tpv(gv,tp,p_ldgf,p_cdgf),
          timestepcounter(0)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        // initial value function
        typedef P_l<GV,RF,TP> P_lType;
        P_lType p_l_initial(gv,tp);
        typedef P_c<GV,RF,TP> P_cType;
        P_cType p_c_initial(gv,tp);
        typedef Dune::PDELab::CompositeGridFunction<P_lType,P_cType> PType;
        PType p_initial(p_l_initial,p_c_initial);
        Dune::PDELab::constraints(p_initial,tpgfs,cg);
        Dune::PDELab::interpolate(p_initial,tpgfs,pold);
        pnew=pold;

        // Newton parameters class
        typedef Dune::PDELab::NewtonParameters NewtonParameters;
        NewtonParameters newtonparameters(param.sub("Newton"));
        // set new parameters
        newtonparameters.set(newton);
        osm.setVerbosityLevel(param.sub("Verbosity").get<int>("Instationary", 0));
      }

      static GeometryType getCube() {
        GeometryType gt;
        gt.makeCube(dim);
        return gt;
      }

      static C getContainer() {
        C cg;
        cg.clear();
        return cg;
      }

      P_lDGF & getP_ldgf(){return p_ldgf;}
      const P_lDGF & getP_ldgf() const {return p_ldgf;}

      P_cDGF & getP_cdgf(){return p_cdgf;}
      const P_cDGF & getP_cdgf() const {return p_cdgf;}

      P_gDGF & getP_gdgf(){return p_gdgf;}
      const P_gDGF & getP_gdgf() const {return p_gdgf;}

      S_lDGF & getS_ldgf(){return s_ldgf;}
      const S_lDGF & getS_ldgf() const {return s_ldgf;}

      S_gDGF & getS_gdgf(){return s_gdgf;}
      const S_gDGF & getS_gdgf() const {return s_gdgf;}

      VELDGF & getVl(){return tpv.getVl();}
      const VELDGF & getVl() const {return tpv.getVl();}

      VELDGF & getVzero(){return tpv.getVzero();}
      const VELDGF & getVzero() const {return tpv.getVzero();}

      VELDGF & getVg(){return tpv.getVg();}
      const VELDGF & getVg() const {return tpv.getVg();}

      P_lDGF & getP_ldgfOld(){return pold_ldgf;}
      const P_lDGF & getP_ldgfOld() const {return pold_ldgf;}

      P_cDGF & getP_cdgfOld(){return pold_cdgf;}
      const P_cDGF & getP_cdgfOld() const {return pold_cdgf;}

      P_gDGF & getP_gdgfOld(){return pold_gdgf;}
      const P_gDGF & getP_gdgfOld() const {return pold_gdgf;}

      S_lDGF & getS_ldgfOld(){return sold_ldgf;}
      const S_lDGF & getS_ldgfOld() const {return sold_ldgf;}

      S_gDGF & getS_gdgfOld(){return sold_gdgf;}
      const S_gDGF & getS_gdgfOld() const {return sold_gdgf;}


      VELDGF & getVlOld(){return tpv.getVlOld();}
      const VELDGF & getVlOld() const {return tpv.getVlOld();}

      VELDGF & getVgOld(){return tpv.getVgOld();}
      const VELDGF & getVgOld() const {return tpv.getVgOld();}

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        if (!graphics) return;
        typedef Dune::PDELab::VTKGridFunctionAdapter<S_lDGF> VTK1;

        pvdwriter.addCellData(std::shared_ptr<VTK1>(new VTK1(s_ldgf,"s_l")));
        if (vtklevel>1)
          {
            typedef Dune::PDELab::VTKGridFunctionAdapter<PF> VTK2;
            typedef Dune::PDELab::VTKGridFunctionAdapter<P_lDGF> VTK3;
            typedef Dune::PDELab::VTKGridFunctionAdapter<P_gDGF> VTK4;
            typedef Dune::PDELab::VTKGridFunctionAdapter<P_cDGF> VTK5;

            pvdwriter.addCellData(std::shared_ptr<VTK2>(new VTK2(pf,"K")));
            pvdwriter.addCellData(std::shared_ptr<VTK3>(new VTK3(p_ldgf,"p_l")));
            pvdwriter.addCellData(std::shared_ptr<VTK4>(new VTK4(p_gdgf,"p_g")));
            pvdwriter.addCellData(std::shared_ptr<VTK5>(new VTK5(p_cdgf,"p_c")));
          }
        if (vtklevel>2)
          {

            typedef Dune::PDELab::VTKGridFunctionAdapter<P_lDGF> VTK6;
            typedef Dune::PDELab::VTKGridFunctionAdapter<P_gDGF> VTK7;
            typedef Dune::PDELab::VTKGridFunctionAdapter<P_cDGF> VTK8;
            typedef Dune::PDELab::VTKGridFunctionAdapter<S_gDGF> VTK9;
            typedef Dune::PDELab::VTKGridFunctionAdapter<S_lDGF> VTK10;
            pvdwriter.addCellData(std::shared_ptr<VTK6>(new VTK6(pold_ldgf,"p_l_old")));
            pvdwriter.addCellData(std::shared_ptr<VTK7>(new VTK7(pold_gdgf,"p_g_old")));
            pvdwriter.addCellData(std::shared_ptr<VTK8>(new VTK8(pold_cdgf,"p_c_old")));
            pvdwriter.addCellData(std::shared_ptr<VTK9>(new VTK9(sold_gdgf,"s_g_old")));
            pvdwriter.addCellData(std::shared_ptr<VTK10>(new VTK10(sold_ldgf,"s_l_old")));
          }
        //if (vtklevel>1)
        tpv.setOutput(pvdwriter,vtklevel);
      }

      void updateVelocity(RF time)
      {
        tpv.interpolate(time);
      }

      bool apply(RF time, RF dt)
      {
        try {
          pold = pnew;

          if (verbosity)
            std::cout << "======= solve two phase problem in class =======\n";
          watch.reset();
          gv.comm().barrier();
          //std::cout << "rank " << gv.comm().rank() << " time " << time << " dt " << dt << std::endl;
          osm.apply(time,dt,pold,pnew);
          // control the difference in saturation
          if(!controlDGFDifference(sold_ldgf,s_ldgf,param.sub("Setup").get<RF>("saturationchange"),0))
            {
              if (gv.comm().rank() == 0)
                std::cout << "TwoPhaseController error, change in saturation is too big" << std::endl;
              pnew = pold;
              return false;
          }

          Dune::PDELab::CopyDataHandle<TPGFS,V> copydh(tpgfs,pnew);
          if (tpgfs.gridView().comm().size()>1) // collect correction
            tpgfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

          tpv.interpolate(time+dt);
          timestepcounter++;
          if (verbosity)
            std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
          return true;
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
          if (verbosity)
            std::cout << "Newton Linesearch Error" << std::endl;
          pnew = pold;
          return false;
        }
        // newton convergence error
        catch (Dune::PDELab::NewtonNotConverged) {
          if (verbosity)
            std::cout << "Newton Convergence Error" << std::endl;
          pnew = pold;
          return false;
        }
        // linear solver error
        catch (Dune::PDELab::NewtonLinearSolverError) {
          if (verbosity)
            std::cout << "Newton Linear Solver Error" << std::endl;
          pnew = pold;
          return false;
        }
        // istl error
        catch (Dune::ISTLError) {
          if (verbosity)
            std::cout << "ISTL Error" << std::endl;
          pnew = pold;
          return false;
        }
        catch (Dune::Exception &e){
          if (verbosity)
          std::cout << "Dune reported error: " << e << std::endl;
          pnew = pold;
          return false;
        }
        catch (...){
          std::cout << "Unknown exception thrown!" << std::endl;
          return false;
        }
      }

      bool applyZeroVelocity(RF time, bool both=true)
      {
        tp.preStep(time,0.1,0);
        tpv.zeroVelocity(both);
        return true;
      }


      bool applyStationary(RF time)
      {
        tp.preStep(time,0.1,0);
        pold = pnew;
        // setup stationary problem
        typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MBE,RF,RF,RF,C,C> GO;
        GO go(tpgfs,cg,tpgfs,cg,lop,mbe);

        typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GO> LS;
        LS ls(tpgfs,1000,param.sub("Newton").get<int>("LinearVerbosity", 0));

        typedef Dune::PDELab::Newton<GO,LS,V> PDESOLVER;
        PDESOLVER newton(go,ls);
        // Newton parameters class
        typedef Dune::PDELab::NewtonParameters NewtonParameters;
        NewtonParameters newtonparameters(param.sub("Newton"));

        newtonparameters.set(newton);
        newton.setLineSearchMaxIterations(60);
        newton.setAbsoluteLimit(1.e-10);
        newton.setMinLinearReduction(1.e-3);

        try {
          if (verbosity)
            std::cout << "======= solve stationary two phase problem =======\n";
          watch.reset();
          // pnew = 0.;
          newton.apply(pnew);
          pold = pnew;

          Dune::PDELab::CopyDataHandle<TPGFS,V> copydh(tpgfs,pnew);
          if (tpgfs.gridView().comm().size()>1) // collect correction
            tpgfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
          tpv.interpolate(time+1.e-6);

          if (verbosity)
            std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
          return true;
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
          if (verbosity)
            std::cout << "Newton Linesearch Error" << std::endl;
          pnew = pold;
          return false;
        }
        // newton convergence error
        catch (Dune::PDELab::NewtonNotConverged) {
          if (verbosity)
            std::cout << "Newton Convergence Error" << std::endl;
          pnew = pold;
          return false;
        }
        // linear solver error
        catch (Dune::PDELab::NewtonLinearSolverError) {
          if (verbosity)
            std::cout << "Newton Linear Solver Error" << std::endl;
          pnew = pold;
          return false;
        }
        // istl error
        catch (Dune::ISTLError) {
          if (verbosity)
            std::cout << "ISTL Error" << std::endl;
          pnew = pold;
          return false;
        }

        catch (...){
          std::cerr << "Unknown exception thrown!" << std::endl;
          return false;
        }

      }


      void resetCounter()
      {
        timestepcounter=0;
      }

      size_t getCounter() const
      {
        return timestepcounter;
      }

      const Dune::ParameterTree & param;

    private:
      const GV& gv;
      TP &tp;
      int verbosity;
      const bool graphics;
      const int vtklevel;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      TPGFS tpgfs;
      P_lSUB p_lsub;
      P_cSUB p_csub;
      LOP lop;
      MLOP mlop;
      MBE mbe;
      C cg;
      GO0 go0;
      GO1 go1;
      IGO igo;
      V pold;
      V pnew;
      P_lDGF p_ldgf;
      P_cDGF p_cdgf;
      P_gDGF p_gdgf;
      P_lDGF pold_ldgf;
      P_cDGF pold_cdgf;
      P_gDGF pold_gdgf;
      S_lDGF s_ldgf;
      S_lDGF sold_ldgf;
      S_gDGF s_gdgf;
      S_gDGF sold_gdgf;
      const PF pf;

      TSMETHOD tsmethod;
      LS ls;
      PDESOLVER newton;
      OSM osm;

      TPV tpv;
      size_t timestepcounter;

    };





    template<typename GV, typename TP, typename VliquidDGF=ConstantVelocity<GV,double> >
    class TwoPhaseSaturatedConstantVelocity
    {

    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};
      typedef ConstADGF<GV,RF> S_lDGF; // water saturation
      typedef Dune::PDELab::PiolaBackwardAdapter<VliquidDGF> VliquidBACK; // apply inverse Piola transformation on each element
      typedef Dune::PDELab::RT0QLocalFiniteElementMap<GV,DF,RF,dim> RT0FEM; // RT0 FEM
      typedef Dune::PDELab::GridFunctionSpace<GV,RT0FEM,Dune::PDELab::NoConstraints,
                                             Dune::PDELab::istl::VectorBackend<> > RT0GFS;
      typedef typename Dune::PDELab::Backend::Vector<RT0GFS,RF>::Type RT0VEC; // vector backend
      typedef Dune::PDELab::DiscreteGridFunctionPiola<RT0GFS,RT0VEC> VELDGF; // discrete grid function for velocity



      /** \brief Constructor
          \param [in] gv_ The grid view
          \param [in] fem_ The finite element map
          \param [in] level_ The level of the grif
      */
      TwoPhaseSaturatedConstantVelocity(const GV &gv_, TP& tp_, RF velocity_x, RF velocity_y, RF velocity_z = 0.)
        : gv(gv_), tp(tp_),
          verbosity(1),
          s_ldgf(gv),
          vliquiddgf(gv,velocity_x,velocity_y, velocity_z),
          vliquiddgf_back(vliquiddgf),
          rt0fem(gv),
          rt0gfs(gv,rt0fem),
          rt0vec_l(rt0gfs,0.0),
          vldgf(rt0gfs,rt0vec_l)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;
        Dune::PDELab::interpolate(vliquiddgf_back,rt0gfs,rt0vec_l);
      }

      S_lDGF & getS_ldgf(){return s_ldgf;}
      const S_lDGF & getS_ldgf() const {return s_ldgf;}

      S_lDGF & getS_ldgfOld(){return s_ldgf;}
      const S_lDGF & getS_ldgfOld() const {return s_ldgf;}

      VELDGF & getVl(){return vldgf;}
      const VELDGF & getVl() const {return vldgf;}

      VELDGF & getVlOld(){return vldgf;}
      const VELDGF & getVlOld() const {return vldgf;}



      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<S_lDGF>(s_ldgf,"s_l"));
        pvdwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<VELDGF>(vldgf,"liquid velocity"));
      }

      bool apply(RF time, RF dt)
      {
        if (verbosity)
          std::cout << "Twophase problem was not solved (Using class TwoPhaseSaturatedConstantVelocity)\n";
        return true;
      }


      bool applyStationary(RF time)
      {
        if (verbosity)
          std::cout << "Twophase problem was not solved (Using class TwoPhaseSaturatedConstantVelocity)\n";
        return true;
      }


    private:
      const GV& gv;
      TP &tp;
      int verbosity;
      Dune::Timer watch;
      S_lDGF s_ldgf;
      VliquidDGF vliquiddgf;
      VliquidBACK vliquiddgf_back;

      RT0FEM rt0fem;
      RT0GFS rt0gfs;

      RT0VEC rt0vec_l;



      VELDGF vldgf;

    };




  } // end namespace Dycap
} // end namespace Dune

#endif
