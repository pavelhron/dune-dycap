// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_MULTICOMPONENTTRANSPORTSOLVER_BOILERPLATE_HH
#define DUNE_DYCAP_MULTICOMPONENTTRANSPORTSOLVER_BOILERPLATE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/backend/istl.hh>

#include<src/newton/newton_reaction.hh>
#include<src/newton/newton_utilities.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include<src/output/output.hh>
#include"../utilities/utilities.hh"
#include"../models/multicomponenttransportop.hh"
#include"../parameters/multitransportparameters.hh"
#include"../utilities/cfltransportcontroller.hh"



namespace Dune {
  namespace Dycap{

    template<typename RF, typename GV, typename CTP, typename RA>
    class MulticomponentImplicitTransportSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      enum { dim = GV::dimension};

      typedef GV GVT;
      typedef CTP TP;
      typedef RF RFT;


      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE0;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;

      // pgfs with blockwise mapper
      typedef Dune::PDELab::ISTLVectorBackend
      <Dune::PDELab::ISTLParameters::static_blocking,CTP::COMPONENTS> VBE;
      typedef Dune::PDELab::PowerGridFunctionSpace<GFS,CTP::COMPONENTS,VBE,
                                                   Dune::PDELab::EntityBlockedOrderingTag> PGFS;
      // grid operator types
      typedef typename PGFS::template ConstraintsContainer<RF>::Type C;



      typedef Dune::PDELab::MulticomponentCCFVSpatialTransportOperator<CTP,RA> LOP;
      typedef Dune::PDELab::MulticomponentCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<PGFS,PGFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<PGFS,PGFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,true> IGO;
      //typedef typename IGO::Traits::Domain V;
      typedef typename Dune::PDELab::BackendVectorSelector<PGFS,RF>::Type V;
      typedef DiscreteGridFunctionFromPGF<PGFS,V,CTP::COMPONENTS> DGF;

      typedef VectorDiscreteGridFunction<PGFS,V> VDGF;
      // choose linear solver
      //typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<PGFS,C> LS;
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<PGFS,C> LS;
      typedef Dune::PDELab::NewtonReaction<IGO,LS,V> CPDESOLVER;

      typedef Dune::PDELab::ImplicitEulerParameter<RF> TimeSteppingMethod;

      typedef Dune::PDELab::OneStepMethod<RF,IGO,CPDESOLVER,V,V> OSM;
      typedef ConcInitial<GV,RF> InitialConcentration;
      typedef CFLTransportController<GV,typename CTP::ComponentType> CFLController;

      static C getContainer() {
        C cg;
        cg.clear();
        return cg;
      }

      static GeometryType getCube() {
        GeometryType gt;
        gt.makeCube(dim);
        return gt;
      }



      MulticomponentImplicitTransportSolver(const GV &gv_, CTP& ctp_, RA& ra, const Dune::ParameterTree & param_, const int verbosity_=1, const int cflcomponent = 0)
        : param(param_), gv(gv_), ctp(ctp_),
          verbosity(verbosity_),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          pgfs(gfs),
          cnew(pgfs,0.),
          cold(pgfs,0.),
          cnull(pgfs,0.),
          lop(ctp,ra,1e-7),
          mlop(ctp),
          cg(getContainer()),
          mbe(5),
          go0(pgfs,cg,pgfs,cg,lop,mbe),
          go1(pgfs,cg,pgfs,cg,mlop,mbe),
          igo(go0,go1),
          //ls(pgfs,cg,5000,0), //superLU
          ls(pgfs,cg,1000,5,param.sub("Verbosity").get<int>("ReactionLinearSolver")),
          timesteppingmethod(),
          cpdesolver(igo,ls),
          osm(timesteppingmethod,igo,cpdesolver),
          cflcontroller(gv,ctp.component(cflcomponent),1.e-7,1),
          timestepcounter(0)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        for (size_t i = 0;i<CTP::COMPONENTS;i++)
          {
            ppdgf.push_back(std::make_shared<DGF>(pgfs,cnew,i));
          }

        //ppdgf.push_back(std::shared_ptr<DGF>(new DGF(pgfs,cnew,0)));

        igo.setMethod(timesteppingmethod);
        osm.setVerbosityLevel(verbosity);


        typedef ConcInitial<GV,RF> IC;
        IC conc_initial(gv,param,ctp.getName(0));
        typedef Dune::PDELab::PowerGridFunction<IC,CTP::COMPONENTS> CType;
        CType c_initial(conc_initial);
        for (size_t i=0;i<CTP::COMPONENTS;++i)
          {
            IC conc_initial(gv,param,ctp.getName(i));
            c_initial.setChild(i,std::make_shared<IC>(conc_initial));
          }

        Dune::PDELab::constraints(c_initial,pgfs,cg);

        if (verbosity)
          std::cout << "constrained dofs=" << cg.size() << " of " << pgfs.globalSize() << std::endl;

        Dune::PDELab::interpolate(c_initial,pgfs,cold);
        cnew = cold;

        typedef Dune::PDELab::NewtonParameters NewtonParameters;
        NewtonParameters newtonparameters(param.sub("NewtonImplicit"));
        // set new parameters
        newtonparameters.set(cpdesolver);
        gfs.update();
      }

      void resetSolution()
      {
                typedef ConcInitial<GV,RF> IC;
        IC conc_initial(gv,param,ctp.getName(0));
        typedef Dune::PDELab::PowerGridFunction<IC,CTP::COMPONENTS> CType;
        CType c_initial(conc_initial);
        for (size_t i=0;i<CTP::COMPONENTS;++i)
          {
            IC conc_initial(gv,param,ctp.getName(i));
            c_initial.setChild(i,std::make_shared<IC>(conc_initial));
          }

        Dune::PDELab::constraints(c_initial,pgfs,cg);

        if (verbosity)
          std::cout << "constrained dofs=" << cg.size() << " of " << pgfs.globalSize() << std::endl;

        Dune::PDELab::interpolate(c_initial,pgfs,cold);
        cnew = cold;
        //  Dune::PDELab::interpolate(conc_initial,gfs,cold);
        // cnew = cold;
        ctp.setTimeTarget(1.e100,1.e100);
      }



      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        VTKGridFunctionImp::FillVTKFunctions<PGFS,V,PVD,PGFS::CHILDREN>::fill("",pgfs,cnew,pvdwriter);
      }

      DGF & getConcDGF(const size_t k){assert(k<PGFS::CHILDREN);return *ppdgf[k];}
      const DGF & getConcDGF(const size_t k) const {assert(k<PGFS::CHILDREN);return *ppdgf[k];}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      V & getOldSolution() {return cold;}
      const V & getOldSolution() const {return cold;}


      void setSolution(V & c) {cnew = c; cold = c;}
      void nullSolution() {cnew = cold = 0.0;}

      GFS & getGFS() {return gfs;}
      const GFS & getGFS() const {return gfs;}

      const FEM & getFem() const {return fem;}

      CTP& getCTP() {return ctp;}
      const CTP& getCTP() const {return ctp;}

      const int getGFSSize(){return gfs.globalSize();}

      inline const GV getGridView () const
      {
        return gv;
      }


      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF cfltimestep = param.get<RF>("Timeloop.implicitcfl")*cflcontroller.suggestTimeStep();
        RF timestep = std::min(dt,cfltimestep);
        RF timetarget = time+dt;
        ctp.setTimeTarget(timetarget,dt);

        while (time<timetarget-1.e-8)
          {
            ctp.setTime(time);

            if (gv.comm().rank() == 0)
              std::cout << "======= composite multicomponent ======= from " << time << " to " << time+timestep << "\n";
            watch.reset();
              try {

              timestep = osm.apply(time,timestep,cold,cnew);
              if (gv.comm().rank() == 0)
                std::cout << "======= composite multicomponent finished ======= time : " << watch.elapsed() << " s" << std::endl;

              if (!controlReactionTimeStep(gv,cnew,1))
                {
                  cnew = cold;
                  timestep/=2.;
                  std::cout << "composite multicomponent timestep was reduced to " << timestep << std::endl;
                  continue;
                }

              Dune::PDELab::CopyDataHandle<PGFS,V> copydh(pgfs,cnew);
              if (pgfs.gridView().comm().size()>1) // collect correction
                pgfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
              cold = cnew; // liquid
              // update  time step
              time += timestep;
              timestep*=1.5;
              timestep = std::min(timestep, timetarget-time);
              timestep = std::min(timestep,cfltimestep);
              ++timestepcounter;
              }
              // newton linear search error
            catch (Dune::PDELab::NewtonLineSearchError) {
              if (verbosity)
                std::cout << "Newton Linesearch Error" << std::endl;
              timestep*=0.5;
              cnew=cold;
              continue;
            }
            // newton convergence error
            catch (Dune::PDELab::NewtonNotConverged) {
              if (verbosity)
                std::cout << "Newton Convergence Error" << std::endl;
              timestep*=0.5;
              cnew=cold;
              continue;
            }
            // linear solver error
            catch (Dune::PDELab::NewtonLinearSolverError) {
              if (verbosity)
                std::cout << "Newton Linear Solver Error" << std::endl;
              timestep*=0.5;
              cnew=cold;
              continue;
            }


            catch (...){
              if (verbosity)
                std::cout << "Unknown exception thrown!" << std::endl;
              timestep*=0.5;
              cnew=cold;
              continue;
            }
          }

        return true;
      }

      bool apply_simply(RF time, RF dt)
      {
        try {

          cold = cnew;
          ctp.setTimeTarget(time+dt,dt);
          ctp.setTime(time);
          if (gv.comm().rank() == 0)
              std::cout << "======= composite multicomponent ======= from " << time << " to " << time+dt << "\n";
          watch.reset();

          osm.apply(time,dt,cold,cnew);
          if (gv.comm().rank() == 0)
            std::cout << "======= composite reaction finished ======= time : " << watch.elapsed() << " s" << std::endl;

          if (!controlReactionTimeStep(gv,cnew))
            {
              cnew = cold;
              return false;
            }

          cold = cnew; // liquid
          return true;
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
          if (verbosity)
            std::cout << "Newton Linesearch Error" << std::endl;
          cnew=cold;
          return false;
        }
        // newton convergence error
        catch (Dune::PDELab::NewtonNotConverged) {
          if (verbosity)
            std::cout << "Newton Convergence Error" << std::endl;
          cnew=cold;
          return false;
        }
        // linear solver error
        catch (Dune::PDELab::NewtonLinearSolverError) {
          if (verbosity)
            std::cout << "Newton Linear Solver Error" << std::endl;
          cnew=cold;
          return false;
        }

        // linear solver error
        catch (...) {
          if (verbosity)
            std::cout << "Error" << std::endl;
          cnew=cold;
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

      OSM & getOSM(){return osm;}

      const Dune::ParameterTree & param;

    private:
      GV gv;
      CTP &ctp;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      PGFS pgfs;
      V cnew;
      V cold;
      V cnull;
      std::vector<std::shared_ptr<DGF> > ppdgf;
      LOP lop;
      MLOP mlop;
      C cg;
      MBE mbe;
      GO0 go0;
      GO1 go1;
      IGO igo;
      LS ls;
      TimeSteppingMethod timesteppingmethod;
      CPDESOLVER cpdesolver;
      OSM osm;
      CFLController cflcontroller;
      size_t timestepcounter;
    };

  } // end namespace Dycap
} // end namespace Dune

#endif
