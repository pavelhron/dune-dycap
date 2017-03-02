// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TRANSPORTSOLVER_BOILERPLATE_HH
#define DUNE_DYCAP_TRANSPORTSOLVER_BOILERPLATE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<dune/common/exceptions.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>

#include"../oldpdelab/pdelab.hh"

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
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include"../models/limiter.hh"
#include<src/output/output.hh>
#include"../utilities/time_step_manager.hh"
#include"../models/componenttransportop.hh"
#include"../parameters/transportparameters.hh"
#include"../utilities/cfltransportcontroller.hh"




namespace Dune {
  namespace Dycap{


    template<bool b, typename T, typename U>
    struct select
    {
      typedef T type;
    };

    template<typename T, typename U>
    struct select<false, T, U>
    {
      typedef U type;
    };

    template<int b, typename T, typename U>
    struct selecti
    {
      typedef T type;
    };

    template<typename T, typename U>
    struct selecti<0, T, U>
    {
      typedef U type;
    };


    template<typename GV, typename CTP, typename ODE, bool reconstruction = false>
    class ExplicitTransportReactionSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      // grid operator types
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type CV;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,CV> DGF;

      typedef Dune::PDELab::LimiterFV<typename CTP::Traits, GFS> Limiter;

      typedef Dune::PDELab::DefaultFluxReconstruction<typename CTP::Traits> DFR;

      typedef typename select<reconstruction,Limiter,DFR>::type FR;
      typedef FluxAdapter<DGF,FR> ConcDGF;
      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP,FR> LOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;
      typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
      typedef Dune::PDELab::SimpleTimeController<RF> TC;
      typedef Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,CV,CV,TC> OSM;

      typedef ConcInitial<GV,RF> InitialConcentration;

      typedef CFLTransportController<GV,CTP,FR> CFLController;
      typedef Dune::PDELab::ExplicitEulerParameter<RF> TimeSteppingMethod;



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

      template< class U>
      typename std::enable_if<std::is_same<U, Limiter>::value, shared_ptr<U> >::type getFlux(DGF& cdgf, GFS& gfs_)
      {
        return shared_ptr<FR>(new FR(gfs_,param.get<RF>("Timeloop.theta"),param.get<std::string>("Timeloop.limiter")));
      }

      template< class U>
      typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(DGF& cdgf, GFS& gfs_)
      {
        return shared_ptr<FR>(new FR());
      }

      ExplicitTransportReactionSolver(const GV &gv_, CTP& ctp_, ODE& ode_, const Dune::ParameterTree & param_, const std::string cname)
        : param(param_), gv(gv_), ctp(ctp_), ode(ode_),
          verbosity(1),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cnull(gfs,0.),
          dgf(gfs,cold),
          flux(getFlux<FR>(dgf,gfs)),
          concdgf(dgf,*flux),
          lop(ctp,*flux),
          cflcontroller(gv,ctp,*flux),
          mlop(ctp),
          cg(getContainer()),
          mbe(5),
          go0(gfs,cg,gfs,cg,lop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo(go0,go1),
          tc(),
          ls(gfs),
          timesteppingmethod(),
          osm(timesteppingmethod,igo,ls,tc)

      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        igo.setMethod(timesteppingmethod);
        timesteppingmethods.setTimestepMethod(osm,param.sub("Timeloop").get<std::string>("explicit"));
        osm.setVerbosityLevel(2);
        InitialConcentration conc_initial(gv,param,cname);
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::constraints(conc_initial,gfs,cg,false);

        ctp.setTimeTarget(1.e100,1.e100);


      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"concentration boilerplate"));
      }

      ConcDGF & getConcDGF(){flux->prestage(cold); return concdgf;}
      //const ConcDGF & getConcDGF() const {return concdgf;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      void setSolution(V & c) {cnew = c; cold = cnew;}


      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF timetarget = time+dt;
        ctp.setTimeTarget(timetarget,dt);
        RF timestep = std::min(cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl"),dt);

        try {

          while (time<timetarget-1.e-8)
            {

              if (verbosity)
                std::cout << "======= solve transport in boilerplate class from time " << time << " to time " << time+timestep << " =======\n";
              watch.reset();
              if (!reconstruction)
                osm.apply(time,timestep,cold,cnew);
              else
                osm.apply(time,timestep,cold,cnew,*flux);

              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s" << std::endl;

              //std::vector<V*> datavector;
              //datavector.push_back(&cnew);
              std::cout << "reaction " << std::endl;
              ode.apply(time,timestep,cnew);
              std::cout << "reaction finished" << std::endl;

              cold = cnew;
              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
            }
          cold = cnew;
          return true;
        }

        // istl error
        catch (Dune::ISTLError) {
          if (verbosity)
            std::cout << "ISTL Error" << std::endl;
          cnew = cold;
          return false;
        }
        catch (Dune::Exception &e){
          std::cerr << "Dune reported error: " << e << std::endl;
          return false;
        }

        catch (...){
          std::cerr << "Unknown exception thrown!" << std::endl;
          return false;
        }

      }

      const RF suggestTimeStep()
      {
        return cflcontroller.suggestTimeStep();
      }

      const Dune::ParameterTree & param;

    private:
      GV gv;
      CTP &ctp;
      ODE &ode;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cnull;
      DGF dgf;
      shared_ptr<FR> flux;
      ConcDGF concdgf;
      LOP lop;
      CFLController cflcontroller;
      MLOP mlop;
      C cg;
      MBE mbe;
      GO0 go0;
      GO1 go1;
      IGO igo;
      TC tc;
      LS ls;
      TimeSteppingMethod timesteppingmethod;
      OSM osm;
      Dune::PDELab::TimeSteppingMethods<RF> timesteppingmethods;
      //Limiter limiter;
    };

    template<typename GV, typename CTP, bool reconstruction = false, typename Initial = ConcInitial<GV,double> >
    class ExplicitTransportSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      // grid operator types
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type CV;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,CV> DGF;

      typedef Dune::PDELab::DefaultFluxReconstruction<typename CTP::Traits> DFR;
      typedef Dune::PDELab::LimiterFV<typename CTP::Traits, GFS> Limiter;

      typedef typename select<reconstruction,Limiter,DFR>::type FR;
      typedef FluxAdapter<DGF,FR> ConcDGF;

      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP,FR> LOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
      typedef typename IGO::Traits::Domain V;


      typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
      typedef Dune::PDELab::SimpleTimeController<RF> TC;
      typedef Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,CV,CV,TC> OSM;

      typedef Initial InitialConcentration;

      typedef CFLTransportController<GV,CTP,FR> CFLController;
      typedef Dune::PDELab::ExplicitEulerParameter<RF> TimeSteppingMethod;


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

      /*
        template <class T>
        class FluxFunction
        {
        public:
        template< class U = T>
        typename std::enable_if<std::is_same<U, SOFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
        {
        return shared_ptr<FR>(new FR(cdgf,2.));
        }

        template< class U = T>
        typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
        {
        return shared_ptr<FR>(new FR());
        }

        };
      */

      template< class U>
      typename std::enable_if<std::is_same<U, Limiter>::value, shared_ptr<U> >::type getFlux(GFS& gfs_)
      {
        return shared_ptr<FR>(new FR(gfs_,param.get<RF>("Timeloop.theta"),param.get<std::string>("Timeloop.limiter")));
      }

      template< class U>
      typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(GFS& gfs_)
      {
        return shared_ptr<FR>(new FR());
      }


      ExplicitTransportSolver(const GV &gv_, CTP& ctp_, const Dune::ParameterTree & param_)
        : param(param_), gv(gv_), ctp(ctp_),
          verbosity(1),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cnull(gfs,0.),
          dgf(gfs,cnew),
          flux(getFlux<FR>(gfs)),
          concdgf(dgf,*flux),
          lop(ctp,*flux),
          cflcontroller(gv,ctp,flux),
          mlop(ctp),
          mbe(5),
          cg(getContainer()),
          go0(gfs,cg,gfs,cg,lop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo(go0,go1),
          tc(),
          ls(gfs),
          timesteppingmethod(),
          osm(timesteppingmethod,igo,ls,tc)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        igo.setMethod(timesteppingmethod);
        osm.setVerbosityLevel(2);
        InitialConcentration conc_initial(gv,param,ctp.getName());
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::constraints(conc_initial,gfs,cg,false);

      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"concentration boilerplate"));
      }

      ConcDGF & getConcDGF(){return concdgf;}
      const ConcDGF & getConcDGF() const {return concdgf;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      void setSolution(V & c) {cnew = c; cold = c;}
      void nullSolution() {cnew = cold = 0.0;}


      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF timestep = std::min(cflcontroller.suggestTimeStep(),dt);
        RF timetarget = time+dt;

        try {

          while (time<timetarget-1.e-8)
            {
              ctp.setTime(time);
              if (verbosity)
                std::cout << "======= solve "<< ctp.getName() << "transport from time " << time << " to time " << time+timestep << " =======\n";
              watch.reset();
              if (!reconstruction)
                osm.apply(time,timestep,cold,cnew);
              else
                osm.apply(time,timestep,cold,cnew,*flux);
              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
              cold = cnew;
              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
            }
          cold = cnew;
          return true;
        }

        // istl error
        catch (Dune::ISTLError) {
          if (verbosity)
            std::cout << "ISTL Error" << std::endl;
          cnew = cold;
          return false;
        }

        catch (...){
          std::cerr << "Unknown exception thrown!" << std::endl;
          return false;
        }

      }

      shared_ptr<FR> getFlux()
      {
        return flux;
      }

      const Dune::ParameterTree & param;

    private:
      GV gv;
      CTP &ctp;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cnull;
      DGF dgf;
      shared_ptr<FR> flux;
      ConcDGF concdgf;
      LOP lop;
      CFLController cflcontroller;
      MLOP mlop;
      MBE mbe;
      C cg;
      GO0 go0;
      GO1 go1;
      IGO igo;
      TC tc;
      LS ls;
      TimeSteppingMethod timesteppingmethod;
      OSM osm;


    };

    /*
    template<typename GV, typename CTP>
    class ImplicitTransportSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      // grid operator types
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type CV;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,CV> ConcDGF;

      typedef Dune::PDELab::ImplicitCCFVSpatialTransportOperator<CTP> LOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,true> IGO;
      typedef typename IGO::Traits::Domain V;

      // choose linear solver
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> LS;
      typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,V> CPDESOLVER;
      typedef Dune::PDELab::ImplicitEulerParameter<RF> TimeSteppingMethod;
      typedef Dune::PDELab::TimeSteppingMethods<RF> TimeSteppingMethods;

      typedef Dune::PDELab::OneStepMethod<RF,IGO,CPDESOLVER,V,V> OSM;
      typedef ConcInitial<GV,RF> InitialConcentration;

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



      ImplicitTransportSolver(GV &gv_, CTP& ctp_, const Dune::ParameterTree & param_)
        : param(param_), gv(gv_), ctp(ctp_),
          verbosity(1),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cnull(gfs,0.),
          concdgf(gfs,cnew),
          lop(ctp,param.sub("Setup").get<bool>("central")),
          mlop(ctp),
          mbe(5),
          cg(getContainer()),
          go0(gfs,cg,gfs,cg,lop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo(go0,go1),
          ls(gfs,cg,5000,0), //superLU
          timesteppingmethod(),
          timesteppingmethods(),
          cpdesolver(igo,ls,1e-8,1e-14),
          osm(timesteppingmethod,igo,cpdesolver)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        timesteppingmethods.setTimestepMethod(osm,param.sub("Timeloop").get<std::string>("implicit"));
        osm.setVerbosityLevel(param.sub("Verbosity").get<int>("ImplicitTransport", 0));
        InitialConcentration conc_initial(gv,param,ctp.getName());
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::constraints(conc_initial,gfs,cg,false);

      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"concentration boilerplate"));
      }

      ConcDGF & getConcDGF(){return concdgf;}
      const ConcDGF & getConcDGF() const {return concdgf;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      void setSolution(V & c) {cnew = c; cold = c;}
      void nullSolution() {cnew = cold = 0.0;}

      void setTimeSteppingMethod(const std::string method_)
      {
        timesteppingmethods.setTimestepMethod(osm,method_);
      }

      const GFS & getGFS() const {return gfs;}

      const FEM & getFem() const {return fem;}

      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF timestep = dt;
        RF timetarget = time+dt;
        ctp.setTimeTarget(timetarget,dt);

        while (time<timetarget-1.e-8)
          {
            ctp.setTime(time);
            // oxygen transport in liquid phase

            if (gv.comm().rank() == 0)
              std::cout << "======= composite reaction =======\n";
            watch.reset();
            try {

              timestep = osm.apply(time,timestep,cold,cnew);
              if (gv.comm().rank() == 0)
                std::cout << "======= composite reaction finished ======= time : " << watch.elapsed() << " s" << std::endl;

              if (!controlReactionTimeStep(gv,cnew))
                {
                  cnew = cold;
                  timestep/=2.;
                  std::cout << "timestep was reduced to " << timestep << std::endl;
                  continue;
                }

              cold = cnew; // liquid
              // update  time step
              time += timestep;
              timestep*=2.;
              timestep = std::min(timestep, timetarget-time);

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
              std::cout << "Unknown exception thrown!" << std::endl;
              timestep*=0.5;
              cnew=cold;
              continue;
            }


          }



        return true;
      }

      const Dune::ParameterTree & param;

    private:
      GV gv;
      CTP &ctp;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cnull;
      ConcDGF concdgf;
      LOP lop;
      MLOP mlop;
      MBE mbe;
      C cg;
      GO0 go0;
      GO1 go1;
      IGO igo;
      LS ls;
      TimeSteppingMethod timesteppingmethod;
      TimeSteppingMethods timesteppingmethods;
      CPDESOLVER cpdesolver;
      OSM osm;
      };
    */

    template<typename GV, typename CTP, typename RA>
    class ImplicitTransportReactionSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      // grid operator types
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type CV;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,CV> ConcDGF;

      typedef Dune::PDELab::ImplicitCCFVSpatialTransportOperator<CTP,RA> LOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,true> IGO;
      typedef typename IGO::Traits::Domain V;

      // choose linear solver
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> LS;
      typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,V> CPDESOLVER;
      typedef Dune::PDELab::ImplicitEulerParameter<RF> TimeSteppingMethod;

      typedef Dune::PDELab::OneStepMethod<RF,IGO,CPDESOLVER,V,V> OSM;
      typedef ConcInitial<GV,RF> InitialConcentration;

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



      ImplicitTransportReactionSolver(const GV &gv_, CTP& ctp_, RA& ra, const Dune::ParameterTree & param_, const std::string cname)
        : param(param_), gv(gv_), ctp(ctp_),
          verbosity(1),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cnull(gfs,0.),
          concdgf(gfs,cold),
          lop(ctp,ra,param.sub("Setup").get<bool>("central")),
          mlop(ctp),
          cg(getContainer()),
          mbe(5),
          go0(gfs,cg,gfs,cg,lop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo(go0,go1),
          ls(gfs,cg,5000,0), //superLU
          timesteppingmethod(),
          cpdesolver(igo,ls,1e-8,1e-14),
          osm(timesteppingmethod,igo,cpdesolver)

      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        igo.setMethod(timesteppingmethod);
        osm.setVerbosityLevel(2);
        InitialConcentration conc_initial(gv,param,cname);
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::constraints(conc_initial,gfs,cg,false);

      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"concentration boilerplate"));
      }

      ConcDGF & getConcDGF(){return concdgf;}
      const ConcDGF & getConcDGF() const {return concdgf;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      void setSolution(V & c) {cnew = c; cold = c;}
      void nullSolution() {cnew = cold = 0.0;}

      const GFS & getGFS() const {return gfs;}

      const FEM & getFem() const {return fem;}

      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF timestep = dt;
        RF timetarget = time+dt;
        ctp.setTimeTarget(timetarget,dt);

        try {

          while (time<timetarget-1.e-12)
            {

              if (verbosity)
                std::cout << "======= solve transport in boilerplate class from time " << time << " to time " << time+timestep << " =======\n";
              watch.reset();
              osm.apply(time,timestep,cold,cnew);

              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
              cold = cnew;
              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
            }
          cold = cnew;
          return true;

        }

        // istl error
        catch (Dune::ISTLError) {
          if (verbosity)
            std::cout << "ISTL Error" << std::endl;
          cnew = cold;
          return false;
        }

        catch (Dune::Exception &e){
          std::cerr << "Dune reported error: " << e << std::endl;
          return false;
        }

        catch (...){
          std::cerr << "Unknown exception thrown!" << std::endl;
          return false;
        }


      }

      const Dune::ParameterTree & param;

    private:
      const GV& gv;
      CTP &ctp;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cnull;
      ConcDGF concdgf;
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
    };


    template<typename GV, typename CTP, bool reconstruction = false, typename Initial = ConcInitial<GV,double> >
    class TransportSolver
    {
    public:
      // choose some types
      typedef typename GV::Grid::ctype DF;
      typedef double RF;

      enum { dim = GV::dimension};

      // common grid function types
      typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
      typedef Dune::PDELab::istl::VectorBackend<> VBE;
      typedef Dune::PDELab::P0ParallelConstraints CON;
      typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;

      // grid operator types
      typedef typename GFS::template ConstraintsContainer<RF>::Type C;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type CV;
      typedef Dune::PDELab::DiscreteGridFunction<GFS,CV> DGF;

      typedef Dune::PDELab::DefaultFluxReconstruction<typename CTP::Traits> DFR;
      typedef Dune::PDELab::LimiterFV<typename CTP::Traits, GFS> Limiter;

      typedef typename select<reconstruction,Limiter,DFR>::type FR;
      typedef DGF ConcDGF;

      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP,FR> LOP;
      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP> ILOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::GridOperator<GFS,GFS,ILOP,MBE,RF,RF,RF,C,C> IGO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> IGO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
      typedef Dune::PDELab::OneStepGridOperator<IGO0,IGO1,true> IIGO;
      typedef typename IGO::Traits::Domain V;


      typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
      // typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> ILS;
      // typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IIGO> ILS;
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,C> ILS;
      //typedef Dune::PDELab::StationaryLinearProblemSolver<IIGO,ILS,V> CPDESOLVER;
      typedef Dune::PDELab::Newton<IIGO,ILS,V> CPDESOLVER;


      typedef Dune::PDELab::SimpleTimeController<RF> TC;
      typedef Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,CV,CV,TC> OSM;

      typedef Initial InitialConcentration;

      typedef CFLTransportController<GV,CTP,FR> CFLController;
      typedef Dune::PDELab::ExplicitEulerParameter<RF> TimeSteppingMethod;
      typedef Dune::PDELab::ImplicitEulerParameter<RF> ITimeSteppingMethod;
      typedef Dune::PDELab::OneStepMethod<RF,IIGO,CPDESOLVER,V,V> IOSM;


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

      /*
        template <class T>
        class FluxFunction
        {
        public:
        template< class U = T>
        typename std::enable_if<std::is_same<U, SOFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
        {
        return shared_ptr<FR>(new FR(cdgf,2.));
        }

        template< class U = T>
        typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
        {
        return shared_ptr<FR>(new FR());
        }

        };
      */

      template< class U>
      typename std::enable_if<std::is_same<U, Limiter>::value, shared_ptr<U> >::type getFlux(GFS& gfs_)
      {
        return shared_ptr<FR>(new FR(gfs_,param.get<RF>("Timeloop.theta"),param.get<std::string>("Timeloop.limiter")));
      }

      template< class U>
      typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(GFS& gfs_)
      {
        return shared_ptr<FR>(new FR());
      }


      TransportSolver(const GV &gv_, CTP& ctp_, const Dune::ParameterTree & param_, const  bool implicit_=false)
        : param(param_), gv(gv_), ctp(ctp_),
          verbosity(1),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cnull(gfs,0.),
          dgf(gfs,cnew),
          flux(getFlux<FR>(gfs)),
          concdgf(gfs,cnew),
          concdgfold(gfs,cold),
          lop(ctp,*flux),
          ilop(ctp,param.sub("Setup").get<RF>("sgmin")),
          cflcontroller(gv,ctp,flux,1.e-7,1),
          mlop(ctp),
          imlop(ctp,param.sub("Setup").get<RF>("sgmin")),
          mbe(5),
          cg(getContainer()),
          go0(gfs,cg,gfs,cg,lop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo0(gfs,cg,gfs,cg,ilop,mbe),
          igo1(gfs,cg,gfs,cg,imlop,mbe),
          igo(go0,go1),
          iigo(igo0,igo1),
          tc(),
          ls(gfs),
          //ils(gfs,cg,5000,0), //superLU
          ils(gfs,cg,1000,5,param.sub("Verbosity").get<int>("LinearSolver")),
          timesteppingmethod(),
          itimesteppingmethod(),
          osm(timesteppingmethod,igo,ls,tc),
          // cpdesolver(iigo,ils,1e-8,1e-14),
          cpdesolver(iigo,ils),
          iosm(itimesteppingmethod,iigo,cpdesolver),
          implicit(implicit_),
          timestepcounter(0)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        if (implicit)
          timesteppingmethods.setTimestepMethod(iosm,param.sub("Timeloop").get<std::string>("implicit"));
        else
          timesteppingmethods.setTimestepMethod(osm,param.sub("Timeloop").get<std::string>("explicit"));
        osm.setVerbosityLevel(2);
        iosm.setVerbosityLevel(2);
        InitialConcentration conc_initial(gv,param,ctp.getName());
        Dune::PDELab::constraints(conc_initial,gfs,cg);
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        cnew=cold;

        //std::cout << "concentration name for initial " << ctp.getName() << std::endl;

        //for (auto i:cold)
        //  std::cout << i <<" ";
        //std::cout << std::endl;

        flux->prestage(cnew);

        typedef Dune::PDELab::NewtonParameters NewtonParameters;
        NewtonParameters newtonparameters(param.sub("Newton"));
        // set new parameters
        newtonparameters.set(cpdesolver);
      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,ctp.getName()));
      }

      ConcDGF & getConcDGF(){return concdgf;}
      const ConcDGF & getConcDGF() const {return concdgf;}

      ConcDGF & getOldConcDGF(){return concdgfold;}
      const ConcDGF & getOldConcDGF() const {return concdgfold;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      GFS & getGFS() {return gfs;}
      const GFS & getGFS() const {return gfs;}

      void setSolution(V & c) {cnew = c; cold = c;}
      void nullSolution() {cnew = cold = 0.0;}

      void update()
      {
        if (!implicit)
          flux->prestage(cnew);
        Dune::PDELab::CopyDataHandle<GFS,V> copydh(gfs,cnew);
        if (gfs.gridView().comm().size()>1) // collect correction
          gfs.gridView().communicate(copydh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
      }


      bool apply(RF time, RF dt)
      {
        update();
        cold=cnew;

        RF timestep = dt;
        if (!implicit)
          timestep = std::min(cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl"),dt);
        RF timetarget = time+dt;



        while (time<timetarget-1.e-8)
          {
            try {

              ctp.setTime(time);
              if (implicit)
                {
                  if (verbosity)
                    std::cout << "======= solve implicit "<< ctp.getName() << " transport from time " << time << " to time " << time+timestep << " =======\n";
                   watch.reset();
                   /*for (auto i:cnew)
                    std::cout << i <<"\n";

                  for (auto i:cold)
                    std::cout << i <<"\n";

                  iosm.apply(time,timestep,cold,cnew);
                  V v(gfs,0.);
                  ilop.setVerbosity(3);
                  igo0.residual(cnew,v);
                  ilop.setVerbosity(0);
                   */
                   iosm.apply(time,timestep,cold,cnew);
                }
              else {
                if (verbosity)
                  std::cout << "======= solve explicit "<< ctp.getName() << " transport from time " << time << " to time " << time+timestep << " =======\n";
                watch.reset();
                if (!reconstruction)
                  osm.apply(time,timestep,cold,cnew);
                else
                  osm.apply(time,timestep,cold,cnew,*flux);
              }
              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s" << std::endl;

              update();
              cold=cnew;
              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
              timestepcounter++;
            }
            // istl error
            catch (Dune::ISTLError) {
              if (verbosity)
                std::cout << "ISTL Error" << std::endl;
              cnew = cold;
              timestep*=0.5;
            }

            catch (Dune::Exception &e){
              if (verbosity)
                std::cout << "Dune reported error: " << e << std::endl;
              cnew = cold;
              timestep*=0.5;
            }
            catch (...){
              if (verbosity)
                std::cout << "Unknown exception thrown!" << std::endl;
              cnew = cold;
              timestep*=0.5;
            }
          }
        update();
        cold=cnew;
        return true;


      }

      shared_ptr<FR> getFlux()
      {
        return flux;
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
      GV gv;
      CTP &ctp;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cnull;
      DGF dgf;
      shared_ptr<FR> flux;
      ConcDGF concdgf;
      ConcDGF concdgfold;
      LOP lop;
      ILOP ilop;
      CFLController cflcontroller;
      MLOP mlop;
      MLOP imlop;
      MBE mbe;
      C cg;
      GO0 go0;
      GO1 go1;
      IGO0 igo0;
      IGO1 igo1;
      IGO igo;
      IIGO iigo;
      TC tc;
      LS ls;
      ILS ils;
      TimeSteppingMethod timesteppingmethod;
      ITimeSteppingMethod itimesteppingmethod;
      OSM osm;
      CPDESOLVER cpdesolver;
      IOSM iosm;
      Dune::PDELab::TimeSteppingMethods<RF> timesteppingmethods;
      bool implicit;
      size_t timestepcounter;
    };




  } // end namespace Dycap
} // end namespace Dune

#endif
