// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:

#ifndef DUNE_DYCAP_ADHESION_ESTIMATION_MODEL_HH
#define DUNE_DYCAP_ADHESION_ESTIMATION_MODEL_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>


#include<src/newton/newton.hh>
#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include<src/output/pvdwriter.hh>
#include<src/utilities/utilities.hh>
#include<src/ode/ode.hh>
#include<src/ode/odegridoperator.hh>
#include<src/utilities/multispantwophase.hh>

#include<src/models/componenttransportop.hh>
#include"transportparameters.hh"
#include<src/utilities/cfltransportcontroller.hh>

#include<src/utilities/utilities.hh>
#include<src/utilities/multispantwophase.hh>
#include<src/output/output.hh>

#include<src/models/componenttransportop.hh>
#include<src/boilerplate/twophasesolver_pcpl.hh>
#include<src/boilerplate/twophasevelocity.hh>
#include"reaction_model.hh"
#include<src/output/basicio.hh>



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


    // class for solving advection-diffusion-reaction problem
    template<typename GV, typename CTP, typename ODE, bool reconstruction = false>
    class CombineTransportSolver
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

      typedef Dune::PDELab::DefaultFluxReconstruction<typename CTP::Traits> DFR;
      typedef Dune::PDELab::SecondOrderFluxReconstruction<typename CTP::Traits, ConcDGF> SOFR;

      typedef typename select<reconstruction,SOFR,DFR>::type FR;

      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP,FR> LOP;
      typedef Dune::PDELab::ModifiedCCFVSpatialTransportOperator<CTP> DLOP;
      typedef Dune::PDELab::ModifiedCCFVTemporalOperator<CTP> MLOP;

      typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;

      typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C> GO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,DLOP,MBE,RF,RF,RF,C,C> DGO0;
      typedef Dune::PDELab::GridOperator<GFS,GFS,MLOP,MBE,RF,RF,RF,C,C> GO1;
      typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
      typedef Dune::PDELab::OneStepGridOperator<DGO0,GO1,true> DIGO;
      typedef typename Dune::PDELab::BackendVectorSelector<GFS,RF>::Type V;

      typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
      typedef Dune::PDELab::SimpleTimeController<RF> TC;
      typedef Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,CV,CV,TC> OSM;

      // choose linear solver
      typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS,C> DLS;
      typedef Dune::PDELab::StationaryLinearProblemSolver<DIGO,DLS,V> CPDESOLVER;
      typedef Dune::PDELab::Alexander2Parameter<RF> DTimeSteppingMethod;

      typedef Dune::PDELab::OneStepMethod<RF,DIGO,CPDESOLVER,V,V> DOSM;

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
      typename std::enable_if<std::is_same<U, SOFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
      {
        return shared_ptr<FR>(new FR(cdgf,1.));
      }

      template< class U>
      typename std::enable_if<std::is_same<U, DFR>::value, shared_ptr<U> >::type getFlux(ConcDGF& cdgf)
      {

        return shared_ptr<FR>(new FR());
      }

      CombineTransportSolver(const GV &gv_, CTP& ctp_, ODE& ode_, const Dune::ParameterTree & param_, const std::string cname)
        : param(param_), gv(gv_), ctp(ctp_), ode(ode_),
          verbosity(param.sub("Input").get<int>("verbosity")),
          fem(getCube()),
          con(),
          gfs(gv,fem,con),
          cnew(gfs,0.),
          cold(gfs,0.),
          cadhesion(gfs,0.),
          conc_initial(gv,param,cname),
          concdgf(gfs,cold),
          concadhesion(gfs,cadhesion),
          flux(getFlux<FR>(concdgf)),
          lop(ctp,*flux),
          dlop(ctp),
          cflcontroller(gv,ctp,*flux),
          mlop(ctp),
          cg(getContainer()),
          mbe(5),
          go0(gfs,cg,gfs,cg,lop,mbe),
          dgo0(gfs,cg,gfs,cg,dlop,mbe),
          go1(gfs,cg,gfs,cg,mlop,mbe),
          igo(go0,go1),
          digo(dgo0,go1),
          tc(),
          ls(gfs),
          dls(gfs,cg,5000,0), //superLU
          timesteppingmethod(),
          dtimesteppingmethod(),
          osm(timesteppingmethod,igo,ls,tc),
          cpdesolver(digo,dls,1e-8,1e-14,0),
          dosm(dtimesteppingmethod,digo,cpdesolver),
          os_scheme(param.sub("Input").get<std::string>("osscheme")),
          tmax(0.0)
      {
        if (gv.comm().rank()>0)
          verbosity = 0;

        igo.setMethod(timesteppingmethod);
        osm.setVerbosityLevel(0);
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::interpolate(conc_initial,gfs,cadhesion);
        Dune::PDELab::constraints(conc_initial,gfs,cg,false);

        ctp.setTimeTarget(1.e100,1.e100);
        dosm.setVerbosityLevel(0);

      }

      template<typename PVD>
      void setOutput(PVD & pvdwriter)
      {
        pvdwriter.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>(concdgf,"concentration boilerplate"));
      }

      ConcDGF & getConcDGF(){return concdgf;}
      const ConcDGF & getConcDGF() const {return concdgf;}

      ConcDGF & getAdhesion(){return concadhesion;}
      const ConcDGF & getAdhesion() const {return concadhesion;}

      V & getSolution() {return cnew;}
      const V & getSolution() const {return cnew;}

      void setSolution(V & c) {cnew = c; cold = cnew;}

      RF getOutflow()
      {
        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename Dune::PDELab::ElementGeometry<Element> EG;
        ElementIterator tit = gv.template begin<0>();
        // loop once over the grid and find last element (in 1D)
        for (ElementIterator it = gv.template begin<0>();
             it!=gv.template end<0>(); ++it)
          {
            tit = it;
          }

        EG eg(*tit);
        // cell geometry
        const Dune::FieldVector<RF,dim>&
          cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
        typename ConcDGF::Traits::RangeType value;
        concdgf.evaluate(eg.entity(),cell_center_local,value);

        value/=ctp.getBgConcentration();

        return value;
      }

      //! set parameters for the model
      void setParamValue(std::string name, RF value)
      {
        if (name=="Dt")
          ctp.setDiffusion(value);
        else if (name=="velocity")
          {
            ctp.setVelocity(value);
            ode.getODEModel().getParameters().setVelocity(value);
          }
        else if (name=="katt")
          ode.getODEModel().getParameters().setKatt(value);
        else if (name=="Smax")
          ode.getODEModel().getParameters().setSmax(value);
        else if (name=="beta")
          ode.getODEModel().getParameters().setBeta(value);
        else if (name=="riso")
          {
            ode.getODEModel().getParameters().setRiso(value);
            ctp.setRiso(value);
          }
        else if (name=="kdet")
          ode.getODEModel().getParameters().setKdet(value);
        else
          std::cerr << "bad parameter set in parameter " << name << std::endl;
      }

      void reset()
      {
        Dune::PDELab::interpolate(conc_initial,gfs,cold);
        Dune::PDELab::interpolate(conc_initial,gfs,cnew);
        Dune::PDELab::interpolate(conc_initial,gfs,cadhesion);
      }


      void printSolution(std::ofstream& solution, RF timesteps = 100.)
      {
        if (verbosity)
          std::cout << "======= Adhesion Model, Solution ======= at time " << tmax << "\n";

        reset();
        RF reactiontime = 0;
        RF timestep = tmax/timesteps;
        timestep = std::max(timestep,cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl"));


        solution << "# time |  concentration" << "\n";
        solution << std::left << std::setw(10) << reactiontime << "  " << getOutflow()  << "\n";

        while (reactiontime<tmax-1e-8)
          {
            timestep = std::min(timestep,tmax-reactiontime);
            apply(reactiontime,timestep);
            reactiontime += timestep;
            solution << std::left << std::setw(10) << reactiontime << "  " << getOutflow()<< "\n";
          }
        solution << "\n\n";

        typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename Dune::PDELab::ElementGeometry<Element> EG;
        // loop once over the grid and find elements with given x_coord
        Dune::PDELab::ElementMapper<GV> cell_mapper(gv);

        for (ElementIterator it = gv.template begin<0>();
             it!=gv.template end<0>(); ++it)
          {
            EG eg(*it);
            // compute unique id
            typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
            solution << std::left << std::setw(10) << eg.geometry().center()[0] << " " << cadhesion.block(n) << "\n";
          }
        solution << "\n";

      }

      bool apply_A_R_D(RF time, RF timetarget)
      {
        ctp.onlyAdvection();
        RF timestep = std::min(cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl"),timetarget-time);

        try {
          while (time<timetarget-1.e-8)
            {
              ctp.reset();
              ctp.onlyAdvection();
              if (verbosity)
                std::cout << "======= solve transport part from " << time << " to " << time+timestep << " =======\n";
              watch.reset();
              ctp.onlyAdvection();
              osm.apply(time,timestep,cold,cnew);

              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;

              if (verbosity)
                std::cout << "reaction from " << time << " to " << time+timestep << " =======\n";
              watch.reset();
              std::vector<V*> datavector;
              datavector.push_back(&cnew);
              datavector.push_back(&cadhesion);
              ode.apply(time,timestep,datavector);
              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;

              cold=cnew;
              if (verbosity)
                std::cout << "======= solve diffusion part from " << time << " to " << time+timestep << " =======\n";
              watch.reset();
              ctp.onlyDiffusion();
              dosm.apply(time,timestep,cold,cnew);
              cold=cnew;
              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;

              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
            }
          ctp.reset();
          if (verbosity)
            std::cout <<"A_R_D finished " << std::endl;
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
          std::cerr << "Unknown exception thrown!"<< std::endl;
          return false;
        }
        return false;
      }

      bool apply_AD_R(RF time, RF timetarget)
      {
        RF timestep = std::min(cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl"),timetarget-time);

        try {
          while (time<timetarget-1.e-8)
            {
              ctp.reset();
              if (verbosity)
                std::cout << "======= solve transport and diffusion part from " << time << " to " << time+timestep << " =======\n";
              watch.reset();
              osm.apply(time,timestep,cold,cnew);

              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;

              if (verbosity)
                std::cout << "reaction from " << time << " to " << time+timestep << " =======\n";
              watch.reset();
              std::vector<V*> datavector;
              datavector.push_back(&cnew);
              datavector.push_back(&cadhesion);
              ode.apply(time,timestep,datavector);
              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;
              cold=cnew;

              if (verbosity)
                std::cout << "... done : " << watch.elapsed() << " s with timestep " << timestep << std::endl;

              time+=timestep;
              timestep = std::min(timestep, timetarget-time);
            }
          ctp.reset();
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
          std::cerr << "Unknown exception thrown!"<< std::endl;
          return false;
        }
        return false;
      }




      bool apply(RF time, RF dt)
      {
        cold = cnew;
        RF timetarget = time+dt;
        if (timetarget>tmax)
          tmax=timetarget;

        ctp.setTimeTarget(timetarget,dt);

        if (os_scheme=="AD_R")
          return apply_AD_R(time,timetarget);
        else if (os_scheme=="A_R_D")
          return apply_A_R_D(time,timetarget);
        else {
          std::cerr << "OS Scheme "<< os_scheme << " is not known";
          return false;
        }
      }



      const RF suggestTimeStep()
      {
        return cflcontroller.suggestTimeStep();
      }

      CTP & getCTP()
      {
        return ctp;
      }

      ODE & getODE()
      {
        return ode;
      }

      const Dune::ParameterTree & param;

    private:
      const GV& gv;
      CTP &ctp;
      ODE &ode;
      int verbosity;
      Dune::Timer watch;
      FEM fem;
      CON con;
      GFS gfs;
      V cnew;
      V cold;
      V cadhesion;
      InitialConcentration conc_initial;
      ConcDGF concdgf;
      ConcDGF concadhesion;
      shared_ptr<FR> flux;
      LOP lop;
      DLOP dlop;
      CFLController cflcontroller;
      MLOP mlop;
      C cg;
      MBE mbe;
      GO0 go0;
      DGO0 dgo0;
      GO1 go1;
      IGO igo;
      DIGO digo;
      TC tc;
      LS ls;
      DLS dls;
      TimeSteppingMethod timesteppingmethod;
      DTimeSteppingMethod dtimesteppingmethod;
      OSM osm;
      CPDESOLVER cpdesolver;
      DOSM dosm;
      std::string os_scheme;
      RF tmax;
    };


  } // End namespace Dycap
} // end namespace Dune

#endif
