// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/paamg/amg.hh>

#include<dune/dycap/newton/newton_utilities.hh>
#include<dune/dycap/newton/newton.hh>
#include<dune/dycap/oldpdelab/pdelab.hh>

#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/constraints/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/typetree/typetree.hh>


#include<dune/pdelab/backend/istlsolverbackend.hh>
#include<dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include<dune/dycap/output/pvdwriter.hh>
#include<dune/dycap/utilities/utilities.hh>

#include<dune/dycap/physics/hydraulicparameters.hh>
#include<dune/dycap/physics/physical_chemistry.hh>

#define PLPC TRUE

#if PLPC
#include<dune/dycap/models/pcpl_twophaseccfv_incompressible.hh>
#include<dune/dycap/parameters/plpc.hh>

#elif PLPG
#include<dune/dycap/models/twophaseccfv.hh>
#include<dune/dycap/parameters/plpg.hh>

#elif PLSG
#include<dune/dycap/models/plsg_twophaseccfv.hh>
#include<dune/dycap/parameters/plsg.hh>
#endif

#include<dune/dycap/parameters/twophaseparameters.hh>
#include<dune/dycap/utilities/multispantwophase.hh>

template<class G>
void run (const G& g, const HeleshawDomain<G::dimension> & domain, Dune::ParameterTree & param)
{
  
  // std::cout << "number of boundary segments " << g.numBoundarySegments() << std::endl;
  typedef typename G::LeafGridView GV;
  GV gv = g.leafGridView();

  // choose some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // verbosity
  const int verbosity = param.sub("Verbosity").get<int>("verbosity", 2);
  const int verbosityInstationary = param.sub("Verbosity").get<int>("Instationary", verbosity);

  // timeloop
  RF timestep = param.sub("Timeloop").get<RF>("timestep");

  // update timesteps according to refinement
  timestep /= 1<<domain.refine;


  /********************************************************************/
  /*********************** twophase problem ***************************/
  /********************************************************************/

  // Make grid function space
  Dune::GeometryType gt;
  gt.makeCube(dim);
  // finite volumes finite element map
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(gt);
  // only parallel constraints
  typedef Dune::PDELab::P0ParallelConstraints CON;
  // vector backend, blocksize 2
  typedef Dune::PDELab::istl::VectorBackend<> VBE0;
  // grid function space
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;
  typedef Dune::PDELab::ISTLVectorBackend
    <Dune::PDELab::ISTLParameters::static_blocking,2> VBE;
  // pgfs with blockwise mapper
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS,2,VBE,Dune::PDELab::EntityBlockedOrderingTag> TPGFS;
  if (gv.comm().rank() == 0)
    std::cout << "=== two phase function space setup\n";
  watch.reset();
  CON con;
  GFS gfs(gv,fem,con);
  TPGFS tpgfs(gfs);

  if (gv.comm().rank() == 0)
    std::cout << "=== multispan twophase parameter object\n";
  watch.reset();
  // make parameter object
  typedef MultiSpanTwoPhase<GV,RF> TimeManager;
  TimeManager timemanager(gv, param);


  // make parameter object
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp(gv, domain, param);
  if (gv.comm().rank() == 0)
    std::cout << "... done : " << watch.elapsed() << " s" << std::endl;

  typedef PermeabilityField<TP> PermField;
  PermField permfield(tp);

#if (PLPC || PLPCFLUX)
  // make subspaces for visualization
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0> > P_lSUB;
  P_lSUB p_lsub(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1> > P_cSUB;
  P_cSUB p_csub(tpgfs);


  // initial value function
  typedef P_l<GV,RF,TP> P_lType;
  P_lType p_l_initial(gv,tp);
  typedef P_c<GV,RF,TP> P_cType;
  P_cType p_c_initial(gv,tp);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,P_cType> PType;
  PType p_initial(p_l_initial,p_c_initial);



#elif PLPG
  // make subspaces for visualization
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0> > P_lSUB;
  P_lSUB p_lsub(tpgfs);
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1> > P_gSUB;
  P_gSUB p_gsub(tpgfs);

  // initial value function
  typedef P_l<GV,RF> P_lType;
  P_lType p_l_initial(gv,tp);
  typedef P_g<GV,RF> P_gType;
  P_gType p_g_initial(gv,tp);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,P_gType> PType;
  PType p_initial(p_l_initial,p_g_initial);

#elif PLSG
  // make subspaces for visualization
  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0> > P_lSUB;
  P_lSUB p_lsub(tpgfs);

  typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1> > S_gSUB;
  S_gSUB s_gsub(tpgfs);


  // initial value function
  typedef P_l<GV,RF> P_lType;
  P_lType p_l_initial(gv,tp);
  typedef S_g<GV,RF> S_gType;
  S_gType s_g_initial(gv,tp);
  typedef Dune::PDELab::CompositeGridFunction<P_lType,S_gType> PType;
  PType p_initial(p_l_initial,s_g_initial);
#endif



  // make constraints map and initialize it from a function
  typedef typename TPGFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  cg.clear();
  Dune::PDELab::constraints(p_initial,tpgfs,cg,false);

#if PLPCFLUX
  // make grid operator space
  typedef Dune::PDELab::TwoPhaseTotalFlux<TP> LOP;
  typedef Dune::PDELab::TwoPhaseTotalFluxTemporal<TP> MLOP;
#else
  typedef Dune::PDELab::TwoPhaseTwoPointFluxOperator<TP> LOP;
  typedef Dune::PDELab::TwoPhaseOnePointTemporalOperator<TP> MLOP;
#endif

  LOP lop(tp);
  MLOP mlop(tp);

  // matrix backend
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5);
  Dune::PDELab::ImplicitEulerParameter<RF> method;
  typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MBE,RF,RF,RF,C,C> GO0;
  GO0 go0(tpgfs,cg,tpgfs,cg,lop,mbe);
  typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,MLOP,MBE,RF,RF,RF,C,C> GO1;

  GO1 go1(tpgfs,cg,tpgfs,cg,mlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics (do not call this for IGO before osm.apply() was called!)
  typename GO0::Traits::Jacobian jac(go0);
  std::cout << jac.patternStatistics() << std::endl;

  // make vector for old time step and initialize
  typedef typename IGO::Traits::Domain V;
  V pold(tpgfs, 0.0);
  Dune::PDELab::interpolate(p_initial,tpgfs,pold);

  // make vector for new time step and initialize
  V pnew(tpgfs);
  pnew = pold;

  if (gv.comm().rank() == 0)
    std::cout << "... done : " << watch.elapsed() << " s" << std::endl;


  /*******************************************************************/
  /**************************   SOLVERS   ***************************/
  /*******************************************************************/
  if (gv.comm().rank() == 0)
    std::cout << "=== solvers setup\n";
  watch.reset();

  //typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  LS ls(tpgfs,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0));

  typedef Dune::PDELab::Newton<IGO,LS,V> PDESOLVER;
  PDESOLVER newton(igo,ls);
  // Newton parameters class
  typedef Dune::PDELab::NewtonParameters NewtonParameters;
  NewtonParameters newtonparameters(param.sub("Newton"));
  // set new parameters
  newtonparameters.set(newton);
  
  // time-stepper
  Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,V,V> osm(method,igo,newton);
  osm.setVerbosityLevel(verbosityInstationary);



  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;



  /*******************************************************************/
  /************************** COMPUTATION  ***************************/
  /*******************************************************************/
  
  int counter = 0;
  // get initial timestep
  timemanager.set_dt(timestep);
  RF tstep = timemanager.getTimeStepSize();
  
  std::cout << "timestep init " << tstep << std::endl;
  while (!timemanager.finalize())
    {
      tstep = timemanager.getTimeStepSize();
      
      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep " << counter  << "\n";

      // initial solution of two phase problem
      if (timemanager.compute("computeInitial")) {
        if (gv.comm().rank() == 0)
          std::cout << "=== compute Initial Solution\n";
        watch.reset();

        //        typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<GO0> LS;
        typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GO0> LS;
        LS ls(tpgfs,500,param.sub("Newton").get<int>("LinearVerbosity", 0));

        typedef Dune::PDELab::Newton<GO0,LS,V> PDESOLVER;
        PDESOLVER newton(go0,ls);

        typedef Dune::PDELab::Newton<GO0,LS,V> TPDESOLVER;
        TPDESOLVER tnewton(go0,ls );
        newtonparameters.set(newton);
        newton.setLineSearchStrategy("hackbuschReuskenAcceptBest");
        newton.setLineSearchMaxIterations(60);
        newton.setAbsoluteLimit(1.e-14);

        if (gv.comm().rank() == 0)
          std::cout << "... apply Newton\n";
        newton.apply(pnew);
        pold = pnew;
        if (gv.comm().rank() == 0)
          std::cout << "... done : " << watch.elapsed() << " s" << std::endl;
      }


      // solve two phase problem
      if (timemanager.compute("computeTwoPhase")) {
        if (gv.comm().rank() == 0)
          std::cout << "======= solve two phase problem =======\n";
        watch.reset();
        try {
          osm.apply(timemanager.getTime(),tstep,pold,pnew);
          //  pold = pnew;
          if (gv.comm().rank() == 0)
            std::cout << "... done\n";
          //  velocity_controller.control(rt0_l);
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
          if (gv.comm().rank() == 0)
            std::cout << "Newton Linesearch Error" << std::endl;
          timemanager.notifyFailure();
          pnew = pold;
          continue;
        }
        catch (Dune::PDELab::NewtonNotConverged) {
          if (gv.comm().rank() == 0)
            std::cout << "Newton Convergence Error" << std::endl;
          timemanager.notifyFailure();
          pnew = pold;
          continue;
        }
        catch (Dune::PDELab::NewtonLinearSolverError) {
          if (gv.comm().rank() == 0)
            std::cout << "Newton Linear Solver Error" << std::endl;
          timemanager.notifyFailure();
          pnew = pold;
          continue;
        }
        catch (Dune::ISTLError) {
          if (gv.comm().rank() == 0)
            std::cout << "ISTL Error" << std::endl;
          timemanager.notifyFailure();
          pnew = pold;
          continue;
        }

        if (gv.comm().rank() == 0)
          std::cout << "... took : " << watch.elapsed() << " s" << std::endl;
      }

      // notify success in this timestep
      timemanager.notifySuccess(5);
      pold = pnew;

      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep done, T= " << timemanager.getTime() << " =============\n\n";
      counter++;
    }

}

int rank;
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      {
        if(helper.rank()==0)
          std::cout << "parallel run on " << helper.size() << " processes" << std::endl;
      }
    rank = helper.rank();

    std::string configfile = argv[0]; configfile += ".conf";

    // parse cmd line parameters
    if (argc > 1 && argv[1][0] != '-')
      configfile = argv[1];
    for (int i = 1; i < argc; i++)
      {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
          {
            if(helper.rank()==0)
              std::cout << "usage: ./twophase_test <configfile> [OPTIONS]" << std::endl;
            return 0;
          }
      }
    double start=MPI_Wtime();

    // read parameters from file
    Dune::ParameterTree param;
    Dune::ParameterTreeParser parser;
    parser.readINITree(configfile, param);
    addConfigFiles(param);

    int dim=param.sub("Domain").get<int>("dim", 2);

    // 2D
    if (dim==2)
      {
        HeleshawDomain<2> domain(param.sub("Domain"));

        // make grid
        Dune::FieldVector<double,2> L;
        L[0] = domain.width;
        L[1] = domain.height;
        Dune::array<int,2> N;
        N[0] = domain.nx;
        N[1] = domain.ny;
        std::bitset<2> periodic(false);
        int overlap = param.get<int>("overlap", 0);;
        Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
        grid.globalRefine(domain.refine);
        // solve problem
        run(grid,domain,param);
      }
    if (dim==3)
      {
        HeleshawDomain<3> domain(param.sub("Domain"));

        // make grid
        Dune::FieldVector<double,3> L;
        L[0] = domain.width; 
        L[1] = domain.depth; 
        L[2] = domain.height;
        Dune::array<int,3> N;
        N[0] = domain.nx;
        N[1] = domain.nz;
        N[2] = domain.ny;
        std::bitset<3> periodic(false);
        int overlap = param.get<int>("overlap", 0);;
        Dune::YaspGrid<3> grid(L,N,periodic,overlap);
        grid.globalRefine(domain.refine);
        // solve problem
        run(grid,domain,param);
      }

    if(helper.rank()==0)
      std::cout<<"Total computation time was "<<MPI_Wtime()-start
               <<" seconds."<<std::endl;

    // test passed
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}


