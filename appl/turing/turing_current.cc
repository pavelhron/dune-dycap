// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve two-component diffusion-reaction system
    with finite volume scheme. This class tests the
    multi domain grid.
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>


#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

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
#include<dune/pdelab/gridfunctionspace/subspace.hh>

#include "local_operator.hh"
#include "timemanager.hh"

#include<src/utilities/utilities.hh>
#include<src/utilities/controller.hh>
#include<src/newton/newton_reaction.hh>
#include<src/newton/newton_utilities.hh>
#include<src/output/output.hh>

#include "componentparameters.hh"
#include "initial_conditions.hh"
#include "reactionadapter.hh"


template<class GV>
void run (const GV& gv, Dune::ParameterTree & param)
{
    const bool verbosity = param.sub("Verbosity").get<bool>("verbosity", false);
    if (verbosity)
    Dune::gridinfo(gv.grid());
    typedef typename GV::Grid::ctype DF;
    typedef double RF;
    const int dim = GV::dimension;
    Dune::Timer watch;


    // output to vtk file (true/false)
    const bool graphics = param.get<bool>("writeVTK", false);

    ConvergenceAdaptiveTimeStepper<RF> timemanager(param);

    // for each component parameters different classes
    typedef DiffusionParameter<GV,RF> DP;
    DP dp1(param, "species_2");
    DP dp2(param, "species_1");


    typedef Dune::PDELab::MulticomponentDiffusion<GV,RF,DP,2> VCT;
    VCT vct(dp1,dp2);

    // <<<2>>> Make grid function space
    typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
    FEM fem(Dune::GeometryType(Dune::GeometryType::cube,dim));
    typedef Dune::PDELab::P0ParallelConstraints CON;
    typedef Dune::PDELab::istl::VectorBackend<> VBE0;
    typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS;

    typedef Dune::PDELab::istl::VectorBackend
      <Dune::PDELab::istl::Blocking::fixed,2> VBE;
    typedef Dune::PDELab::PowerGridFunctionSpace<GFS,VCT::COMPONENTS,VBE,
            Dune::PDELab::EntityBlockedOrderingTag> TPGFS;
    watch.reset();
    CON con;
    GFS gfs(gv,fem,con);
    TPGFS tpgfs(gfs);
    if (verbosity) std::cout << "=== function space setup " <<  watch.elapsed() << " s" << std::endl;

    // some informations
    gfs.update();
    if (verbosity) std::cout << "number of DOF =" << gfs.globalSize() << std::endl;


    // <<<2b>>> make subspaces for visualization
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<0> > U0SUB;
    U0SUB u0sub(tpgfs);
    typedef Dune::PDELab::GridFunctionSubSpace<TPGFS,Dune::TypeTree::TreePath<1> > U1SUB;
    U1SUB u1sub(tpgfs);


    // make constraints map and initialize it from a function (boundary conditions)
    typedef typename TPGFS::template ConstraintsContainer<RF>::Type CC;
    CC cg;
    cg.clear();
    Dune::PDELab::constraints(tpgfs,cg,false);

    typedef ReactionAdapter<RF> RA;
    RA ra(param);

    typedef Dune::PDELab::MulticomponentCCFVSpatialDiffusionOperator<VCT,RA> LOP;
    LOP lop(vct,ra); //local operator including reaction adapter
    typedef Dune::PDELab::MulticomponentCCFVTemporalOperator<VCT> TLOP;
    TLOP tlop(vct);
    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().
    // grid operators
    typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,LOP,MBE,RF,RF,RF,CC,CC> CGO0;
    CGO0 cgo0(tpgfs,cg,tpgfs,cg,lop,mbe);
    typedef Dune::PDELab::GridOperator<TPGFS,TPGFS,TLOP,MBE,RF,RF,RF,CC,CC> CGO1;
    CGO1 cgo1(tpgfs,cg,tpgfs,cg,tlop,mbe);
    // one step grid operator for instationary problems

    typedef Dune::PDELab::OneStepGridOperator<CGO0,CGO1,true> IGO;
    IGO igo(cgo0,cgo1);


    // <<<7>>> make vector for old time step and initialize
    typedef typename IGO::Traits::Domain V;
    V uold(tpgfs,0.0);
    V unew(tpgfs,0.0);

    // initial conditions
    typedef U0Initial<GV,RF> U0InitialType;
    U0InitialType u0initial(gv,param);
    typedef U1Initial<GV,RF> U1InitialType;
    U1InitialType u1initial(gv,param);


    typedef Dune::PDELab::CompositeGridFunction<U0InitialType,U1InitialType> UInitialType; //new
    UInitialType uinitial(u0initial,u1initial); //new
    Dune::PDELab::interpolate(uinitial,tpgfs,uold);
    unew = uold;

  // <<<5>>> Select a linear solver backend
#if HAVE_MPI
#if HAVE_SUPERLU
  typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS;
  LS ls(tpgfs,param.sub("Newton").get<int>("LSMaxIterations", 100),param.sub("Newton").get<int>("LinearVerbosity", 0));
#else
  typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
  LS ls(tpgfs,cc,5000,5,param.sub("Newton").get<int>("LinearVerbosity", 0));
#endif
#else //!parallel
#if HAVE_SUPERLU
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
    LS ls (param.sub("Newton").get<int>("LinearVerbosity", 0));
#else
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,param.sub("Newton").get<int>("LinearVerbosity", 0));
#endif
#endif

    // <<<6>>> Solver for non-linear problem per stage
    typedef Dune::PDELab::NewtonReaction<IGO,LS,V> PDESOLVER;
    PDESOLVER pdesolver(igo,ls);
    pdesolver.setLineSearchStrategy("hackbuschReuskenAcceptBest"); //strategy for linesearch
    typedef Dune::PDELab::NewtonParameters NewtonParameters;
    NewtonParameters newtonparameters(param.sub("Newton"));
    newtonparameters.set(pdesolver);

    // <<<7>>> time-stepper
    Dune::PDELab::Alexander2Parameter<RF> method;
    Dune::PDELab::OneStepMethod<RF,IGO,PDESOLVER,V,V> osm(method,igo,pdesolver);
    Dune::PDELab::TimeSteppingMethods<RF> tsmethods;
    tsmethods.setTimestepMethod(osm,param.get<std::string>("timesolver"));

    osm.setVerbosityLevel(param.sub("Verbosity").get<int>("Instationary", 0));

    char basename[255];
    sprintf(basename,"%s-%01d",param.get<std::string>("VTKname","").c_str(),param.sub("Domain").get<int>("refine"));

    typedef Dune::PVDWriter<GV> PVDWriter;
    PVDWriter pvdwriter(gv,basename,Dune::VTK::conforming);
    // discrete grid functions
    typedef Dune::PDELab::DiscreteGridFunction<U0SUB,V> U0DGF;
    U0DGF u0dgf(u0sub,uold);
    pvdwriter.addVertexData(std::make_shared<typename Dune::PDELab::VTKGridFunctionAdapter<U0DGF>>(u0dgf,"species_2"));
    typedef Dune::PDELab::DiscreteGridFunction<U1SUB,V> U1DGF;
    U1DGF u1dgf(u1sub,uold);
    pvdwriter.addVertexData(std::make_shared<typename Dune::PDELab::VTKGridFunctionAdapter<U1DGF>>(u1dgf,"species_1"));



    // <<<9>>> time loop

    if (graphics)
        pvdwriter.write(0);


    double dt = timemanager.getTimeStepSize();
    while (!timemanager.finalize())
    {
        dt = timemanager.getTimeStepSize();
        if (gv.comm().rank() == 0)
        {
            if (verbosity)
            {
                std::cout << "======= solve reaction problem ======= from " << timemanager.getTime() << "\n";
                std::cout << "time " << timemanager.getTime() << " dt " << timemanager.getTimeStepSize() << std::endl;
            }
        }

        watch.reset();
        try
        {
            // do time step
            osm.apply(timemanager.getTime(),dt,uold,unew);

            if (!controlReactionTimeStep(gv, unew))
            {
                timemanager.notifyFailure();
                unew = uold;
                continue;

            }
            uold = unew;

            if (gv.comm().rank() == 0 && verbosity)
                std::cout << "... done\n";
        }
        // newton linear search error
        catch (Dune::PDELab::NewtonLineSearchError) {
            if (gv.comm().rank() == 0 && verbosity)
                std::cout << "Newton Linesearch Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::PDELab::NewtonNotConverged) {
            if (gv.comm().rank() == 0 && verbosity)
                std::cout << "Newton Convergence Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::PDELab::NewtonLinearSolverError) {
            if (gv.comm().rank() == 0 && verbosity )
                std::cout << "Newton Linear Solver Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }
        catch (Dune::ISTLError) {
            if (gv.comm().rank() == 0 && verbosity )
                std::cout << "ISTL Error" << std::endl;
            timemanager.notifyFailure();
            unew = uold;
            if (!verbosity)
            {
              std::cout << "x" << std::flush;
            }
            continue;
        }

        if (gv.comm().rank() == 0 && verbosity)
            std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

        // notify success in this timestep
        timemanager.notifySuccess(5);

        if (graphics && timemanager.isTimeForOutput())
        {
          pvdwriter.write(timemanager.getTime());
          if (!verbosity)
            std::cout << "o" << std::flush;
        }
        else if (!verbosity)
        {
          std::cout << "." << std::flush;
        }
    }


    std::cout << std::endl;



}


bool isInside(const Dune::FieldVector<double, 2>& point, const Dune::FieldVector<double, 2>& dimension)
{

    const auto& x = point[0];
    const auto& y = point[1];
    const auto& width = dimension[0];
    const auto& height = dimension[1];

    bool inside=(((x >= 1.00000000000000000E+000 && x <= 4.90000000000000000E+001 && y >= 1.00000000000000000E+000 && y <= 4.90000000000000000E+001)) ? 1.00000000000000000E+000 : 0.00000000000000000E+000);
    return inside;


}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
    try{
        //Maybe initialize Mpi
        Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

        std::string configfile = argv[0]; configfile += ".conf";

        // parse cmd line parameters
        if (argc > 1 && argv[1][0] != '-')
            configfile = argv[1];
        for (int i = 1; i < argc; i++)
        {
            if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
            {
                if(helper.rank()==0)
                    std::cout << "usage: ./turing_current <configfile> [OPTIONS]" << std::endl;
                return 0;
            }
        }

        // read parameters from file
        Dune::ParameterTree param;
        Dune::ParameterTreeParser parser;
        parser.readINITree(configfile, param);

        int dim=param.sub("Domain").get<int>("dim", 2);

        if (dim==2)
        {
            // make grid
            Dune::FieldVector<double,2> L;
            L[0] = param.get<double>("Domain.height");
            L[1] = param.get<double>("Domain.width");

            std::array<int,2> N;
            N[0] = param.get<int>("Domain.nx");
            N[1] = param.get<int>("Domain.ny");

            std::bitset<2> periodic(false);
            int overlap = param.get<int>("overlap");
            Dune::YaspGrid<2> grid(L,N,periodic,overlap,Dune::MPIHelper::getCollectiveCommunication());
            grid.globalRefine(param.get<int>("Domain.refine",0));

            typedef Dune::YaspGrid<2>::LeafGridView GV;
            const GV& gv=grid.leafGridView();

            std::cout << "running problem " << std::endl;
            // solve problem
            run(gv,param);
        }
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    //catch (...){
    //  std::cerr << "Unknown exception thrown!" << std::endl;
    //  return 1;
    //}
}
