#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#ifdef _PARALLEL_
#include<mpi.h>
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <dune/common/timer.hh>
#include<dune/common/exceptions.hh>
#include<src/output/output.hh>

#define TWOCOMPONENTS FALSE

#include"model.hh"
#include <src/estimation/utilities.hh>
#include"fitdataclass.hh"
#include"fitmodel.hh"

#ifdef _GAUSSNEWTON_
#include <src/estimation/gaussnewton.hh>
#else
#include <src/estimation/levemberg_marquardt.hh>
#endif

using namespace std;

// main program for parameter estimation
int main (int argc, char **argv)
{

  //Maybe initialize Mpi
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
  if(Dune::MPIHelper::isFake)
    std::cout<< "This is a sequential program." << std::endl;
  else
    {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " processes" << std::endl;
    }

  std::string configfile = argv[0]; configfile += ".conf";

  // parse cmd line parameters
  if (argc > 1 && argv[1][0] != '-')
    configfile = argv[1];
  for (int i = 1; i < argc; i++)
    {
      if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
          if(helper.rank()==0)
            std::cout << "usage: ./parameter_estimation_anaerobic <configfile> [OPTIONS]" << std::endl;
          return 0;
        }
    }

  // read parameters from file
  Dune::ParameterTree param;
  Dune::ParameterTreeParser parser;
  parser.readINITree(configfile, param);

  // 1D domain
  HeleshawDomain<1> domain(param.sub("Domain"));

  // make grid
  Dune::FieldVector<double,1> L;
  L[0] = domain.width;

  Dune::array<int,1> N;
  N[0] = domain.nx;

  std::bitset<1> periodic(false);
  int overlap = param.get<int>("overlap", 0);
  typedef typename Dune::YaspGrid<1> Grid;
  Grid grid(L,N,periodic,overlap);
  grid.globalRefine(domain.refine);

  typedef typename  Grid::LeafGridView GV;
  // choose some types
  //typedef typename GV::Grid::ctype DF;
  typedef double RF;

  GV gv=grid.leafGridView();

  if(argc < 2)
    {
      cout << "Usage: parameter_estimation_anaerobic <datafilename> " << endl;
      cout << "   (e.g. parameter_estimation_aerobic data.dat )" << endl;
      throw error_abort();
    }

  std::cout << "read main input file" << std::endl;
  string dataFileName;
  if (argc>2)
    dataFileName=argv[2];
  else
    dataFileName = param.sub("Input").get<std::string>("datafile");
  std::cout << "data " << dataFileName << std::endl;


#if TWOCOMPONENTS
  //sort_experiment_data(param,dataFileName);
#endif

  string parameterFileName;
  if (argc>3)
    parameterFileName=argv[3];
  else
    parameterFileName = param.sub("Input").get<std::string>("parameterfile");
  std::cout << "parameter file " << dataFileName << std::endl;

  int maxIterations=param.sub("Init").get<int>("maxiterations");
  std::cout << "maxit" << maxIterations << std::endl;
  double initialLambda=param.sub("Init").get<double>("lambda");

  double precision=param.sub("Init").get<double>("precision");
  std::cout << "maxit " << maxIterations << " lambda " << initialLambda << " precision " << precision << std::endl;



  int verbosity =  param.sub("Input").get<int>("verbosity");

  std::cout << "construct fitData class" << std::endl;
  // read data
  DycapFitDataClass fitData(dataFileName,1);
  std::cout << "construct fitData class finished" << std::endl;

  /*******************************************************************/
  /************************** reactive problem****************************/
  /*******************************************************************/

  std::cout << "construct reactionparameter " << std::endl;
  // make reaction parameter object
  typedef ReactionParameter<GV,RF> RP;
  RP rp(gv, param);

  // setup ODE problem
  typedef AdhesionModel<RP> ODEProblem;
  ODEProblem rmodel(rp,1);
  rmodel.setVerbosityLevel(verbosity);

  std::cout << "construct odeproblem " << std::endl;
  typedef OdeBase<ODEProblem> ODESolver;
  // setup ODE solver
  std::shared_ptr<ODESolver> odesolverp = std::shared_ptr<ODESolver>(new DummyOdeBase<ODEProblem>(rmodel));

  typedef EcoliGrowthGridOperator<RF,GV, ODEProblem, ODESolver> ODEGOS;
  ODEGOS odegos(gv,rmodel,*odesolverp,verbosity);

  OdeTimeSteppingMethods<ODEProblem> odemethods(rmodel);
  odemethods.setTimestepMethod(odegos, param.sub("Setup").get<std::string>("odesolver"));

  /********************************************************************/
  /*********************** transport problem  *************************/
  /********************************************************************/


  std::cout << "construct ctp " << std::endl;
  //  make grid operator space for water transport problem
  typedef TransportLiquid<GV,RF> CTP;
  CTP ctp1(gv, param,"TransportC1","intervalsc1");
  std::cout << "construct solver " << std::endl;

#if TWOCOMPONENTS
  CTP ctp2(gv, param,"TransportC2", "intervalsc2");

  ODEGOS odegos2(gv,rmodel,*odesolverp,verbosity);
  odemethods.setTimestepMethod(odegos2, param.sub("Setup").get<std::string>("odesolver"));
  typedef typename Dune::Dycap::CombineTransportSolver<GV,CTP,ODEGOS> Model;
  Model model1(gv,ctp1,odegos,param,"TransportC1");
  Model model2(gv,ctp2,odegos2,param,"TransportC2");

  typedef typename Dune::Dycap::TwoComponentSolver<Model> TCModel;
  TCModel tcmodel(model1,model2);
#else

  // create object which is able to solve advection-reaction problem
  // set new parameters, compare with data and compute again
  typedef typename Dune::Dycap::CombineTransportSolver<GV,CTP,ODEGOS> TCModel;
  TCModel tcmodel(gv,ctp1,odegos,param,"TransportC1");
#endif


  // parameter estimation procedure
  try {
    //
    std::cout << "construct fitModel class" << std::endl;
    // Init fit routine
    FitModel<TCModel> fitModel(tcmodel,parameterFileName,dataFileName);



#ifdef _GAUSSNEWTON_
    std::cout << "Construct GaussNewtonClass" << std::endl;
    GaussNewtonClass fitObject(fitData,fitModel);
#else
    std::cout << "Construct LevMarqClass" << std::endl;
    LevMarqClass fitObject(fitData,fitModel);
#endif
    fitObject.SetMaxIterations(maxIterations);
    fitObject.SetAbsLimit(1e-15);
#ifndef _GAUSSNEWTON_
    // set parameter for levenberg marquardt
    fitObject.SetInitialLambda(initialLambda);
#endif

    time_t rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    std::cout << "Starting optimization at " << asctime(timeinfo) << std::endl;
    Dune::Timer timer;
    timer.reset();

    // fit with a given precision
    FLOAT bestResiduum=fitObject.fit(precision);
    double elapsed = timer.elapsed();
    int hours = floor(elapsed/3600);
    int mins = floor((elapsed-hours*3600)/60);
    int secs = round(elapsed - 3600*hours - 60*mins);


    // output results
    cout << fitObject;
    cout << "Number of iterations used: " << fitObject.getNumberOfIterations()  << " with residuum " << bestResiduum << endl;
    cout << "Time elapsed: " << hours << "h " << mins << "min "
         << secs << "sec"<< std::endl;


    parameterFileName+=".new";
    fitModel.WriteParameterFile(parameterFileName);


    FLOAT timesteps =  param.sub("Input").get<FLOAT>("timesteps");
    fitModel.SolutionOutput(fitData,timesteps);
    //visualise(fitObject,dataFileName);

  }
  catch (Dune::Exception e) {
    std::cout << e << std::endl;
  }



#if TWOCOMPONENTS
  GnuplotSolution<GV> gnuplot(gv);
  std::ostringstream s1;
  s1 << "simulation1_result_end.dat";
  gnuplot.write(s1.str(), tcmodel.getConcDGF(1), tcmodel.getAdhesion(1));
  std::ostringstream s2;
  s2 << "simulation2_result_end.dat";
  gnuplot.write(s2.str(), tcmodel.getConcDGF(2), tcmodel.getAdhesion(2));
#else
  GnuplotSolution<GV> gnuplot(gv);
  std::ostringstream s1;
  s1 << "simulation1_result_end.dat";
  gnuplot.write(s1.str(), tcmodel.getConcDGF(), tcmodel.getAdhesion());
#endif

}
