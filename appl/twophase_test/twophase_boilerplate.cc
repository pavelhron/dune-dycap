// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>

#include<src/boilerplate/twophasesolver_pcpl.hh>

//==============================================================================
// driver
//==============================================================================

template<class G, class Domain>
void run (const G& g, Domain & domain, Dune::ParameterTree & param)
{

  typedef typename G::LeafGridView GV;
  GV gv = g.leafGridView();

  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;

  // verbosity
  //const int verbosity = param.sub("Verbosity").get<int>("verbosity", 2);
  const bool graphics = param.get<bool>("writeVTK", false);

  // timeloop
  RF timestep = param.sub("Timeloop").get<RF>("timestep");

  // update timesteps according to refinemen
  timestep /= 1<<domain.refine;

  /********************************************************************/
  /*********************** twophase problem ***************************/
  /********************************************************************/

  if (gv.comm().rank() == 0)
    std::cout << "=== multispan time manager object\n";
  watch.reset();
  // make parameter object
  typedef MultiSpanTwoPhase<GV,RF> TimeManager;
  TimeManager timemanager(gv, param);

  if (gv.comm().rank() == 0)
    std::cout << "... done : " << watch.elapsed() << " s" << std::endl;


  if (gv.comm().rank() == 0)
    std::cout << "=== multispan twophase parameter object\n";
  watch.reset();
  // make parameter object
  typedef TwoPhaseParameter<GV,RF> TP;
  TP tp(gv, domain, param);
  if (gv.comm().rank() == 0)
    std::cout << "... done : " << watch.elapsed() << " s" << std::endl;

  if (gv.comm().rank() == 0)
    std::cout << "=== twophase solver object\n";
  watch.reset();
  typedef Dune::Dycap::TwoPhaseSolver<GV,TP> TPS;
  TPS tps(gv,tp,param);
  if (gv.comm().rank() == 0)
    std::cout << "... done : " << watch.elapsed() << " s" << std::endl;



  /*******************************************************************/
  /************************** VTK OUTPUT  ****************************/
  /*******************************************************************/

  if (gv.comm().rank() == 0)
    std::cout << "=== absolute discrete grid function setup\n";
  watch.reset();


  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  if (gv.comm().rank() == 0)
    std::cout << "=== output setup\n";
  watch.reset();

  char basename[255];
  std::string material = param.get<std::string>("Setup.material");
  sprintf(basename,"%s-%s-%s-%01dl%01dd",param.get<std::string>("VTKname","").c_str(),"boilerplate",param.sub(material).get<std::string>("model","").c_str(),domain.refine,dim);



  typedef typename  Dune::SubsamplingVTKWriter<GV> SVTKWriter;
  typedef Dune::PVDWriter<GV, SVTKWriter > PVDWriter;
  PVDWriter pvdwriter(gv,0,basename,Dune::VTK::conforming);
  tps.setOutput(pvdwriter);

  // object to store timesteps
  TimeStepWriter<RF> timestepwriter(basename, gv.comm().rank());

  RF contentfactor= param.sub(material).get<RF>("porosity")*1.e6; // in ml!!
  if (dim<3)
    contentfactor*=domain.depth;
  RF watercontent = computeMass(tps.getS_ldgf(),contentfactor);

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;



  /*******************************************************************/
  /************************** COMPUTATION  ***************************/
  /*******************************************************************/

  if (graphics)
    pvdwriter.write(0);

  // VTK output every mod time steps
  int counter = 1;

  // get initial timestep
  timemanager.set_dt(timestep);
  RF tstep = timemanager.getTimeStepSize();

  RF newwc = computeMass(tps.getS_ldgf(),contentfactor);
  if (gv.comm().rank() == 0)
    std::cout << "watercontent at time " << timemanager.getTime() << " is " << newwc<< " ml " <<  newwc-watercontent<<" ml\n";

  while (!timemanager.finalize())
    {
      tstep = timemanager.getTimeStepSize();
      timestepwriter.add_timestep(timemanager.getTime(),tstep,timemanager.getDtMin(), timemanager.getDtMax());

      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep " << counter  << "\n";

      // initial solution of two phase problem
      if (timemanager.compute("computeInitial")) {
        tps.applyStationary(timemanager.getTime());
      }


      // solve two phase problem
      else if (timemanager.compute("computeTwoPhase")) {
        if(tps.apply(timemanager.getTime(),tstep))
          {
            if (gv.comm().rank() == 0)
              {
                std::cout << "twophase problem solved" << std::endl;
              }
          }
        else
          {
            timemanager.notifyFailure();
            if (gv.comm().rank() == 0)
              std::cout << "twophase problem NOT solved " << std::endl;
            continue;
          }
      }
      else
        {
          tps.applyZeroVelocity(timemanager.getTime());
          if (gv.comm().rank() == 0)
            std::cout << "twophase with ZERO velocity" << std::endl;
        }


      // notify success in this timestep
      timemanager.notifySuccess(5);

      timestepwriter.add_timestep(timemanager.getTime(),tstep,timemanager.getDtMin(), timemanager.getDtMax());

      if (graphics && timemanager.isTimeForOutput())
        {
          pvdwriter.write(timemanager.getTime());
        }

      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep done, T= " << timemanager.getTime() << " =============\n\n";
      counter++;
    }

}

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

    std::string configfile = argv[0]; configfile += ".ini";

    // parse cmd line parameters
    if (argc > 1 && argv[1][0] != '-')
      configfile = argv[1];
    for (int i = 1; i < argc; i++)
      {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
          {
            if(helper.rank()==0)
              std::cout << "usage: ./twophase_boilerplate <configfile> [OPTIONS]" << std::endl;
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
        Dune::YaspGrid<2> grid(L,N,periodic,overlap,Dune::MPIHelper::getCollectiveCommunication());
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
