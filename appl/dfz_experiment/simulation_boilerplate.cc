// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include<vector>
#include<map>

#include<src/newton/newton.hh>
#include<src/output/pvdwriter.hh>
#include<src/output/basicio.hh>
#include<src/output/output.hh>

#include<src/utilities/utilities.hh>
#include<src/utilities/multispantwophase.hh>
#include<src/utilities/cfltransportcontroller.hh>

#include<src/parameters/transportparameters.hh>
#include<src/parameters/multitransportparameters.hh>
#include<src/parameters/reactionparameters.hh>

#include<src/ode/ode.hh>
#include<src/ode/odegridoperator.hh>

#include<src/boilerplate/twophasesolver_pcpl.hh>
#include<src/boilerplate/twophasevelocity.hh>
#include<src/boilerplate/transportsolver.hh>
#include<src/boilerplate/multicomponenttransportsolver.hh>


// used reaction model which was
#include"reaction_model.hh"

#define SECONDORDER true
#if SECONDORDER
#warning ***** WARNING ***** second order ******
#else
#warning ***** WARNING ***** first order  ******
#endif

template<typename  Sl, typename O2water, typename O2air, typename RF>
class O2
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename Sl::Traits::GridViewType,
                                                                           typename Sl::Traits::RangeFieldType,Sl::Traits::dimRange,
                                                                           typename Sl::Traits::RangeType>,
                                          O2<Sl,O2water,O2air,RF> >
{
  const Sl& sl;
  const O2water& o2water;
  const O2air& o2air;
  const RF K_H;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename Sl::Traits::GridViewType,
                                           typename Sl::Traits::RangeFieldType,Sl::Traits::dimRange,
                                           typename Sl::Traits::RangeType> Traits;


  O2 (const Sl& sl_,  const O2water& o2water_, const O2air& o2air_, const RF K_H_) : sl(sl_), o2water(o2water_), o2air(o2air_), K_H(K_H_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename Sl::Traits::RangeType sl_value, sg_value;
    typename O2air::Traits::RangeType o2air_value;
    typename O2water::Traits::RangeType o2water_value;
    sl.evaluate(e,x,sl_value);
    sg_value = 1.-sl_value;
    o2air.evaluate(e,x,o2air_value);
    o2water.evaluate(e,x,o2water_value);
    const RF Mo = 32.; // g/mol
    o2air_value*=Mo;
    Dune::FieldMatrix<RF,2,2> A;
    Dune::FieldVector<RF,2> b;
    Dune::FieldVector<RF,2> xx(0.);

    if (sg_value < 0.00001)
      sg_value = 0.00001;

        xx[0] = o2water_value;
        xx[1] = o2air_value;

        A[0][0] =sl_value;
        A[0][1] =sg_value;

        A[1][0] = 0.0;
        A[1][1] = 0.0;
        b = 0.0;
        A.umv(xx,b);
        A[1][0] = 1.0;
        A[1][1] = -K_H;
        A.solve(xx,b);
        y = xx[0];
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return sl.getGridView();
  }
};


void timeFormat(int reatime)
{
  int hour = reatime/3600;
  reatime = reatime%3600;
  int min = reatime/60;
  reatime = reatime%60;
  int sec = reatime;

  std::cout<<"!!!time in interval HH:MM:SS format is: "<<hour<<" h "
           <<min<<" m "<<sec<<" s\n";
}

template<class G>
void run (const G& g, const HeleshawDomain<G::dimension> & domain, Dune::ParameterTree & param)
{
  typedef typename G::LeafGridView GV;
  GV gv = g.leafGridView();

  // choose some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  const int dim = GV::dimension;
  Dune::Timer watch;
  double ttwophase(0), ttransport(0), treaction(0);

  typename Dune::Dycap::DycapTimer computationwatch;
  computationwatch.reset();


  // verbosity
  const int verbosity = param.sub("Verbosity").get<int>("verbosity", 2);

  std::string resultpath = param.get<std::string>("resultpath","");
  std::string cutpath = "cut";
  std::string vtkpath = "vtk";
  std::string solutionpath = "solution";
  if (!(resultpath==""))
    {
      cutpath = Dune::concatPaths(resultpath,cutpath);
      vtkpath = Dune::concatPaths(resultpath,vtkpath);
      solutionpath = Dune::concatPaths(resultpath,solutionpath);
      mkdir(resultpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

  if (gv.comm().rank()==0)
    {
      mkdir(vtkpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(cutpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(solutionpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

      std::ofstream s;
      s.open( Dune::concatPaths(resultpath,"param.conf").c_str());
      param.report(s);
      s.close();
    }
  gv.comm().barrier();


  // timeloop
  RF timestep = param.sub("Timeloop").get<RF>("timestep");
  RF subtimestepmax = param.sub("Timeloop").get<RF>("subtimestepmax");

  // at this gas saturation will be the concentration in gas phase zero
  const RF sgmin = param.sub("Setup").get<RF>("sgmin");

  // update timesteps according to refinement
  timestep /= 1<<domain.refine;

  /********************************************************************/
  /*********************** twophase problem ***************************/
  /********************************************************************/
  if (gv.comm().rank() == 0)
    std::cout << "=== multispan twophase parameter object\n";
  watch.reset();
  typedef MultiSpanTwoPhase<GV,RF> TimeManager;
  TimeManager timemanager(gv, param);

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

  //velocity for visualization
  typedef PoreVelocity<typename TPS::VELDGF,TP> PoreVelocity;
  PoreVelocity porevelocity(tps.getVl(),tp);

  /********************************************************************/
  /*********************** transport problem  **************************/
  /********************************************************************/

  // Make grid function space
  Dune::GeometryType gt;
  gt.makeCube(dim);
  // finite volumes finite element map
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(gt);
  typedef Dune::PDELab::P0ParallelConstraints CON;
  CON con;

  // grid function space
  typedef Dune::PDELab::istl::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> CGFS;

  if (gv.comm().rank() == 0)
    std::cout << "=== transport function space setup\n";
  watch.reset();
  CGFS cgfs(gv,fem,con);

  /********************/
  // for each phase are transport parameters different
  typedef TransportParameterBase<GV,RF> TPB;
  typedef TransportLiquid<GV,RF,TP,typename TPS::S_lDGF,typename TPS::S_lDGF,typename TPS::VELDGF> C1TP;
  typedef TransportGas<GV,RF,TP,typename TPS::S_gDGF,typename TPS::S_gDGF, typename TPS::VELDGF,typename TPS::P_gDGF,typename TPS::P_gDGF> C2TP;
  // transport parameters objects for each component (with different properties like diffusion etc.)
  TPB* tpb1 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportMicroorganism");
  TPB* tpb2 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportDOC");
  // TPB* tpb2 = new C2TP(tp,domain,tps.getS_gdgfOld(),tps.getS_gdgf(),tps.getVg(),tps.getP_gdgf(),tps.getP_gdgf(),"TransportOxygenGas");
  TPB* tpb3 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportOxygenWater");
  TPB* tpb4 = new C2TP(tp,domain,tps.getS_gdgfOld(),tps.getS_gdgf(),tps.getVg(),tps.getP_gdgfOld(),tps.getP_gdgf(),"TransportOxygenGas");
  TPB* tpb5 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVzero(),tps.getVzero(),"PoreMicroorganism");


  typedef typename Dune::Dycap::TransportSolver<GV,TPB,SECONDORDER> ETS;
  ETS ets1(gv,*tpb1,param);
  ETS ets2(gv,*tpb2,param);
  ETS ets3(gv,*tpb3,param);
  ETS ets4(gv,*tpb4,param,true);

  typedef Dune::PDELab::MulticomponentTransport<GV,RF,TPB,5> CTP;
  CTP ctp(tpb1,tpb2,tpb3,tpb4,tpb5); // ecoli,doc,o2g,o2w

  //last component for adhesion
  typedef typename ETS::CV CV;
  CV c5new(cgfs,0.0);

  // concentration grid functions
  typedef Dune::PDELab::DiscreteGridFunction<CGFS,CV> ConcDGF;
  ConcDGF c5dgf(cgfs,c5new);

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;
  /*********************/

  /********************************************************************/
  /*********************** reactive problem  **************************/
  /********************************************************************/
  if (gv.comm().rank() == 0)
    std::cout << "======= reactive function space setup  =======\n";
  watch.reset();

  // make reaction parameter object
  typedef ReactionParameterEcoli<GV,RF,TP,typename TPS::S_lDGF, typename TPS::VELDGF> RP;
  RP rp(gv, tp, tps.getS_ldgf(), tps.getS_ldgf(), tps.getVl());
  // setup ODE problem
  // exchange term
  typedef EcoliGrowthAdhesionModel<RP> ODEProblem;
  //typedef EcoliGrowthModel<RP> ODEProblem;
  ODEProblem growth(rp,param.sub("Microorganism").get<int>("verbosity",0));
  typedef EcoliGrowthModel<RP> ODEProblemOriginal;
  ODEProblemOriginal growth_original(rp,param.sub("Microorganism").get<int>("verbosity",0));

  typedef OdeBase<ODEProblem> ODESolver;
  std::shared_ptr<ODESolver> odesolver(setOdeSolver(growth,param.get<std::string>("odemethod","ExplicitEuler")));
  growth.setVerbosityLevel(verbosity);

  //grid operator for ode problem
  typedef EcoliGrowthGridOperator<RF,GV, ODEProblem, ODESolver> ODEGOS;
  ODEGOS odegos(gv,growth,*odesolver,param.get<int>("Verbosity.Reaction"));

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;


  /*******************************************************************/
  /************************** DGF OUTPUT  ****************************/
  /*******************************************************************/

  if (gv.comm().rank() == 0)
    std::cout << "=== absolute discrete grid function setup\n";
  watch.reset();

  const RF weight = param.sub("Microorganism").get<RF>("weight");
  const RF porosity = param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("porosity");



  typedef ProductGridFunction<GV,typename ETS::ConcDGF,typename TPS::S_lDGF> PoreConcentrationEcoli;
  PoreConcentrationEcoli poreEcoli_dgf(ets1.getConcDGF(),tps.getS_ldgf(),porosity/weight); // [g/l pore volume]


  typedef ProductGridFunction<GV,ConcDGF,typename TPS::S_lDGF> PoreConcentrationEcoliAdhesion;
  PoreConcentrationEcoliAdhesion poreEcoliadhesion_dgf(c5dgf,tps.getS_ldgf(),porosity/weight); // [g/l pore volume]

  typedef PlusGridFunction<GV,typename ETS::ConcDGF,ConcDGF> EcoliAbs;
  EcoliAbs ecoliabs_dgf(ets1.getConcDGF(),c5dgf);

  typedef ProductGridFunction<GV,EcoliAbs,typename TPS::S_lDGF> PoreConcentrationEcoliAbs;
  PoreConcentrationEcoliAbs poreEcoliabs_dgf(ecoliabs_dgf,tps.getS_ldgf(),porosity/weight); // [g/l pore volume]

  typedef ProductGridFunction<GV,typename ETS::ConcDGF,typename TPS::S_lDGF> AbsO2Water;
  AbsO2Water abso2water_dgf(ets3.getConcDGF(),tps.getS_ldgf());

  typedef ConstantDiscreteGridFunction<GV,RF> ConstDGF;
  ConstDGF molecularmassO2(gv,32.);

  typedef ProductGridFunction<GV,typename ETS::ConcDGF,ConstDGF> O2Air;
  O2Air o2air_dgf(ets4.getConcDGF(),molecularmassO2);

  typedef ProductGridFunction<GV,O2Air,typename TPS::S_gDGF> AbsO2Air;
  AbsO2Air abso2air_dgf(o2air_dgf,tps.getS_gdgf());

  typedef PlusGridFunction<GV,AbsO2Water,AbsO2Air> AbsO2;
  AbsO2 abso2_dgf(abso2water_dgf,abso2air_dgf);

  typedef O2<typename TPS::S_lDGF, typename ETS::ConcDGF, typename ETS::ConcDGF, RF> O2combination;
  O2combination o2combination(tps.getS_ldgf(), ets3.getConcDGF(), ets4.getConcDGF(), param.sub("PhaseExchange").get<RF>("K_H"));

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  /*******************************************************************/
  /********************* VTK OUTPUT for instationary problems ********/
  /*******************************************************************/
  RF watercontent = computeMass(tps.getS_ldgf());

  if (gv.comm().rank() == 0)
    std::cout << "=== VTK output setup\n";
  watch.reset();

  // graphics for initial value
  bool graphics = param.get<bool>("writeVTK",false);
  bool graphicsdebug = param.get<bool>("writeVTKdebug",false);
  bool graphics_subtime = param.get<bool>("writeVTKsubtime",false);

  // graphics for initial value
  std::string basename = param.get<std::string>("VTKname","")+ "refine_"+std::to_string(static_cast<long long int>(domain.refine));

  typedef Dune::PVDWriter<GV> PVDWriter;
  PVDWriter pvdwriter(gv,basename,Dune::VTK::conforming,Dune::VTK::appendedraw,vtkpath);
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename ETS::ConcDGF>>(ets1.getConcDGF(),"Ecoli [g/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoli>>(poreEcoli_dgf,"Ecoli [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename ETS::ConcDGF>>(ets2.getConcDGF(),"DOC [g/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename ETS::ConcDGF>>(ets3.getConcDGF(),"O2 [mg/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename ETS::ConcDGF>>(ets4.getConcDGF(),"O2 [mol/m3 air]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<O2combination>>(o2combination,"O2 [mg/l water] at equilibrium"));

  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<ConcDGF>>(c5dgf,"Ecoli adhesion [g/l water]"));
  pvdwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::VELDGF>>(tps.getVl(),"darcy flux velocity [m/s]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreVelocity>>(porevelocity,"pore velocity [m/d]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::VELDGF>>(tps.getVg(),"gas velocity"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::S_lDGF>>(tps.getS_ldgf(),"water saturation"));

  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoli>>(poreEcoli_dgf,"Ecoli (nonadhese) [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoliAdhesion>>(poreEcoliadhesion_dgf,"Ecoli (adhese) [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoliAbs>>(poreEcoliabs_dgf,"Ecoli (total) [cells/ml pore volume]"));

  ets4.setOutput(pvdwriter);

  using V0 = Dune::PDELab::Backend::Vector<CGFS,RF>;
  V0 partition(cgfs,gv.comm().rank());
  Dune::PDELab::AddDataHandle<CGFS,V0> pdh(cgfs,partition);
  if (cgfs.gridView().comm().size()>1)
    cgfs.gridView().communicate(pdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  typedef Dune::PDELab::DiscreteGridFunction<CGFS,V0> DGF0;
  DGF0 pdgf(cgfs,partition);
  pvdwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF0>>(pdgf,"processor decomposition"));


  // where should we made a cut
  RF cut1 =  param.get<RF>("cut1");
  RF cut2 =  param.get<RF>("cut2");
  RF cut3 =  param.get<RF>("cut3");

  std::string bbasenamecut1  = Dune::concatPaths(cutpath ,"cut1_");
  std::string bbasenamecut2  = Dune::concatPaths(cutpath ,"cut2_");
  std::string bbasenamecut3  = Dune::concatPaths(cutpath ,"cut3_");
  // class to make cut in the solution for arbitrary DGF
  GnuplotCutSolutionInTime<GV> bgnuplot1(gv,cut1,bbasenamecut1, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");
  GnuplotCutSolutionInTime<GV> bgnuplot2(gv,cut2,bbasenamecut2, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");
  GnuplotCutSolutionInTime<GV> bgnuplot3(gv,cut3,bbasenamecut3, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  /*******************************************************************/
  /************************** COMPUTATION  ***************************/
  /*******************************************************************/


  if (graphics)
    pvdwriter.write(0);

  typedef CFLTransportController<GV,TPB,typename ETS::FR> CFLController;
  CFLController cflcontroller(gv,*tpb3,ets3.getFlux() ,1.e-7,1);


  // VTK output every mod time steps
  //  int modGnuplot = param.get<int>("GNUPLOTcounter",20);
  int counter = 0;

  // get initial timestep
  timemanager.set_dt(timestep);
  RF tstep = timemanager.getTimeStepSize();

  if (gv.comm().rank() == 0)
    std::cout << "timestep init " << tstep << std::endl;


  int nrofspans=timemanager.getNumberOfSpans();

  std::vector<int> twostepnumber;
  std::vector<int> transportstepnumber;
  std::vector<double> twotime;
  std::vector<double> transporttime;
  std::vector<double> rtime;
  twostepnumber.resize(nrofspans);
  transportstepnumber.resize(nrofspans);
  twotime.resize(nrofspans);
  transporttime.resize(nrofspans);
  rtime.resize(nrofspans);


  int oldspannr=timemanager.getSpanNumber();
  while (!timemanager.finalize())
    {
      tstep = timemanager.getTimeStepSize();

      gv.comm().barrier();
      int spannr=timemanager.getSpanNumber();
      if (spannr>oldspannr)
        {
          twostepnumber[oldspannr]=tps.getCounter();
          tps.resetCounter();
          transportstepnumber[oldspannr]=ets4.getCounter();
          ets4.resetCounter();
        }
      gv.comm().barrier();

      if (gv.comm().rank() == 0)
        std::cout << "span number " << spannr << " from " << nrofspans << "\n";

      watch.reset();
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
              std::cout << "twophase problem solved" << std::endl;
          }
        else
          {
            timemanager.notifyFailure();
            if (gv.comm().rank() == 0)
              std::cout << "twophase problem NOT solved " << std::endl;
            continue;
          }
      }
      else if (timemanager.compute("computeZeroVelocity"))
        {
          tps.applyZeroVelocity(timemanager.getTime());
          if (gv.comm().rank() == 0)
            std::cout << "twophase with ZERO velocity" << std::endl;
        }
      ttwophase+=watch.elapsed();
      twotime[spannr]+=watch.elapsed();
      watch.reset();

      // advance concentrations (subtimesteps)
      RF subtime = timemanager.getTime();
      RF subtimetarget = subtime+tstep;

      ctp.setTime(subtime);
      ctp.setTimeTarget(subtimetarget,tstep);

      if (gv.comm().rank() == 0)
        std::cout << "compute cfl number" << std::endl;

      RF subtimestep = std::min(cflcontroller.suggestTimeStep()*param.sub("Setup").get<RF>("cfl")-1.e-8,subtimetarget-subtime);
      subtimestep=std::min(subtimestep,subtimestepmax);


      if (gv.comm().rank() == 0)
        std::cout << "subtimestep after cfl number is " << subtimestep << std::endl;

      // subtime loop (subtimesteps are smaller than timesteps for two phase flow)

      while (subtime<subtimetarget)
        {
          ctp.setTime(subtime);


          // control timestep
          RF timestepO2liquid = subtimestep;

          // ensure timestep constraints
          if (gv.comm().rank() == 0)
            std::cout << "transport subtime " << subtime << " with subtimestep " << subtimestep << std::endl;

          replaceSolutionUnderLimit(ets4.getGFS(),tps.getS_gdgf(),ets4.getSolution(),sgmin,0.0);
          replaceSolutionUnderLimit(ets4.getGFS(),ets4.getConcDGF(),ets4.getSolution(),0.0,0.0,1);
          ets4.update();

          watch.reset();
          // Ecoli transport in liquid phase
          if (timemanager.compute("computeTransportMicroorganism"))
            ets1.apply(subtime,subtimestep);

          // DOC transport in liquid phase
          if (timemanager.compute("computeTransportDOC"))
            ets2.apply(subtime,subtimestep);

          if (timemanager.compute("computeTransportOxygenWater"))
            ets3.apply(subtime,subtimestep);

          if (timemanager.compute("computeTransportOxygenGas"))
            ets4.apply(subtime,subtimestep);
          ttransport+=watch.elapsed();

          if (graphics_subtime)
            pvdwriter.write(subtime+subtimestep-1);

          ets1.update();

          transporttime[spannr]+=watch.elapsed();

          gv.comm().barrier();
          // compute microbiological growth
          if (timemanager.compute("computeReaction"))
            {


              RF reactiontimestep = subtimestep;
              RF reactiontime = subtime;

              RF reactiontimetarget = subtime+subtimestep;
              if (gv.comm().rank() == 0)
                std::cout << "======= Monod Model ======= from " << reactiontime << " to " << reactiontimetarget << "\n";
              watch.reset();

              std::vector<CV*> datavector;

              datavector.push_back(&ets1.getSolution()); // bacteriums
              datavector.push_back(&ets2.getSolution()); // DOC
              datavector.push_back(&ets3.getSolution()); // oxygen in liquid
              datavector.push_back(&ets4.getSolution()); // oxygen in air
              datavector.push_back(&c5new); // oxygen in air


              while (reactiontime<reactiontimetarget-1.e-9)
                {
                  odegos.apply(reactiontime,reactiontimestep,datavector);
                  reactiontime += reactiontimestep;
                  reactiontimestep = std::min(reactiontimestep, subtimestepmax);
                  reactiontimestep = std::min(reactiontimestep,reactiontimetarget-reactiontimestep);

                  if (graphics_subtime)
                    {
                      pvdwriter.write(reactiontime);
                      std::cout << "write subtimestep solution at " << reactiontime << std::endl;
                    }
                }
              if (gv.comm().rank() == 0)
                std::cout << "======= monod model finished ======= time : " << watch.elapsed() << " s" << std::endl;
              treaction+=watch.elapsed();
              rtime[spannr]+=watch.elapsed();
            }


          //for (auto i : ets1.getSolution())
          //  std::cout << i << " " << std::endl;

          if (graphicsdebug)
            pvdwriter.write(subtime+subtimestep);


          // update sub time step
          subtime += subtimestep;
          subtimestep=std::min(subtimestep,subtimestepmax);
          subtimestep = std::min(subtimestep, subtimetarget-subtime);

          //   pvdwriter.write(subtime);
          ctp.setTime(subtime);
        }

      if (gv.comm().rank()==0)
        std::cout << "times: twophase " << ttwophase << " transport " << ttransport << " treaction " << treaction << std::endl;


      // notify success in this timestep
      timemanager.notifySuccess(5);
      //tstep=timemanager.getTimeStepSize();

      if (timemanager.isTimeForOutput() && graphics)
        {
        pvdwriter.write(timemanager.getTime());
        bgnuplot1.output(timemanager.getTime(), tps.getS_ldgf(), ets1.getConcDGF(),ets2.getConcDGF(),ets3.getConcDGF(),ets4.getConcDGF(), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        bgnuplot2.output(timemanager.getTime(), tps.getS_ldgf(), ets1.getConcDGF(),ets2.getConcDGF(),ets3.getConcDGF(),ets4.getConcDGF(), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        bgnuplot3.output(timemanager.getTime(), tps.getS_ldgf(), ets1.getConcDGF(),ets2.getConcDGF(),ets3.getConcDGF(),ets4.getConcDGF(), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        }

      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep done, T= " << timemanager.getTime() << " ============= with ct " << computationwatch.elapsed() << " s \n\n";

      gv.comm().barrier();
      counter++;
      oldspannr=timemanager.getSpanNumber();
      twostepnumber[oldspannr]=tps.getCounter();
      transportstepnumber[oldspannr]=ets4.getCounter();

    }

  RF water_content = computeMass(tps.getS_ldgf());
  if (gv.comm().rank()==0)
    std::cout << "water content is " <<  water_content - watercontent << std::endl;
  if (gv.comm().rank()==0)
    std::cout << "twophase time step nr.  " << tps.getCounter() << " transport time step nr. " << ets2.getCounter() << std::endl;


  if (gv.comm().rank()==0)
    for (size_t i=0;i<rtime.size();++i)
    {
      RF reartime1=twotime[i];
      std::cout<<"\n!!phasetime " << i  << " is "<<reartime1  << " with " << twostepnumber[i] << std::endl;
      RF reartime2=transporttime[i];
      std::cout<<"\n!!transporttime " << i  << " is "<<reartime2 << " with " << transportstepnumber[i] << std::endl;
      RF reartime3=rtime[i];
      std::cout<<"\n!!reactiontime " << i  << " is "<<reartime3 << " with " << transportstepnumber[i] << std::endl;
      std::cout<<"\n!!reactiontime " << i  << " is "<<reartime2+reartime3 << " with "<< transportstepnumber[i] << std::endl;

      timeFormat(static_cast<int>(reartime1));
      timeFormat(static_cast<int>(reartime2));
      timeFormat(static_cast<int>(reartime3));
    }

}





template<class G>
void run_implicit(const G& g, const HeleshawDomain<G::dimension> & domain, Dune::ParameterTree & param)
{
  typedef typename G::LeafGridView GV;
  GV gv = g.leafGridView();

  // choose some types
  typedef typename GV::Grid::ctype DF;
  typedef double RF;
  //const int dim = GV::dimension;
  Dune::Timer watch;
  double ttwophase(0), ttransport(0), treaction(0);

  typename Dune::Dycap::DycapTimer computationwatch;
  computationwatch.reset();


  // verbosity
  // const int verbosity = param.sub("Verbosity").get<int>("verbosity", 2);

  std::string resultpath = param.get<std::string>("resultpath","");
  std::string cutpath = "cut";
  std::string vtkpath = "vtk";
  std::string solutionpath = "solution";
  if (!(resultpath==""))
    {
      cutpath = Dune::concatPaths(resultpath,cutpath);
      vtkpath = Dune::concatPaths(resultpath,vtkpath);
      solutionpath = Dune::concatPaths(resultpath,solutionpath);
      mkdir(resultpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

  if (gv.comm().rank()==0)
    {
      mkdir(vtkpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(cutpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(solutionpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

      std::ofstream s;
      s.open( Dune::concatPaths(resultpath,"param.conf").c_str());
      param.report(s);
      s.close();
    }
  gv.comm().barrier();


  // timeloop
  RF timestep = param.sub("Timeloop").get<RF>("timestep");
  RF subtimestepmax = param.sub("Timeloop").get<RF>("subtimestepmax");

  // at this gas saturation will be the concentration in gas phase zero
  const RF sgmin = param.sub("Setup").get<RF>("sgmin");

  // update timesteps according to refinement
  timestep /= 1<<domain.refine;

  /********************************************************************/
  /*********************** twophase problem ***************************/
  /********************************************************************/
  if (gv.comm().rank() == 0)
    std::cout << "=== multispan twophase parameter object\n";
  watch.reset();
  typedef MultiSpanTwoPhase<GV,RF> TimeManager;
  TimeManager timemanager(gv, param);

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

  //velocity for visualization
  typedef PoreVelocity<typename TPS::VELDGF,TP> PoreVelocity;
  PoreVelocity porevelocity(tps.getVl(),tp);


  /********************/
  // for each phase are transport parameters different
  typedef TransportParameterBase<GV,RF> TPB;
  typedef TransportLiquid<GV,RF,TP,typename TPS::S_lDGF,typename TPS::S_lDGF,typename TPS::VELDGF> C1TP;
  typedef TransportGas<GV,RF,TP,typename TPS::S_gDGF,typename TPS::S_gDGF, typename TPS::VELDGF,typename TPS::P_gDGF,typename TPS::P_gDGF> C2TP;
  // transport parameters objects for each component (with different properties like diffusion etc.)
  TPB* tpb1 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportMicroorganism");
  TPB* tpb2 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportDOC");
  // TPB* tpb2 = new C2TP(tp,domain,tps.getS_gdgfOld(),tps.getS_gdgf(),tps.getVg(),tps.getP_gdgf(),tps.getP_gdgf(),"TransportOxygenGas");
  TPB* tpb3 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVlOld(),tps.getVl(),"TransportOxygenWater");
  TPB* tpb4 = new C2TP(tp,domain,tps.getS_gdgfOld(),tps.getS_gdgf(),tps.getVg(),tps.getP_gdgfOld(),tps.getP_gdgf(),"TransportOxygenGas");
  TPB* tpb5 = new C1TP(tp,domain,tps.getS_ldgfOld(),tps.getS_ldgf(),tps.getVzero(),tps.getVzero(),"PoreMicroorganism");


  typedef Dune::PDELab::MulticomponentTransport<GV,RF,TPB,5> CTP;
  CTP ctp(tpb1,tpb2,tpb3,tpb4,tpb5); // ecoli,doc,o2g,o2w


  /********************************************************************/
  /*********************** reactive problem  **************************/
  /********************************************************************/
  if (gv.comm().rank() == 0)
    std::cout << "======= reactive function space setup  =======\n";
  watch.reset();

  // make reaction parameter object
  typedef ReactionParameterEcoli<GV,RF,TP,typename TPS::S_lDGF, typename TPS::VELDGF> RP;
  RP rp(gv, tp, tps.getS_ldgf(), tps.getS_ldgf(), tps.getVl());
  // setup ODE problem
  // exchange term
  typedef EcoliGrowthAdhesionModel<RP> ODEProblem;
  ODEProblem growth(rp,param.sub("Microorganism").get<int>("verbosity",0));

   if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  typedef ReactionAdapter<ODEProblem> RAdapter;
  RAdapter radapter(growth);

  //typedef ReactionBaseAdapter RAdapter;
  //RAdapter radapter;

  typedef Dune::Dycap::MulticomponentImplicitTransportSolver<RF,GV,CTP,RAdapter> MCTS;
  MCTS mcts(gv,ctp,radapter,param,0,2);
  typedef typename MCTS::DGF DGF;

  /*******************************************************************/
  /************************** DGF OUTPUT  ****************************/
  /*******************************************************************/

  if (gv.comm().rank() == 0)
    std::cout << "=== absolute discrete grid function setup\n";
  watch.reset();

  const RF weight = param.sub("Microorganism").get<RF>("weight");
  const RF porosity = param.sub(param.sub("Setup").get<std::string>("material")).get<RF>("porosity");



  typedef ProductGridFunction<GV,DGF,typename TPS::S_lDGF> PoreConcentrationEcoli;
  PoreConcentrationEcoli poreEcoli_dgf(mcts.getConcDGF(0),tps.getS_ldgf(),porosity/weight); // [g/l pore volume]


  typedef ProductGridFunction<GV,DGF,typename TPS::S_lDGF> PoreConcentrationEcoliAdhesion;
  PoreConcentrationEcoliAdhesion poreEcoliadhesion_dgf(mcts.getConcDGF(4),tps.getS_ldgf(),porosity/weight); // [g/l pore volume]

  typedef PlusGridFunction<GV,DGF,DGF> EcoliAbs;
  EcoliAbs ecoliabs_dgf(mcts.getConcDGF(0),mcts.getConcDGF(4));

  typedef ProductGridFunction<GV,EcoliAbs,typename TPS::S_lDGF> PoreConcentrationEcoliAbs;
  PoreConcentrationEcoliAbs poreEcoliabs_dgf(ecoliabs_dgf,tps.getS_ldgf(),porosity/weight); // [g/l pore volume]

  typedef ProductGridFunction<GV,DGF,typename TPS::S_lDGF> AbsO2Water;
  AbsO2Water abso2water_dgf(mcts.getConcDGF(2),tps.getS_ldgf());

  typedef ConstantDiscreteGridFunction<GV,RF> ConstDGF;
  ConstDGF molecularmassO2(gv,32.);

  typedef ProductGridFunction<GV,DGF,ConstDGF> O2Air;
  O2Air o2air_dgf(mcts.getConcDGF(3),molecularmassO2);

  typedef ProductGridFunction<GV,O2Air,typename TPS::S_gDGF> AbsO2Air;
  AbsO2Air abso2air_dgf(o2air_dgf,tps.getS_gdgf());

  typedef PlusGridFunction<GV,AbsO2Water,AbsO2Air> AbsO2;
  AbsO2 abso2_dgf(abso2water_dgf,abso2air_dgf);

  typedef O2<typename TPS::S_lDGF, DGF, DGF, RF> O2combination;
  O2combination o2combination(tps.getS_ldgf(), mcts.getConcDGF(2), mcts.getConcDGF(3), param.sub("PhaseExchange").get<RF>("K_H"));

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  /*******************************************************************/
  /********************* VTK OUTPUT for instationary problems ********/
  /*******************************************************************/
  RF watercontent = computeMass(tps.getS_ldgf());

  if (gv.comm().rank() == 0)
    std::cout << "=== VTK output setup\n";
  watch.reset();

  // graphics for initial value
  bool graphics = param.get<bool>("writeVTK",false);

  // graphics for initial value
  std::string basename = param.get<std::string>("VTKname","")+ "refine_"+std::to_string(static_cast<long long int>(domain.refine));

  typedef Dune::PVDWriter<GV> PVDWriter;
  PVDWriter pvdwriter(gv,basename,Dune::VTK::conforming,Dune::VTK::appendedraw,vtkpath);
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(mcts.getConcDGF(0),"Ecoli [g/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoli>>(poreEcoli_dgf,"Ecoli [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(mcts.getConcDGF(1),"DOC [g/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(mcts.getConcDGF(2),"O2 [mg/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(mcts.getConcDGF(3),"O2 [mol/m3 air]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<O2combination>>(o2combination,"O2 [mg/l water] at equilibrium"));


  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(mcts.getConcDGF(4),"Ecoli adhesion [g/l water]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::VELDGF>>(tps.getVl(),"darcy flux velocity [m/s]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreVelocity>>(porevelocity,"pore velocity [m/d]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::VELDGF>>(tps.getVg(),"gas velocity"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<typename TPS::S_lDGF>>(tps.getS_ldgf(),"water saturation"));

  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoli>>(poreEcoli_dgf,"Ecoli (nonadhese) [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoliAdhesion>>(poreEcoliadhesion_dgf,"Ecoli (adhese) [cells/ml pore volume]"));
  pvdwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<PoreConcentrationEcoliAbs>>(poreEcoliabs_dgf,"Ecoli (total) [cells/ml pore volume]"));

  /*
  typedef typename Dune::PDELab::BackendVectorSelector<CGFS,RF>::Type V0;
  V0 partition(cgfs,gv.comm().rank());
  Dune::PDELab::AddDataHandle<CGFS,V0> pdh(cgfs,partition);
  if (cgfs.gridView().comm().size()>1)
    cgfs.gridView().communicate(pdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  typedef Dune::PDELab::DiscreteGridFunction<CGFS,V0> DGF0;
  DGF0 pdgf(cgfs,partition);
  pvdwriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF0>(pdgf,"processor decomposition"));
  */

 // where should we made a cut
  RF cut1 =  param.get<RF>("cut1");
  RF cut2 =  param.get<RF>("cut2");
  RF cut3 =  param.get<RF>("cut3");

  std::string bbasenamecut1  = Dune::concatPaths(cutpath ,"cut1_");
  std::string bbasenamecut2  = Dune::concatPaths(cutpath ,"cut2_");
  std::string bbasenamecut3  = Dune::concatPaths(cutpath ,"cut3_");
  // class to make cut in the solution for arbitrary DGF
  GnuplotCutSolutionInTime<GV> bgnuplot1(gv,cut1,bbasenamecut1, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");
  GnuplotCutSolutionInTime<GV> bgnuplot2(gv,cut2,bbasenamecut2, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");
  GnuplotCutSolutionInTime<GV> bgnuplot3(gv,cut3,bbasenamecut3, "h, Sl, ecoli, doc[g/l], o2water[mg/l], o2air[mol/m3], ecoli [cell/cm^3 pore volume], o2air[mg/l], o2air[mg/l liquid]");

  if (gv.comm().rank() == 0)
    std::cout << "... took : " << watch.elapsed() << " s" << std::endl;

  /*******************************************************************/
  /************************** COMPUTATION  ***************************/
  /*******************************************************************/


  if (graphics)
    pvdwriter.write(0);

  typedef typename Dune::Dycap::TransportSolver<GV,TPB,false> ETS;
  ETS ets3(gv,*tpb3,param);
  typedef CFLTransportController<GV,TPB,typename ETS::FR> CFLController;
  CFLController cflcontroller(gv,*tpb3,ets3.getFlux() ,1.e-7,1);


  // VTK output every mod time steps
  //  int modGnuplot = param.get<int>("GNUPLOTcounter",20);
  int counter = 0;

  // get initial timestep
  timemanager.set_dt(timestep);
  RF tstep = timemanager.getTimeStepSize();

  if (gv.comm().rank() == 0)
    std::cout << "timestep init " << tstep << std::endl;


  int nrofspans=timemanager.getNumberOfSpans();

  std::vector<int> twostepnumber;
  std::vector<int> transportstepnumber;
  std::vector<double> twotime;
  std::vector<double> transporttime;
  twostepnumber.resize(nrofspans);
  transportstepnumber.resize(nrofspans);
  twotime.resize(nrofspans);
  transporttime.resize(nrofspans);

  int oldspannr=timemanager.getSpanNumber();


  while (!timemanager.finalize())
    {
      tstep = timemanager.getTimeStepSize();

      gv.comm().barrier();
      int spannr=timemanager.getSpanNumber();
      if (spannr>oldspannr)
        {
          twostepnumber[oldspannr]=tps.getCounter();
          tps.resetCounter();
          transportstepnumber[oldspannr]=mcts.getCounter();
          mcts.resetCounter();
        }
      gv.comm().barrier();

      if (gv.comm().rank() == 0)
        std::cout << "span number " << spannr << " from " << nrofspans << "\n";

      watch.reset();
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
              std::cout << "twophase problem solved" << std::endl;
          }
        else
          {
            timemanager.notifyFailure();
            if (gv.comm().rank() == 0)
              std::cout << "twophase problem NOT solved " << std::endl;
            continue;
          }
      }

      else if (timemanager.compute("computeZeroVelocity"))
        {
          tps.applyZeroVelocity(timemanager.getTime());
          if (gv.comm().rank() == 0)
            std::cout << "twophase with ZERO velocity" << std::endl;
        }
      ttwophase+=watch.elapsed();
      twotime[spannr]+=watch.elapsed();
      watch.reset();

      if (timemanager.compute("computeFullImplicit"))
        mcts.apply(timemanager.getTime(),tstep);
      treaction+=watch.elapsed();
      transporttime[spannr]+=watch.elapsed();

      if (gv.comm().rank()==0)
        std::cout << "times: twophase " << ttwophase << " transport " << ttransport << " treaction " << treaction << std::endl;

      // notify success in this timestep
      timemanager.notifySuccess(5);

      if (timemanager.isTimeForOutput() && graphics)
        {
        pvdwriter.write(timemanager.getTime());
        bgnuplot1.output(timemanager.getTime(), tps.getS_ldgf(), mcts.getConcDGF(0),mcts.getConcDGF(1),mcts.getConcDGF(2),mcts.getConcDGF(3), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        bgnuplot2.output(timemanager.getTime(), tps.getS_ldgf(), mcts.getConcDGF(0),mcts.getConcDGF(1),mcts.getConcDGF(2),mcts.getConcDGF(3), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        bgnuplot3.output(timemanager.getTime(), tps.getS_ldgf(), mcts.getConcDGF(0),mcts.getConcDGF(1),mcts.getConcDGF(2),mcts.getConcDGF(3), poreEcoli_dgf, o2air_dgf,abso2_dgf,poreEcoli_dgf, poreEcoliadhesion_dgf, poreEcoliabs_dgf);
        }

      if (gv.comm().rank() == 0)
        std::cout << "============= Timestep done, T= " << timemanager.getTime() << " ============= with ct " << computationwatch.elapsed() << " s \n\n";
      counter++;

      gv.comm().barrier();
      oldspannr=timemanager.getSpanNumber();
      twostepnumber[oldspannr]=tps.getCounter();
      transportstepnumber[oldspannr]=mcts.getCounter();
    }

  RF water_content = computeMass(tps.getS_ldgf());
  if (gv.comm().rank()==0)
    std::cout << "water content is " <<  water_content - watercontent << std::endl;
  if (gv.comm().rank()==0)
     std::cout << "twophase time step nr.  " << tps.getCounter() << " transport time step nr. " << mcts.getCounter() << std::endl;



  if (gv.comm().rank()==0)
    for (size_t i=0;i<twotime.size();++i)
    {
      RF reartime1=twotime[i];
      std::cout<<"\n!!phasetime " << i  << " is "<<reartime1  << " with " << twostepnumber[i] << std::endl;
      RF reartime2=transporttime[i];
      std::cout<<"\n!!transporttime " << i  << " is "<<reartime2 << " with " << transportstepnumber[i] << std::endl;

      timeFormat(static_cast<int>(reartime1));
      timeFormat(static_cast<int>(reartime2));

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
    std::string configfile = argv[0]; configfile += ".conf";

    // parse cmd line parameters
    if (argc > 1 && argv[1][0] != '-')
      configfile = argv[1];
    for (int i = 1; i < argc; i++)
      {
        if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
          {
            if(helper.rank()==0)
              std::cout << "usage: ./simulation <configfile> [OPTIONS]" << std::endl;
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
        int overlap = param.get<int>("overlap");
        Dune::YaspGrid<2> grid(L,N,periodic,overlap,Dune::MPIHelper::getCollectiveCommunication());
        grid.globalRefine(domain.refine);
        // solve problem
        bool fullimplicit = param.get<bool>("fullimplicit");
        if (fullimplicit)
          run_implicit(grid,domain,param);
          else
        run(grid,domain,param);
      }

    if(helper.rank()==0)
      std::cout<<"!!Total computation time was "<<MPI_Wtime()-start
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
