// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<iostream>
#include <iomanip>
#include<vector>
#include<map>
#include<gsl/gsl_sf_gamma.h>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/timer.hh>
#include<dune/common/parametertreeparser.hh>
#include<src/physics/hydraulicparameters.hh>
#include<src/physics/physical_chemistry.hh>

template<typename RF>
class TwoPhaseParameter
{

public:

  //! constructor
  TwoPhaseParameter(const Dune::ParameterTree & param_) :

    param(param_),  // parameter class
    material(param.sub("Setup").get<std::string>("material")), // material to use
    porosity(param.sub(material).get<RF>("porosity"))
  {

    // hydraulic parameters
    std::string model = param.sub(material).get<std::string>("model");
    bool linear_interpolation = param.sub(material).get<bool>("linearinterpolation",false);
    bool cubic_interpolation = param.sub(material).get<bool>("cubicinterpolation",false);

    if (linear_interpolation && cubic_interpolation)
      DUNE_THROW(Dune::Exception, "Only one interpolation can be used");
    if (model == "regBrooksCorey")
      {
        Dune::PM::RegBrooksCoreyParam<RF>* brooksCorey = new Dune::PM::RegBrooksCoreyParam<RF>;
        hydrParam = brooksCorey;
        brooksCorey->Setdelta(param.sub(material).get<RF>("delta",0.) );
        brooksCorey->SetLambda( param.sub(material).get<RF>("lambda",0.) );
        brooksCorey->SetTau( param.sub(material).get<RF>("tau", 0.5) );
        brooksCorey->SetEntryPressure( param.sub(material).get<RF>("pentry") );
        brooksCorey->InitRegularisation(param.sub(material).get<RF>("deltareg", 1e-3) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "BrooksCorey")
      {
        Dune::PM::BrooksCoreyParam<RF>* brooksCorey = new Dune::PM::BrooksCoreyParam<RF>;
        hydrParam = brooksCorey;
        brooksCorey->Setdelta(param.sub(material).get<RF>("delta",0.) );
        brooksCorey->SetLambda( param.sub(material).get<RF>("lambda",0.) );
        brooksCorey->SetTau( param.sub(material).get<RF>("tau", 0.5) );
        brooksCorey->SetEntryPressure( param.sub(material).get<RF>("pentry") );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "vanGenuchten")
      {
        Dune::PM::VanGenuchtenParam<RF>* vanGenuchten = new Dune::PM::VanGenuchtenParam<RF>;
        hydrParam = vanGenuchten;
        vanGenuchten->Setdelta(param.sub(material).get<RF>("delta",0.) );
        vanGenuchten->SetAlpha( param.sub(material).get<RF>("alpha") );
        vanGenuchten->SetN( param.sub(material).get<RF>("n") );
        if (param.sub(material).hasKey("m"))
          vanGenuchten->SetM( param.sub(material).get<RF>("m") );
        vanGenuchten->SetTau( param.sub(material).get<RF>("tau", 0.5) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "modVanGenuchten")
      {
        Dune::PM::ModVanGenuchtenParam<RF>* vanGenuchten = new Dune::PM::ModVanGenuchtenParam<RF>;
        hydrParam = vanGenuchten;
        vanGenuchten->Setdelta(param.sub(material).get<RF>("delta",0.) );
        vanGenuchten->SetAlpha( param.sub(material).get<RF>("alpha") );
        vanGenuchten->SetN( param.sub(material).get<RF>("n") );
        if (param.sub(material).hasKey("m"))
          vanGenuchten->SetM( param.sub(material).get<RF>("m") );
        vanGenuchten->SetTau( param.sub(material).get<RF>("tau", 0.5) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }
      }
    else
      DUNE_THROW(Dune::Exception, "Unknown hydraulic model " + model);
    assert(hydrParam);
  }


  //! porosity
  RF phi () const
  {
    return porosity;
  }

  RF pc (RF s_l) const
  {
    return hydrParam->Pc(s_l);
  }

  //! inverse capillary pressure function
  RF s_l (RF pc) const
  {
    return hydrParam->Sl(pc);
  }

  //! liquid surface tension
  RF sigma_l () const
  {
    return Dune::PM::PhysicalChemistry<RF>::RelSurfaceTensionWaterAir(293.16);
  }


  //! liquid phase relative permeability
  RF kr_l (RF s_l) const
  {
    return hydrParam->KRelL(s_l);
  }

  inline const Dune::ParameterTree & getParam() const
  {
    return param;
  }

  const Dune::PM::HydrParamBase<RF>* getHydrParam() const
  {
    return hydrParam;
  }


  const Dune::ParameterTree & param;
  const std::string material;


private:
  const RF porosity;
  Dune::PM::HydrParamBase<RF>* hydrParam;

};

template<typename RF, typename TP>
class PhaseExchangeModels{

  TP & tp;
  const Dune::ParameterTree & param;
  std::string fname;
  std::string comment;
  std::ofstream s;


public:
  //! Construct
  /**
   * \param fname     Name of the file to write to.
   * \param gv_       Grid View used for the vectors passed in.

   * \note The filename passed should be the same on all ranks.  Only rank
   *       0 will actually write to the file.
   */
  PhaseExchangeModels(TP & tp_, const std::string fname_, const std::string comment_="") :
    tp(tp_), param(tp.getParam()), fname(fname_), comment(comment_)
  {
  }

  //! Write a record in the output file (Element index)
  /**
   * \param t        time
   * \param x_coord  x-coordinate
   * \param ... dgf  discrete grid function ellipses
   */
  void output() {

    std::string sname = fname;
    s.exceptions(std::ios_base::badbit | std::ios_base::eofbit
                   | std::ios_base::failbit);
      s.open(sname.c_str());
      s << std::setprecision(14) << std::scientific;
      if(comment.size() > 0)
        s << "# " << comment << "\n";

      std::vector<std::string> models = {"niemet", "geistlinger"};
      const RF md = param.sub("PhaseExchange").get<RF>("md");
      const RF kappa = param.sub("PhaseExchange").get<RF>("kappa");

      const std::string material(param.sub("Setup").get<std::string>("material"));

      // alpha = alpha' / 100 * rho * g
      const RF alpha =  param.sub(material).get<RF>("alpha");
        const RF n = param.sub(material).get<RF>("n");
        RF m;
        if (param.sub(material).hasKey("m"))
          m =  param.sub(material).get<RF>("m");
        else m = 1.-1./n;

      const RF w = m-1./n;
      const RF z = 1.+1./n;

      std::cout << "pc " << tp.pc(1.) << " " << tp.pc(0.01)<< std::endl;

      // loop once over the grid and find elements with given x_coord
    for (RF sl = 1; sl>.0; sl-=0.0001)
      {

        RF agw;
        RF sg = 1.0-sl;
        s << sl << "\t";

        RF phi = tp.phi();
        RF sigma = tp.sigma_l();
        RF pc=tp.pc(sl);

        // Gvirtzmann, Roberts 1991 (modified version)
        agw = 6*(1-phi)/md*sg*kappa;
        s << agw  << "\t";

        // Niemet 2002
        const RF u = std::pow(sl,1/m);
        agw = 3*phi/(2.*sigma)*m/alpha *gsl_sf_beta(w,z)*(1.-gsl_sf_beta_inc(w,z,u));
        s << agw  << "\t";

        // Miller et al. (1990) S = pc phi / sigma
        agw = pc*phi/sigma*sg;
        s << agw  << "\t";

        // Joekar-Niasar 2008
        agw = (234. + 3858.*sl - 0.224*pc - 3992*sl*sl + 0.006*sl*pc + 1.283e-5*pc*pc);
        s << agw  << "\t";

        // Geistlinger
        agw = 6*phi/md*sg*kappa;
        s << agw  << "\t";

        // quadratic
        agw = (1699.4 - 1969.8*sl - 1.4*pc +274*sl*sl + 1.219*sl*pc +0.000289799*pc*pc);
        s << agw  << "\t";

        s << "\n";
      }

    s.close();
  }
};


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
              std::cout << "usage: ./twophase_onecomponent <configfile> [OPTIONS]" << std::endl;
            return 0;
          }
      }
    double start=MPI_Wtime();

    // read parameters from file
    Dune::ParameterTree param;
    Dune::ParameterTreeParser parser;
    parser.readINITree(configfile, param);

    typedef TwoPhaseParameter<double> TP;
    TP tp(param);

    std::string output = "result.dat";
    PhaseExchangeModels<double, TP> phaseexchange(tp,output);
    phaseexchange.output();

    int r = system("gnuplot exchange.gpl");


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
