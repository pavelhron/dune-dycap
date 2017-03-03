// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TIMESTEPOUTPUT_HH
#define DUNE_DYCAP_TIMESTEPOUTPUT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iomanip>
#include"gnuplot.hh"
#include<dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

class TimeStepHelper
{
public:
  TimeStepHelper()
  {
  }

  inline int get_newton_iteration()
  {
    return 0;
  }

};

/**
 * \brief Writes timesteps for twophase problem in a txt file.
 *
 * Can be used as a control for the computation effort of twophase problem.
 * Because the greatest problem is bad Newton convergence and therefore
 * small times steps, this class allows to compare arbitrary methods.
 *
 * \tparam RF      C++ type of the floating point parameters
 * \tparam Helper  Helper class which should store the number of newton iteration in each step.
 **/
template <class RF, class Helper=TimeStepHelper>
class TimeStepWriter
{

public:
  TimeStepWriter(std::string dataFile_, int rank_, Helper & helper_, const bool pdfoutput_=false):
    dataFile(dataFile_), rank(rank_),helper(Dune::stackobject_to_shared_ptr(helper_)), pdfoutput(pdfoutput_)
  {
    open_file();
  }

  TimeStepWriter(std::string dataFile_, int rank_, const bool pdfoutput_=false):
    dataFile(dataFile_), rank(rank_), pdfoutput(pdfoutput_)
  {
    helper = Dune::make_shared<Helper>();
    open_file();
  }

  void open_file()
  {
    dataFile.append("_timesteps.txt");
    if (rank==0){
      s.open(dataFile.c_str());
      if (!s.is_open())
        {
          DUNE_THROW(Dune::Exception, dataFile + " not opened");
        }

      s << std::setprecision(14) << std::scientific;
    }
  }

  static TimeStepHelper & getHelper()
  {
    TimeStepHelper timestephelper;
    return timestephelper;
  }

  ~TimeStepWriter()
  {
    if (rank==0)
      {
        s.close();
        GNUplot plotter;

        std::ostringstream os;
        std::string soutput = dataFile;
        soutput.replace(soutput.length()-3,soutput.length(),"eps");

        if (pdfoutput) {
          os << "set terminal postscript eps color \"Helvetica\" 20 ;\n"
             << "set output \"" <<  soutput << "\" ; \n" ;
        }

        os << "set y2tics 2 nomirror ;\n"
           << "set ytics nomirror ;\n"
           << "set grid xtics y2tics nopolar;\n"
           << "set offsets 0,0,100,100;\n"
           << "solution_output=\"" << dataFile.c_str() << "\";\n"
           << "plot solution_output u 1:2 w l ti \"timesteps\", solution_output u 1:5 w p axes x1y2 ti \"NI\", solution_output u 1:3 w l ti \"min\", solution_output u 1:4 w l ti \"max\";\n ";

        plotter(os.str());
        os.clear();
      }
  }

  //! Add new timesteps
  void add_timestep(RF time, RF timestep, RF dtmin, RF dtmax)
  {
    if (rank==0)
      s << time << "  " << timestep << "  " << dtmin << "  " << dtmax << " " << (*helper).get_newton_iteration()<< "\n";
  }

private:
  std::string dataFile;
  int rank;
  Dune::shared_ptr<Helper> helper;
  std::ofstream s;
  const bool pdfoutput;
};

#endif
