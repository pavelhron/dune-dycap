// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_VTK_PVDWRITER_DEBUG_HH
#define DUNE_VTK_PVDWRITER_DEBUG_HH

#include <vector>
#include <fstream>
#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/pdelab/instationary/onestep.hh>

#include <stdlib.h>

namespace Dune
{

  template< class GridView, class VTK = VTKWriter<GridView> >
  class PVDWriterDebug : public VTK
  {
    typedef std::vector<unsigned int> Int1D;
    typedef std::vector<Int1D> Int2D;
    const GridView & gv;
    std::string basename;
    bool output;
    PDELab::FilenameHelper fn;
    std::string path;
    Int2D steps;

    Dune::VTK::OutputType outputtype;

  public:

    PVDWriterDebug(const GridView & gv_, std::string basename_, bool output_ = true,
	      Dune::VTK::DataMode datamode_ = Dune::VTK::conforming,
	      Dune::VTK::OutputType outputtype_ = Dune::VTK::appendedraw) :
      VTK(gv_,datamode_), gv(gv_),
      basename(basename_), output(output_), fn(basename_),
      path("vtk"),
      steps( Int2D(0, Int1D(0,0))),
      outputtype(outputtype_) {}

    void set_time(int time)
    {
      Int1D  ts(Int1D(0));
      steps.push_back(ts);
      //  std::cout << "set time " << time << " stepsize " << steps.size() << std::endl;
      set_iteration(0);
    }

    void set_iteration(int iteration)
    {
      steps.back().push_back(iteration);
      //  std::cout << "iteration " << iteration << " iteration vector length " << steps.back().size()<< std::endl;
    }

    void set_output(bool output_)
    {
      output = output_;
    }


    void write(int ls)
    {
      if (!output)
        return;
      /* remember current time step */
      steps.back().back()=ls;
      //  std::cout << "line search " << ls << std::endl;
      /* make sure the directory exists */
      // mkdir("vtk", 777);
      /* write VTK file */
      VTK::pwrite(fn.getName(),path,"",outputtype);
      /* write pvd file */
      std::string pvdname = basename + ".pvd";
      std::ofstream pvd(pvdname.c_str());
      if (gv.comm().rank() == 0)
        std::cout << "WRITE PVD FILE " << pvdname << " " << (steps.size()-1)*10000+((steps[steps.size()-1].size()-1)*100)+ls<< std::endl;
      assert(pvd.is_open());
      pvd << std::fixed;
      pvd << "<?xml version=\"1.0\"?>\n"
	  << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
	  << "<Collection>\n";
      PDELab::FilenameHelper fnloop(basename);
      for (unsigned int i=0; i<steps.size(); i++)
	{
          for (unsigned int j=0; j<steps[i].size();++j)
	    {
              for (unsigned int k=0; k<=steps[i][j];++k)
                {
                  std::string fname = this->getParallelHeaderName(fnloop.getName(), path, gv.comm().size());
	      pvd << "  <DataSet timestep=\"" << 10000*i+100*j+k
		  << "\" file=\"" << fname << "\"/>\n";
                  fnloop.increment();
                }
            }
	}
      pvd << "</Collection>\n"
	  << "</VTKFile>\n";
      pvd.close();

      /* increment counter */
       fn.increment();
    }
  };

} // end namespace Dune

#endif // DUNE_VTK_PVDWRITER_HH
