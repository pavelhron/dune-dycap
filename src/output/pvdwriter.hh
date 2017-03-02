// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_VTK_PVDWRITER_HH
#define DUNE_VTK_PVDWRITER_HH

#include <vector>
#include <fstream>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

namespace Dune
{
    class FilenameHelper
    {
    public:
        FilenameHelper(const char *basename_, int i_=0)
            : i(i_)
        {
            sprintf(basename,"%s",basename_);
        }

        FilenameHelper(const std::string & basename_, int i_=0)
            : i(i_)
        {
            sprintf(basename,"%s",basename_.c_str());
        }

        const char *getName (int i_)
        {
            sprintf(fname,"%s-%05d",basename,i_);
            return fname;
        }

        const char *getName ()
        {
            sprintf(fname,"%s-%05d",basename,i);
            return fname;
        }

        void setName(const std::string & basename_)
        {
            sprintf(basename,"%s",basename_.c_str());
        }

        void increment ()
        {
            i++;
        }

    private:
        char fname[255];
        char basename[255];
        int i;
    };


    template< class GridView, class VTK = VTKWriter<GridView> >
    class PVDWriter : public VTK
    {
        const GridView & gv;
        std::string basename;
        FilenameHelper fn;
        std::string path;
        std::vector<double> timesteps;
        Dune::VTK::OutputType outputtype;
        unsigned int offset;

    public:

        PVDWriter(const GridView & gv_, unsigned int level_, std::string basename_,
                  Dune::VTK::DataMode datamode_ = Dune::VTK::conforming,
                  Dune::VTK::OutputType outputtype_ = Dune::VTK::appendedraw,
                  std::string path_="vtk", unsigned int offset_=0) :
            VTK(gv_,level_), gv(gv_),
            basename(basename_), fn(basename_,offset_),
            path(path_), outputtype(outputtype_),
            offset(offset_){}

        PVDWriter(const GridView & gv_, std::string basename_,
                  Dune::VTK::DataMode datamode_ = Dune::VTK::conforming,
                  Dune::VTK::OutputType outputtype_ = Dune::VTK::appendedraw,
                  std::string path_="vtk", unsigned int offset_=0) :
            VTK(gv_,datamode_), gv(gv_),
            basename(basename_), fn(basename_,offset_),
            path(path_), outputtype(outputtype_),
            offset(offset_){}

        //it works
        void write(double time)
        {
            /* remember current time step */
            timesteps.push_back(time);

            unsigned found = path.find_last_of("/\\");
            std::string pres = "";
            if (found)
                pres = path.substr(0,found);

            if (path.find("/") == std::string::npos)
                pres = "";

                //std::string pos = path.substr(found+1);

            //fn.setName(pos);

            // fn.setName(pres);
            /* make sure the directory exists */
            mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            /* write VTK file */
            VTK::pwrite(fn.getName(),path,"",outputtype);
            /* write pvd file */

            //std::cout << pres  << " " << basename << std::endl;
            std::string pvdname =  Dune::concatPaths(pres,basename) + ".pvd";
            std::ofstream pvd(pvdname.c_str());
            if (gv.comm().rank()==0)
                std::cout << "WRITE PVD FILE " << pvdname << std::endl;
            assert(pvd.is_open());
            pvd << std::fixed;
            pvd << "<?xml version=\"1.0\"?>\n"
                << "<VTKFile type=\"Collection\" version=\"0.1\">\n"
                << "<Collection>\n";
            FilenameHelper fnloop(basename,offset);
            for (unsigned int i=0; i<timesteps.size(); i++)
                {
                    //std::cout << pres << " " << path << " " << Dune::relativePath(pres,path) << " " <<fnloop.getName() << std::endl;
                    std::string fname = this->getParallelHeaderName(fnloop.getName(),Dune::relativePath(pres,path), gv.comm().size());
                    pvd << "  <DataSet timestep=\"" << timesteps[i]
                        << "\" file=\"" << fname << "\"/>\n";
                    fnloop.increment();
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
