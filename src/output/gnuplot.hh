// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_GNUPLOT_HH
#define DUNE_DYCAP_GNUPLOT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <fstream>
#include <stdio.h>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <src/utilities/gridinfo.hh>

/**
 * \brief Class to execute Gnuplot instructions through pipeline
 **/
class GNUplot
{

public:
  GNUplot() throw(std::string)
  {
    gnuplotpipe=popen("gnuplot -persist","w");
    if (!gnuplotpipe) {
      throw("Gnuplot not found !");
    }
  }
  ~GNUplot()
  {
    fprintf(gnuplotpipe,"exit\n");
    pclose(gnuplotpipe);
  }

  void operator ()(const std::string& command)
  {
    fprintf(gnuplotpipe,"%s\n",command.c_str());
    fflush(gnuplotpipe);
    // flush is necessary, nothing gets plotted else
  }
  // send any command to gnuplot
protected:
  FILE *gnuplotpipe;
};


template<typename GV>
class Gnuplot1DSolution {

  typedef typename GV::Grid::ctype DF;
  // get some types
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::template Codim<0>::Geometry Geometry;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef std::size_t size_type;

  const GV &gv;
  unsigned int qorder;
  DF h;

  enum {dim = GV::dimension};

  typedef Dune::QuadratureRule<DF,dim> QR;
  typedef Dune::QuadratureRules<DF,dim> QRs;
  typedef typename QR::const_iterator QIterator;

public:
  //! Construct
  /**
   * \param fname     Name of the file to write to.
   * \param gv_       Grid View used for the vectors passed in.

   * \note The filename passed should be the same on all ranks.  Only rank
   *       0 will actually write to the file.
   */
  Gnuplot1DSolution(const GV &gv_, unsigned qorder_ = 10) :
    gv(gv_), qorder(qorder_)
  {
    h = computeH<GV,true>(gv);

    if (dim>2)
      DUNE_THROW(Dune::Exception,"Gnuplot1DSolution works only for dim<=2, your current dim is" + dim);
    if (gv.comm().rank()>0)
      DUNE_THROW(Dune::Exception,"Gnuplot1DSolution works only sequantiel!!");
  }


  //! Write a record in the output file
  /**
   * \param filename output filename
   * \param ... dgf  discrete grid function ellipses
   */
  template<typename... DGF>
  void write(const std::string& filename, const DGF&... dgf) {

    // open file
    std::ofstream file(filename.c_str());
    // write all column names
    file << "# coord\t";
    file << "\n";

    if (dim ==1)
      {
        // loop once over the grid and find elements with given x_coord
        for (ElementIterator it = gv.template begin<0>();
             it!=gv.template end<0>(); ++it)
          {
            const Geometry& geo = it->geometry();
            Dune::GeometryType gt = geo.type();
            const QR& rule = QRs::rule(gt,qorder);
            const QIterator qend = rule.end();

            Dune::FieldVector<DF, dim>
              positionl(0.01);
            Dune::FieldVector<DF, dim>
              positionl_global = geo.global(positionl);
            file << positionl_global << "\t";
            DGFoutput(file,*it,positionl,dgf...);
            file << "\n";
            //  EG eg(*it);
            //Dune::FieldVector<DF, dim>
            //  cell_center_global = eg.geometry().center();

            //file << cell_center_global[0]-h/2. << "\t";
            for (QIterator qit=rule.begin(); qit != qend; ++qit)
              {
                Dune::FieldVector<DF, dim>
                  position = qit->position();
                Dune::FieldVector<DF, dim>
                  position_global = geo.global(position);

                file << position_global << "\t";
                // evaluate the given grid functions at integration point
                DGFoutput(file,*it,position,dgf...);
                file << "\n";
              }

            Dune::FieldVector<DF, dim>
              positionr(0.99);
            Dune::FieldVector<DF, dim>
              positionr_global = geo.global(positionr);
            file << positionr_global << "\t";
            DGFoutput(file,*it,positionr,dgf...);
            file << "\n";
            //file << cell_center_global[0]+h/2. << "\t";
            //DGFoutput(file,eg,dgf...);
            //file << "\n";

          }// end it
      }
    /*   else
         {
         ElementIterator it = gv.template begin<0>();
         EG eg(*it);
         Dune::FieldVector<DF, dim>
         cell_center_global = eg.geometry().center();
         DF hmax = cell_center_global[dim-1]+1.e-6;
         // loop once over the grid and find elements with given x_coord
         for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
         {
         EG eg(*it);
         Dune::FieldVector<DF, dim>
         cell_center_global = eg.geometry().center();

         if (cell_center_global[dim-1]<hmax)
         {
         file << cell_center_global[0]-h/2. << "\t";
         DGFoutput(file,eg,dgf...);
         file << cell_center_global[0]+h/2. << "\t";
         DGFoutput(file,eg,dgf...);
         file << "\n";
         }
         }// end it

         }*/
  }

private:

  //! output from DGF (only rank 0), ellipses
  template<typename EG, typename POS, typename HDGF, typename... DGF>
  void DGFoutput(std::ostream & file, const EG& eg, const POS& pos, const HDGF& hdgf, const DGF&... dgf) {

    typename HDGF::Traits::RangeType value;
    hdgf.evaluate(eg,pos,value);
    file << value << "\t";
    DGFoutput(file,eg,pos,dgf...);
  }


  //! output from DGF
  template<typename EG, typename POS, typename HDGF>
  void DGFoutput(std::ostream& file, const EG& eg, const POS& pos, const HDGF& hdgf) {
    typename HDGF::Traits::RangeType value;
    hdgf.evaluate(eg,pos,value);
    file << value << "\t";
  }

};



#endif
