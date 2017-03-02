// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_SOLUTIONCUT_HH
#define DUNE_DYCAP_SOLUTIONCUT_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <fstream>
#include <dune/grid/common/gridenums.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <src/utilities/utilities.hh>


/** \brief Base Class to write a solution.
 *  output for different discrete grid function
 *  with ellipses
 */
template<typename GV>
class GnuplotSolutionBase {

public:
  typedef typename GV::Grid::ctype DF;
  enum {dim = GV::dimension};


  GnuplotSolutionBase(const GV &gv_)
  {
    h = computeH<GV,true>(gv_);
  }

protected:
  //! output from DGF (only rank 0), ellipses
  template<typename EG, typename HDGF, typename... DGF>
  void DGFoutput(std::ostream & file, const EG& eg, const HDGF& hdgf, const DGF&... dgf) {

    typename HDGF::Traits::RangeType value(0);
    // cell geometry
    const Dune::FieldVector<DF,dim>&
      cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);

    file << value << "\t";
    DGFoutput(file,eg,dgf...);
  }

  //! output from DGF
  template<typename EG, typename HDGF>
  void DGFoutput(std::ostream & file, const EG& eg, const HDGF& hdgf) {
    typename HDGF::Traits::RangeType value(0);
    // cell geometry
    const Dune::FieldVector<DF,dim>&
      cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);
    file << value << "\t";
  }

  DF h; // the smallest element edge size
};

/** \brief Class to write a 1D solution to the gnuplot file.
 *
 */

template<typename GV>
class GnuplotSolution: public GnuplotSolutionBase<GV> {

  typedef  GnuplotSolutionBase<GV> Base;
  typedef typename GV::Grid::ctype DF;
  // get some types
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef std::size_t size_type;

  const GV &gv;
  enum {dim = GV::dimension};

public:

  GnuplotSolution(const GV &gv_) :
    Base(gv_), gv(gv_)
  {
    if (gv.comm().rank()>0)
      DUNE_THROW(Dune::Exception,"GnuplotSolution class works only sequentiell!!");
  }

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
            EG eg(*it);

            // cell geometry
            const Dune::FieldVector<DF,dim>&
              cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

            Dune::FieldVector<DF, dim>
              cell_center_global = eg.geometry().global(cell_center_local);

            file << cell_center_global[0]-this->h/2. << "\t";
            this->DGFoutput(file,eg,dgf...);
            file << "\n";
            file << cell_center_global[0]+this->h/2. << "\t";
            this->DGFoutput(file,eg,dgf...);
            file << "\n";

          }// end it
      }
    else // in 2D only the first row of elements
      {
        ElementIterator it = gv.template begin<0>();
        EG eg(*it);

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

        Dune::FieldVector<DF, dim>
          cell_center_global = eg.geometry().global(cell_center_local);
        DF hmax = cell_center_global[dim-1]+1.e-6;
        // loop once over the grid and find elements with given x_coord
        for (ElementIterator it = gv.template begin<0>();
             it!=gv.template end<0>(); ++it)
          {
            EG eg(*it);

            // cell geometry
            const Dune::FieldVector<DF,dim>&
              cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

            Dune::FieldVector<DF, dim>
              cell_center_global = eg.geometry().global(cell_center_local);

            if (cell_center_global[dim-1]<hmax)
              {
                file << cell_center_global[0]-this->h/2. << "\t";
                this->DGFoutput(file,eg,dgf...);
                file << cell_center_global[0]+this->h/2. << "\t";
                this->DGFoutput(file,eg,dgf...);
                file << "\n";
              }
          }// end it
      }
  }
};

/** \brief Class to write a solution CUT to a gnuplot file.
 *
 * Writes a cut through the solution (1D cut in 2D grid or 2D cut in 3D grid) to a file.
 * It works for arbitrary number of DiscreteGridFunctions, sequentiel and parallel
 *
 * \tparam GV GridView
 */

template<typename GV>
class GnuplotCutSolutionInTime: public GnuplotSolutionBase<GV> {

  typedef  GnuplotSolutionBase<GV> Base;
  typedef typename GV::Grid::ctype DF;
  // get some types
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef std::size_t size_type;
  typedef typename GV::Grid::LeafIndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  const GV &gv;
  using Base::h;
  std::vector<ElementIterator> egstorage;
  std::string fname;
  const std::string &comment;
  bool horizontal;
  std::ofstream s;

  enum {dim = GV::dimension};

public:
  //! Construct
  /**
   * \param fname     Name of the file to write to.
   * \param gv_       Grid View used for the vectors passed in.

   * \note The filename passed should be the same on all ranks.  Only rank
   *       0 will actually write to the file.
   */
  GnuplotCutSolutionInTime( const GV &gv_, const DF x_coord, const std::string &fname_, const std::string &comment_, bool horizontal_=false) : Base(gv_),
                                                                                                                                               gv(gv_), fname(fname_), comment(comment_), horizontal(horizontal_)
  {
    // loop once over the grid and find elements with given x_coord
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {
        // skip ghost and overlap
        if (it->partitionType()!=Dune::InteriorEntity)
          continue;

        EG eg(*it);
        const size_type nr_corners = eg.geometry().corners();
        size_type nr_corners_half = nr_corners/2;

        size_type nrx = 0;

        if (!horizontal)
          for (unsigned int i=0; i<nr_corners; i++)
          {
            Dune::FieldVector<DF,dim> corner = eg.geometry().corner(i);

            if (std::abs(corner[0]-x_coord)<1.e-6)
              {
                Dune::FieldVector<DF, dim> center = eg.geometry().center();
                DUNE_THROW(Dune::Exception,"Element wit coordinates " << center << " is on the x-coordinate " << x_coord << " for cut output, corner " << corner);

              }
            if (corner[0] < x_coord)
              ++nrx;
          }
        else if (horizontal && dim>1)
          for (unsigned int i=0; i<nr_corners; i++)
          {
            Dune::FieldVector<DF,dim> corner = eg.geometry().corner(i);

            if (std::abs(corner[1]-x_coord)<1.e-6)
              {
                Dune::FieldVector<DF, dim> center = eg.geometry().center();
                DUNE_THROW(Dune::Exception,"Element wit coordinates " << center << " is on the x-coordinate " << x_coord << " for cut output, corner " << corner);

              }
            if (corner[1] < x_coord)
              ++nrx;
          }
        // stores appropriate element iterators
        if (nrx == nr_corners_half)
          egstorage.push_back(it);

      }// end it
  }


  //help structure
  struct mergedList
  {
    double coordvalue;
    std::string text;
    friend std::ostream& operator <<(std::ostream& os, mergedList A)
    {
      os << A.coordvalue << "\t" << A.text;
      return os;
    }

    friend std::istream& operator >>(std::istream& is, mergedList& A)
    {
      is >> A.coordvalue;
      std::getline(is, A.text);
      return is;
    }

  };


  struct my_sorter {
    bool operator() (mergedList one, mergedList two) { return one.coordvalue < two.coordvalue ; }
  };

  // function reads a line from a given file and skipes all lines starting with an '#' character
  void skipComment(std::ifstream &infile)
  {
    if (infile.eof())
      return;
    do
      {
        char c=0;
        if (infile.get(c))
          {
            if (c=='#')
              {
                std::string dummy;
                std::getline(infile,dummy);
              }
            else if ((c!='\a')&&(c!='\b')&&(c!='\f')&&(c!='\n')&&(c!='\r')&&(c!='\t')&&(c!='\v')&&(c!=' ')&&(c!='\'')&&(c!='\"'))
              {
                infile.unget(); // decrement get pointer
                break;
              }
          }
      }
    while (!infile.eof());
  }


  //! Write a record in the output file (Element index)
  /**
   * \param t        time
   * \param x_coord  x-coordinate
   * \param ... dgf  discrete grid function ellipses
   */
  template<class Time, typename... DGF>
  void output(Time t, const DGF&... dgf) {

    std::string sname =  fname + "t" + std::to_string(static_cast<long long int>(t));
    std::string psname =   sname + "_p" + std::to_string(static_cast<long long int>(gv.comm().rank())) + ".dat";

    //std::cout << "opening file " << psname.c_str() << std::endl;
    // open file for each processor
    s.exceptions(std::ios_base::badbit | std::ios_base::eofbit
                 | std::ios_base::failbit);
    s.open(psname.c_str());
    if (!s.is_open())
      DUNE_THROW(Dune::IOError, "Could not open file " << psname.c_str()  << "!");
    //else
    //  std::cout << "File " << psname.c_str() << " was opened" << std::endl;
    s << "#test";
    s << std::setprecision(14) << std::scientific;

    // loop over all elements from egstorage
    for(auto pit=egstorage.begin();pit!=egstorage.end();++pit)
      {
        EG eg(*(*pit));
        Dune::FieldVector<DF, dim>
          cell_center_global = eg.geometry().center();
        if (dim>2)
          s << cell_center_global[1] << "\t";
        if (!horizontal)
          s << cell_center_global[dim-1] << "\t";
        else
          s << cell_center_global[0] << "\t";
        this->DGFoutput(s,eg,dgf...);
        s << "\n";
      }
    s.close();
    gv.comm().barrier();

    if (gv.comm().rank()==0)
      {

        std::ifstream indata;
        std::ofstream outdata;

        std::string outputname = sname + ".dat";
        outdata.open(outputname.c_str());

        // This can be a vector. No need for array here.
        std::vector<mergedList> D1;
        for (auto i=0;i<gv.comm().size();i++)
          {

            std::string psname =   sname + "_p" + std::to_string(static_cast<long long int>(i))+ ".dat";
            //std::cout << "processing file " << i << " with name " << psname.c_str()<< std::endl;
            indata.open(psname.c_str());
            skipComment(indata);

            mergedList tmp;
            while (indata >> tmp)
              D1.push_back(tmp);
            indata.close();
            //remove(psname.c_str());
            //  DUNE_THROW(Dune::Exception,"File " << psname << " was not deleted.");
          }


        //std::cout << "Before sorting" << std::endl;
        //std::copy(D1.begin(), D1.end(), std::ostream_iterator<mergedList>(std::cout, "\n"));


        // Sort the vector using the std::sort algorithm.
        // http://www.cplusplus.com/reference/algorithm/sort/ for an example
        std::sort(D1.begin(), D1.end(), my_sorter());

        //std::cout << "After sorting" << std::endl;
        // std::copy(D1.begin(), D1.end(), std::ostream_iterator<mergedList>(std::cout, "\n"));

        // Write the sorted list to the output file
        std::copy(D1.begin(), D1.end(), std::ostream_iterator<mergedList>(outdata, "\n"));
        outdata.close();
      }

  }

};





template<typename GV>
class GnuplotCutSolution :public GnuplotSolutionBase<GV> {

  typedef  GnuplotSolutionBase<GV> Base;
  typedef typename GV::Grid::ctype DF;
  // get some types
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef std::size_t size_type;
  typedef typename GV::Grid::LeafIndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  const GV &gv;
  std::ofstream s;
  bool enabled;

  enum {dim = GV::dimension};

public:
  //! Construct
  /**
   * \param fname     Name of the file to write to.
   * \param gv_       Grid View used for the vectors passed in.

   * \note The filename passed should be the same on all ranks.  Only rank
   *       0 will actually write to the file.
   */
  GnuplotCutSolution(const GV &gv_, const std::string &fname, const std::string &comment="")  DUNE_DEPRECATED : Base(gv_),
    gv(gv_)
  {
    if (dim>2)
      DUNE_THROW(Dune::Exception,"GnuplotCutSolution works only for dim<=2, your current dim is" + dim);
    if (gv.comm().rank()>0)
      DUNE_THROW(Dune::Exception,"GnuplotCutSolution works only sequantiel!!");

    enabled = fname.size() > 0;
    if(enabled) {

      std::ostringstream sname;
      sname  << fname;
      if(gv.comm().rank() == 0) {
        s.exceptions(std::ios_base::badbit | std::ios_base::eofbit
                     | std::ios_base::failbit);
        s.open(sname.str());
        s << std::setprecision(14) << std::scientific;
        if(comment.size() > 0)
          s << "# " << comment << "\n";
        s << "#results, cut in y direction of the domain\n";
      }
    }

  }

  //! Write a record in the output file (Element index)
  /**
   * \param t        time
   * \param x_coord  x-coordinate
   * \param ... dgf  discrete grid function ellipses
   */
  template<class Time, typename... DGF>
  void output(Time t, DF x_coord, const DGF&... dgf) {

    s << "\n";
    s << "#time " << t << " cut in y direction, x= " << x_coord << "\n";
    std::vector<ElementIterator> egstorage;

    // loop once over the grid and find elements with given x_coord
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {
        EG eg(*it);

        size_type nr_corners = eg.geometry().corners();
        size_type nr_corners_half = nr_corners/2;

        size_type nrx = 0;
        size_type nry = 0;

        for (unsigned int i=0; i<nr_corners; i++)
          {
            Dune::FieldVector<DF,dim> corner = eg.geometry().corner(i);
            if (corner[0] < x_coord)
              ++nrx;
            if (corner[0] > x_coord)
              ++nry;
          }

        // stores appropriate element iterators
        if (nrx == nry  && nry== nr_corners_half)
          egstorage.push_back(it);
      }// end it

    //  const IndexSet& indexset = gv.indexSet();

    // loop over all elements from egstorage
    for(auto pit=egstorage.begin();pit!=egstorage.end();++pit)
      {
        EG eg(*(*pit));
        Dune::FieldVector<DF, dim>
          cell_center_global = eg.geometry().center();


        s << cell_center_global[dim-1] << "\t";
        this->DGFoutput(s,eg,dgf...);
        s << "\n";
      }
  }


};




#endif
