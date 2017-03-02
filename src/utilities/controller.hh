// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_CONTROLLER_HH
#define DUNE_DYCAP_CONTROLLER_HH

#include <vector>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/typetree/typetree.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include"gridfunction_utilities.hh"

/** \brief Control time step after reaction.

    If some concentration is negative, then it returns false.
    Otherwise it returns true

    \tparam GV          The grid view
    \tparam V           Vector backend (containing DOF)
*/
template<class GV, class V>
bool controlReactionTimeStep (GV& gv, V& v, int verbosity=0)
{
  int passed = 1;
  const typename V::ElementType eps = 1.e-8;
  for (auto it=v.begin();it!=v.end();++it)
    if (*it < -eps)
      {
        if (verbosity)
          std::cout << "negative value is " << *it << std::endl;
        passed = 0;
        break;
      }
    else if (*it < 0)
      *it=0;

  if (gv.comm().size()>1)
    passed =  gv.comm().min(passed);
  if (passed) return true;
  else return false;
}

/** \brief control

    If some concentration is negative, then it returns false.
    Otherwise it returns true

    \tparam DGF          Discrete Grid Function

    \param dgf1          discrete grid function to be compared to dgf2
    \param border        absolute limit for the difference
*/
template<typename DGF>
const bool controlDGFDifference(const DGF& dgf1, const DGF& dgf2, double border = 0.2, int verbosity=0)
{
  typedef typename DGF::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef MinusGridFunction<GV,DGF,DGF> MGF;
  typedef typename GV::Grid::ctype RF;
  typedef typename DGF::Traits::RangeType RangeType;

  enum { dim = GV::dimension };

  // leaf entity geometry
  //typedef typename Element::Geometry Geometry;


  MGF mgf(dgf1,dgf2);

  RangeType value;
  RF maxvalue=0.0;
  int passed = 1;

  // loop once over the grid
  for (ElementIterator it = dgf1.getGridView().template begin<0>();
       it!=dgf1.getGridView().template end<0>(); ++it)
    {
      Dune::PDELab::ElementGeometry<Element> eg(*it);

      // cell geometry
      const Dune::FieldVector<RF,dim>&
        cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

      mgf.evaluate(eg.entity(),cell_center_local,value);
      if (std::abs(value)>border)
        {
          if (verbosity)
            std::cout << "Change in twophase step is too big, processor " <<  dgf1.getGridView().comm().rank()<< std::endl;
          passed = 0;
          if (std::abs(maxvalue)>std::abs(value))
            maxvalue = value;
          break;
        }
    }
  if (verbosity)
    std::cout << "max difference is " << maxvalue <<  " rank " <<  dgf1.getGridView().comm().rank() << " passed " << passed << std::endl;

  if (dgf1.getGridView().comm().size()>1)
    passed =  dgf1.getGridView().comm().min(passed);

  if (passed) return true;
  else return false;
}



//! Class to control solution.
template<class GFS>
class VelocityController {
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  enum { dim = GV::dimension };

  // leaf entity geometry
  typedef typename Element::Geometry Geometry;

  // intersection iterator type
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef typename IntersectionIterator::Intersection Intersection;


public:
  VelocityController(const GFS &gfs_) :
    gfs(gfs_)
  {
  }

  template <typename V>
  bool control (V & v)
  {
    typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
    typedef typename V::Traits::RangeFieldType RF;
    typedef typename V::Traits::RangeType RangeType;

    LFS lfs(gfs);
    // map each cell to unique id
    Dune::PDELab::ElementMapper<GV> cell_mapper(gfs.gridView());


    // loop once over the grid
    for (ElementIterator it = gfs.gridView().template begin<0>();
         it!=gfs.gridView().template end<0>(); ++it)
      {
        Dune::PDELab::ElementGeometry<Element> eg(*it);

        // cell geometry
        const Dune::FieldVector<RF,dim>&
          cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

        Dune::FieldVector<RF, dim>
          cell_center_global = eg.geometry().global(cell_center_local);

        // compute unique id
        const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

        std::cout << ids << std::endl;
        // std::cout << "cell center local " << cell_center_local << std::endl;
        //  std::cout << "cell center global " << cell_center_global << std::endl;

        // bind local function space to element
        lfs.bind(*it);
        RF velocity = 0.;

        // go through all intersections with neighbors and boundary
        IntersectionIterator isend = gfs.gridView().iend(*it);
        for (IntersectionIterator is = gfs.gridView().ibegin(*it); is!=isend; ++is)
          {
            Dune::PDELab::IntersectionGeometry<Intersection> ig(*is,0);

            // face geometry
            const Dune::FieldVector<RF,dim-1>&
              face_local = Dune::ReferenceElements<RF,dim-1>::general(ig.geometry().type()).position(0,0);
            RF face_volume = ig.geometry().volume();

            //  std::cout << "face volume " << face_volume << std::endl;
            // std::cout << "face global " << ig.geometry().global(face_local) << std::endl;

            Dune::FieldVector<RF,dim> face_center_in_element = ig.geometryInInside().global(face_local);
            //  std::cout << "face center in element " << face_center_in_element << std::endl;

            RangeType velo;
            v.evaluate(*(ig.inside()),face_center_in_element,velo);
            RF vn = velo*ig.centerUnitOuterNormal();
            std::cout << "normal velocity " << vn << std::endl;
            velocity+=vn*face_volume;

          }
        std::cout << "outflow flux " << velocity << "\n"<<std::endl;

      }

    return true;
  }

private:
  const GFS &gfs;

};



//! Class to control solution.
template<class TP, class GFS, class SL>
class SolutionController {
  typedef typename GFS::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  typedef typename SL::Traits::RangeFieldType RF;

  enum { dim = GV::dimension };

  // leaf entity geometry
  typedef typename Element::Geometry Geometry;

  // intersection iterator type
  typedef typename GV::IntersectionIterator IntersectionIterator;

  typedef typename IntersectionIterator::Intersection Intersection;


public:
  SolutionController(const GFS &gfs_, const TP& tp_, const SL &sl_old_, const SL &sl_new_) :
    gfs(gfs_), tp(tp_), sl_old(sl_old_), sl_new(sl_new_)
  {
  }

  template <typename V>
  bool control (V & v, RF deltat)
  {
    typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
    typedef typename V::Traits::RangeType RangeType;

    LFS lfs(gfs);
    // map each cell to unique id
    Dune::PDELab::ElementMapper<GV> cell_mapper(gfs.gridView());


    // loop once over the grid
    for (ElementIterator it = gfs.gridView().template begin<0>();
         it!=gfs.gridView().template end<0>(); ++it)
      {
        Dune::PDELab::ElementGeometry<Element> eg(*it);

        // cell geometry
        const Dune::FieldVector<RF,dim>&
          cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

        Dune::FieldVector<RF, dim>
          cell_center_global = eg.geometry().global(cell_center_local);

        RF cell_volume = eg.geometry().volume();

        // compute unique id
        const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);

        std::cout << ids << std::endl;
        // std::cout << "cell center local " << cell_center_local << std::endl;
        //  std::cout << "cell center global " << cell_center_global << std::endl;

        // bind local function space to element
        lfs.bind(*it);
        RF velocity = 0.;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isend = gfs.gridView().iend(*it);
        for (IntersectionIterator is = gfs.gridView().ibegin(*it); is!=isend; ++is)
          {
            Dune::PDELab::IntersectionGeometry<Intersection> ig(*is,0);

            // face geometry
            const Dune::FieldVector<RF,dim-1>&
              face_local = Dune::ReferenceElements<RF,dim-1>::general(ig.geometry().type()).position(0,0);
            RF face_volume = ig.geometry().volume();

            //  std::cout << "face volume " << face_volume << std::endl;
            // std::cout << "face global " << ig.geometry().global(face_local) << std::endl;

            Dune::FieldVector<RF,dim> face_center_in_element = ig.geometryInInside().global(face_local);
            //  std::cout << "face center in element " << face_center_in_element << std::endl;

            RangeType velo;
            v.evaluate(*(ig.inside()),face_center_in_element,velo);
            RF vn = velo*ig.centerUnitOuterNormal();
            std::cout << "normal velocity " << vn << std::endl;
            velocity+=vn*face_volume;

          }
        typename SL::Traits::RangeType sl_oldvalue, sl_newvalue;
        sl_old.evaluate(*it,cell_center_local,sl_oldvalue);
        sl_new.evaluate(*it,cell_center_local,sl_newvalue);
        RF phi = tp.phi(*it,cell_center_local);

        std::cout << "outflow flux " << velocity <<std::endl;
        std::cout << "change " << (sl_newvalue - sl_oldvalue)/deltat*phi*cell_volume << std::endl;
        std::cout << "difference " << velocity + (sl_newvalue - sl_oldvalue)/deltat*phi*cell_volume << "\n"<< std::endl;

      }

    return true;
  }

private:
  const GFS& gfs;
  const TP& tp;
  const SL& sl_old;
  const SL& sl_new;

};


//! Class to control solution.
//! it works only for FV!
template <typename DGF>
typename DGF::Traits::RangeType computeMass (DGF & dgf, typename DGF::Traits::RangeType weight=1.0)
{
  typedef typename DGF::Traits::RangeType RF;
  typedef typename DGF::Traits::DomainFieldType DF;

  typedef typename DGF::Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  enum { dim = GV::dimension };

  RF mass(0.);
  RF value;
  // loop once over the grid
  for (ElementIterator it = dgf.getGridView().template begin<0, Dune::Interior_Partition>();
       it!=dgf.getGridView().template end<0, Dune::Interior_Partition>(); ++it)
    {
      Dune::PDELab::ElementGeometry<Element> eg(*it);

      static_assert(dim==2,"dimension is not 2");

      // cell geometry
      const Dune::FieldVector<DF,dim>&
        cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

      RF cell_volume = eg.geometry().volume();
      dgf.evaluate(eg.entity(),cell_center_local,value);
      mass+=value*cell_volume;

    }

  if (dgf.getGridView().comm().size()>1)
    mass =  dgf.getGridView().comm().sum(mass);

  return mass*weight;
}




template <typename GFS, typename DGF, typename V, typename RF>
void replaceSolutionUnderLimit(const GFS& gfs, const DGF& dgf, V& v, const RF minvalue, const RF setvalue=0, int verbosity = 0)
{

  typedef typename GFS::Traits::GridViewType GV;
  const int dim = GV::dimension;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  // iterate over all cells
  typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;

  typedef typename DGF::Traits::RangeType RangeType;

  // map each cell to unique id
  Dune::PDELab::ElementMapper<GV> cell_mapper(gfs.gridView());
  ElementIterator it = gfs.gridView().template begin<0,Dune::Interior_Partition>();
  ElementIterator endit = gfs.gridView().template end<0,Dune::Interior_Partition>();
  for (; it != endit; ++it)
    {
      typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
      // cell center local
      Dune::PDELab::ElementGeometry<Element> eg(*it);

      // cell geometry
      const Dune::FieldVector<RF,dim>&
        cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

      RangeType value;
      dgf.evaluate(eg.entity(),cell_center_local,value);

      if (value<minvalue)
        {
          if (verbosity)
            std::cout << "value of " << value << " on element " << n << " was changed " << std::endl;
          v.block(n)=setvalue;
        }

    }
  // communicate limited function
  Dune::PDELab::CopyDataHandle<GFS,V> dh(gfs,v);
  if (gfs.gridView().comm().size()>1)
    gfs.gridView().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
}

template <typename GFS, typename V, typename V2, typename RF>
void copySolution(const GFS& gfs, const V& v, V2& v2, RF factor)
{

  typedef typename GFS::Traits::GridViewType GV;
  // const int dim = GV::dimension;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  // iterate over all cells
  typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;


  // map each cell to unique id
  Dune::PDELab::ElementMapper<GV> cell_mapper(gfs.gridView());
  ElementIterator it = gfs.gridView().template begin<0,Dune::Interior_Partition>();
  ElementIterator endit = gfs.gridView().template end<0,Dune::Interior_Partition>();
  for (; it != endit; ++it)
    {
      typename GV::IndexSet::IndexType n = cell_mapper.map(*it);
      v2.block(n)=v.block(n)/factor*(factor-1.);

    }
  // communicate limited function
  Dune::PDELab::CopyDataHandle<GFS,V2> dh(gfs,v2);
  if (gfs.gridView().comm().size()>1)
    gfs.gridView().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
}



#endif
