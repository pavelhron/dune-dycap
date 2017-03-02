// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_COMMON_COMPUTEH_HH
#define DUNE_DYCAP_COMMON_COMPUTEH_HH

#include<iostream>
#include<limits>

#include <dune/common/fvector.hh>

/** \brief Compute the grid resolution h for a given grid

    For unstructured grids it computes the supremum over all cells of
    the diameter of the smallest sphere containing the cell.

    For structured grids it will just give the biggest side length of
    a cell.

    \tparam GV          The grid view
    \tparam structured  Signals a structured grid (this info can not
                        be extracted from GV).
 */
template < typename GV , bool structured>
double computeH( const GV & gv, const bool min = true )
{
  const int dim = GV::dimension;
  typedef typename GV::ctype DF;
  DF h = min ? std::numeric_limits<DF>::max() : 0.0;

  typedef typename GV::template Codim<0>::Geometry::GlobalCoordinate Coordinate;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  ElementIterator it  = gv.template begin<0>();

  if(structured){

    if(!it->geometry().type().isCube()){
      std::cout << std::endl << std::endl;
      std::cout << "#################################################" << std::endl;
      std::cout << "### WARNING: ASSUMING A STRUCTURED CUBE GRID! ###" << std::endl;
      std::cout << std::endl << std::endl;
    }

    ElementIterator eit = gv.template end<0>();
    for (; it != eit; ++it) {
      const Coordinate center_s = it->geometry().center();
      typedef typename GV::IntersectionIterator IntersectionIterator;
      IntersectionIterator iit = gv.ibegin(*it);
      IntersectionIterator eiit = gv.iend(*it);
      for (; iit != eiit; ++iit){
        if(iit->boundary())
          continue;

        DF cell_h = 0;
        const Coordinate center_n = iit->outside().geometry().center();
        for(int d=0; d<dim; ++d)
          cell_h = std::max(std::abs(center_n[d]-center_s[d]), cell_h);
        h = min ? std::min(cell_h,h) : std::max(cell_h,h);
      }
    }
  }
  else{

    ElementIterator eit = gv.template end<0>();
    for (; it != eit; ++it) {
      for(int i = 0; i < it->geometry().corners(); ++i) {
        for(int j = 0; j < i; ++j) {
          const Coordinate dist
            = it->geometry().corner(i) - it->geometry().corner(j);

          h = std::max( dist.two_norm(), h );
        } // j
      } // i
    } // eit

  }
  if (gv.comm().size()>1)
    h = gv.comm().min(h);
  return h;
}

/** \brief Compute the grid domain size

    \tparam GV          The grid view
*/
template < typename GV >
void computeGridViewBoundingBox
( const GV & gv ,
  Dune::FieldVector<typename GV::ctype, GV::dimensionworld> & min,
  Dune::FieldVector<typename GV::ctype, GV::dimensionworld> & max
  )
{
  typedef typename GV::ctype DF;
  typedef typename GV::template Codim<0>::Geometry::GlobalCoordinate Coordinate;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  ElementIterator it  = gv.template begin<0>();

  min = std::numeric_limits<DF>::max();
  max = -std::numeric_limits<DF>::max();

  ElementIterator eit = gv.template end<0>();
  for (; it != eit; ++it) {
    for(int i = 0; i < it->geometry().corners(); ++i) {
      for(int j = 0; j < i; ++j) {
        const Coordinate corner = it->geometry().corner(i);
        min = std::min( min, corner );
        max = std::max( max, corner );
      } // j
    } // i
  } // eit
  if (gv.comm().size()>1)
    {
      min = gv.comm().min(min);
      max = gv.comm().max(max);
    }
}

#endif
