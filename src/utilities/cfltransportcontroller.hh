// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_CFLTRANSPORTCONTROLLER_HH
#define DUNE_DYCAP_CFLTRANSPORTCONTROLLER_HH

#include <dune/common/shared_ptr.hh>
#include <dune/common/typetraits.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include<src/models/fluxreconstruction.hh>

#include<dune/pdelab/gridfunctionspace/localfunctionspace.hh>

using namespace Dune::PDELab;

/**
 * \todo{test how fast is it - we have to go through the whole domain
 *       test it in parallel
 *       delete text output if necessary
 *       is the capacity function correct?
 *       collect the min time step doesn't work in parallel}
 */


/**
 * \brief class to control the timestep size used for the component transport

 * \tparam GV        grid function space
 * \tparam TP         twophase parameters class
 * \tparam FR         flux reconstruction class (for FV)
 */
template<typename GV, typename TP, typename FR=DefaultFluxReconstruction<typename TP::Traits> >
class CFLTransportController
{
  //! \brief export type for range field
  typedef typename TP::Traits::RangeFieldType RF;

  //! \brief export type for doman field
  typedef typename TP::Traits::DomainFieldType DF;

  typedef typename TP::BCType BCType;

  //! \brief elementiterator
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;

  //! \brief export type of entity
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  //! \brief intersectioniterator
  typedef typename GV::IntersectionIterator IntersectionIterator;

  //! \brief export type of intersection
  typedef typename IntersectionIterator::Intersection Intersection;

  //! \brief Enum for domain dimension
  enum {dim = GV::dimension};

  //! \brief flux reconstruction class
  typedef FR FluxReconstruction;

public:

  //! A constructor with default flux
  CFLTransportController(const GV & gv_, TP & tp_, RF z_=1.e-7, int verbosity_ = 0)
    : gv(gv_), tp(tp_), flux(std::make_shared<FR>(fluxDefault())), zero(z_), verbosity(verbosity_)
  {
    if(gv.comm().rank()>0)
      verbosity=0;
  }

  //! A constructor with flux f_
  CFLTransportController(const GV & gv_, TP & tp_, FluxReconstruction & f_, RF z_=1.e-7, int verbosity_ = 0)
    : gv(gv_), tp(tp_), flux(std::make_shared<FR>(f_)), zero(z_), verbosity(verbosity_)
  {
    if(gv.comm().rank()>0)
      verbosity=0;
  }

  //! A constructor with flux f_
  CFLTransportController(const GV & gv_, TP & tp_, std::shared_ptr<FluxReconstruction> f_, RF z_=1.e-7, int verbosity_ = 0)
    : gv(gv_), tp(tp_), flux(f_), zero(z_), verbosity(verbosity_)
  {
    if(gv.comm().rank()>0)
      verbosity=0;
  }


  RF computeNewTimeStep()
  {
    // set min timestep
    RF dtmin = 1.e6;

    // unique id mapper
    ElementMapper<GV> cell_mapper(gv);

    // Traverse grid view
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {

        // element geometry for element with *it
        ElementGeometry<Element> eg(*it);

#ifdef CFLCONTROL_DEBUG
        const typename GV::IndexSet::IndexType ids = cell_mapper.map(*it);
#endif
        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // initialize cell outflux
        RF celloutflux = 0.0;

        // get saturation
        RF snew = tp.snew(eg.entity(),inside_local);

        // find out if the cell is active
        bool active_cell =  snew>zero;

        // if cell is not active, skip, no fluxes
        if (!active_cell)
          {
#ifdef CFLCONTROL_DEBUG
            // get unique id
            std::cout << "cell " << ids << " is not active, new saturation " << tp.snew(eg.entity(),inside_local) << " old saturation" << tp.sold(eg.entity(),inside_local) << std::endl;
#endif
            continue;
          }

        // Traverse intersections
        unsigned int intersection_index = 0;
        IntersectionIterator endit = gv.iend(*it);
        IntersectionIterator iit = gv.ibegin(*it);
        for(; iit!=endit; ++iit, ++intersection_index)
          {
            typedef IntersectionGeometry<Intersection> IG;
            IG ig(*iit,intersection_index);

            // face geometry
            const Dune::FieldVector<DF,IG::dimension-1>&
              face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            RF face_volume = ig.geometry().volume();
            // if (dim==1)
            // face_volume*=ig.inside()->geometry().volume();
            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside().type()).position(0,0);

            // face center in element coordinates, corresponds to ig.geometryInInside().center()
            Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

            // evaluate velocity
            typename TP::Traits::RangeType v(tp.v(ig.inside(),face_center_in_element));

            // the normal velocity
            RF vn = v*ig.centerUnitOuterNormal();

            // outflux from this cell
            if (vn>=0)
              celloutflux += vn*face_volume/flux->timestepFactor();

            // evaluate diffusion coefficients
            typename TP::Traits::RangeFieldType D_inside = tp.D(ig.inside(),inside_local);

            Dune::FieldVector<DF,IG::dimension>
              inside_global = ig.inside().geometry().center();

            if (iit->neighbor())
              {

                const Dune::FieldVector<DF,IG::dimension>&
                  outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside().type()).position(0,0);

                typename TP::Traits::RangeFieldType D_outside = tp.D(ig.outside(),outside_local);
                typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));


                // distance between cell centers in global coordinates
                Dune::FieldVector<DF,IG::dimension>
                  outside_global = ig.outside().geometry().center();
                inside_global -= outside_global;
                RF distance = inside_global.two_norm();


                // add diffusion to the outflux
                celloutflux += D_avg*face_volume/distance;
              }

            else
              {
                // evaluate boundary condition type
                BCType bc = tp.bctype(ig.intersection(),face_local);
                if (bc==BCType::Dirichlet) // Dirichlet boundary
                  {
                    Dune::FieldVector<DF,IG::dimension>
                      outside_global = ig.geometry().center();
                    inside_global -= outside_global;
                    RF distance = inside_global.two_norm();
                    celloutflux += D_inside*face_volume/distance;
                    // std::cout << "face volume is " << face_volume << " volume " << ig.inside()->geometry().volume() <<std::endl;
                  }
              }

          } // iit

        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = tp.c(eg.entity(),inside_local)*eg.geometry().volume();
        //if (dim==1)
        //  cellcapacity*=eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(celloutflux+1E-40);

        if (celldt < 0)
           DUNE_THROW(Dune::Exception, "cell capacity " << cellcapacity << " celloutflux " << celloutflux);


#ifdef CFLCONTROL_DEBUG
        if (dtmin < celldt)
        std::cout << "cell " << ids
                  << "\ncell capacity: " << cellcapacity
                  << "\ncell outflux:  " << celloutflux  << " cell dt " << celldt << std::endl;
#endif

        if (dtmin < celldt)
          {
            //std::cout << "cell capacity " << cellcapacity << " celloutflux " << celloutflux << std::endl;
          }
        // get the minimum
        dtmin = std::min(dtmin,celldt);

      } // it
    return dtmin;
  }


  //! find the maximal suggested time step
  RF suggestTimeStep()
  {

    RF suggested_dt = computeNewTimeStep();

    //gv.comm().barrier();
    // std::cout << "minimal timestep is " << dtmin << " with fluxfactor " <<  flux.timestepFactor() << " rank " << gv.comm().rank()<< std::endl;
    // gv.comm().barrier();
    // min time step has to be multiplied by reconstruction factor
    if (gv.comm().size()>1)
          suggested_dt =  gv.comm().min(suggested_dt);


    if (verbosity)
      std::cout << "cfl controller timestep is " << suggested_dt*flux->timestepFactor() << " with fluxfactor " << flux->timestepFactor() << std::endl;
    return suggested_dt;
  }

private:
  static FluxReconstruction & fluxDefault()
  {
    static FluxReconstruction f;
    return f;
  }


  GV gv;
  TP & tp;
  std::shared_ptr<FluxReconstruction> flux;
  RF zero;
  int verbosity;
};



#endif
