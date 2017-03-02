// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_COMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_COMPONENTTRANSPORTOP_HH

#include<src/utilities/rt0qfem.hh>
#include<dune/localfunctions/lagrange/q1.hh>
#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include"transporttraits.hh"
#include "fluxreconstruction.hh"

/*
  need to test if it works properly!!

 */

#define COMPTRANST_JACOBIAN false

namespace Dune {
  namespace PDELab {


    /** a local operator for a cell-centered finite folume scheme for
        the transport equation

        \nabla \cdot \{v u - D \nabla u \} = q in \Omega
        u = g on \Gamma_D
        \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
        outflow on \Gamma_O

        Modified version for the case

        d_t (c(x,t)u(x,t)) + \nabla \cdot \{v u - D \nabla u \} = q in \Omega

        where c(x,t) may become zero. We assume that the following holds:

        c(x,t+dt) <= eps  ==>  c(x,t) <= eps
        m
        \tparam TP  parameter class implementing ComponentTransportParameterInterface
    */
    template<typename TP, typename FR=DefaultFluxReconstruction<typename TP::Traits>, typename RA = ReactionBaseAdapter>
    class ModifiedCCFVSpatialTransportOperator :
      public NumericalJacobianSkeleton<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianBoundary<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianVolume<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianApplySkeleton<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianApplyBoundary<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianApplyVolume<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianApplyVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public NumericalJacobianVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP,FR,RA> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      typedef FR FluxReconstruction;
      typedef typename TP::BCType BCType;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaVolumePostSkeleton = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume    = true };

      enum { doSkeletonTwoSided = true }; // need to see face from both sides for CFL calculation

      ModifiedCCFVSpatialTransportOperator (TP& tp_, typename TP::Traits::RangeFieldType z=1e-7, bool neighbor_computation_=false)
        : tp(tp_), flux(fluxDefault()), ra(raDefault()), zero(z), neighbor_computation(neighbor_computation_)
      {
      }

      ModifiedCCFVSpatialTransportOperator (TP& tp_, FluxReconstruction & f, typename TP::Traits::RangeFieldType z=1e-7, bool neighbor_computation_=false)
        : tp(tp_), flux(f), ra(raDefault()), zero(z), neighbor_computation(neighbor_computation)
      {
      }

      ModifiedCCFVSpatialTransportOperator (TP& tp_, RA& ra_, typename TP::Traits::RangeFieldType z=1e-7, bool neighbor_computation_=false)
        : tp(tp_), flux(fluxDefault()), ra(ra_), zero(z), neighbor_computation(neighbor_computation_)
      {
      }

      ModifiedCCFVSpatialTransportOperator (TP& tp_, RA& ra_, FluxReconstruction & f, typename TP::Traits::RangeFieldType z=1e-7, bool neighbor_computation_=false)
        : tp(tp_), flux(f), ra(ra_), zero(z), neighbor_computation(neighbor_computation)
      {
      }


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        snew = tp.snew(eg.entity(),inside_local);
        active_cell = snew>zero;
        if (active_cell)
          active_cell_count++;
        else
          inactive_cell_count++;

        celloutflux = 0.0; // prepare dt computation

        if (tp.computeReaction())
          ra.evaluate(eg.entity(),lfsu,x,lfsv,r);

        // evaluate source term
        typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local);
        if (q<0)
          {
            r.accumulate(lfsv,0,-q*eg.geometry().volume()*x(lfsu,0));
            //std::cout << "alpha value is " << q << std::endl;
          }
#ifdef COMPTRANS_DEBUG
        std::cout << "alpha_volume: time=" << time
                  << " pos=" << eg.geometry().center()
                  << " snew=" << snew
                  << " inactive_cell_count=" << inactive_cell_count
                  << " active_cell_count=" << active_cell_count
                  << " active_cell=" << active_cell
                  << " x=" << x(lfsu,0)
                  << std::endl;
#endif // COMPTRANS_DEBUG

      }
      /*
#if COMPTRANST_JACOBIAN
      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        snew = tp.snew(eg.entity(),inside_local);
        active_cell = snew>zero;

      }
#endif // COMPTRANST_JACOBIAN
      */
      // skeleton integral depending on test and ansatz functions
      // each face is only visited TWICE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside().type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside().type()).position(0,0);

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(ig.inside(),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        // activity of neighbor

        // skip if the neighbor is not active
        if (neighbor_computation)
          {
            bool active_neighbor = tp.snew(ig.outside(),outside_local)>zero;

            if (!active_neighbor)
              return;
          }

        // convective flux
        RF u_upwind=0.0;
        if (vn>=0)
          {
            // evaluate face flux on inside
            flux.evaluate(ig.inside(), ig.geometryInInside().center(), x_s(lfsu_s,0), u_upwind);
          }
        else
          {
            // evaluate face flux on outside
            flux.evaluate(ig.outside(), ig.geometryInOutside().center(), x_n(lfsu_n,0), u_upwind);
          }

        r_s.accumulate(lfsu_s,0,(u_upwind*vn)*face_volume);

        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

        // evaluate diffusion coefficients
        typename TP::Traits::RangeFieldType D_inside = tp.D(ig.inside(),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.D(ig.outside(),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside().geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside().geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // diffusive flux
        // note: we do only one-sided evaluation here
        r_s.accumulate(lfsu_s,0,-(D_avg*(x_n(lfsu_n,0)-x_s(lfsu_s,0))/distance)*face_volume);

        celloutflux += D_avg*face_volume/distance;
      }


#if COMPTRANST_JACOBIAN
      // jacobian of skeleton term
      // each face is only visited TWICE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();


        // activity of neighbor
        bool active_neighbor = tp.snew(*(ig.outside()),outside_local)>zero;

        // skip if the neighbor is not active
        if (neighbor_computation)
          {

            if(active_neighbor)
              return;
          }

        // convective flux
        if (vn>=0)
          {
            mat_ss.accumulate(lfsu_s,0,lfsu_s,0,vn*face_volume);
          }
        else
          {
            if(active_neighbor)
              mat_sn.accumulate(lfsu_n,0,lfsu_n,0,vn*face_volume);
          }

        // evaluate diffusion coefficients
        typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // diffusive flux
        // note: we do only one-sided evaluation here
        mat_ss.accumulate(lfsu_s,0,lfsu_s,0,D_avg/distance*face_volume);

        if(active_neighbor)
          mat_sn.accumulate(lfsu_n,0,lfsu_n,0,-D_avg/distance*face_volume);

      }
#endif // COMPTRANST_JACOBIAN


      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        BCType bctype = tp.bctype(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bctype==BCType::Neumann) // Neumann boundary
          {
            typename TP::Traits::RangeFieldType j = tp.j(ig.intersection(),face_local);
            r_s.accumulate(lfsu_s,0,j*face_volume);
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(ig.inside(),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        vmax = std::max(vmax, std::abs(vn));

        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

        if (bctype==BCType::Outflow) // Outflow boundary
          {
            /*  typename TP::Traits::RangeFieldType g;
            g=tp.g(ig.intersection(),face_local);

            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension>
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension>
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            r_s.accumulate(lfsu_s,0,(vn*x_s(lfsu_s,0) - D_inside*(g-x_s(lfsu_s,0))/distance)*face_volume);
            return;
            */
            typename TP::Traits::RangeFieldType g;
            g=tp.g(ig.intersection(),face_local);
            r_s.accumulate(lfsu_s,0,vn*x_s(lfsu_s,0)*face_volume);
            return;
          }


        if (bctype==BCType::Dirichlet) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            g=tp.g(ig.intersection(),face_local);

            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside().type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.D(ig.inside(),inside_local);
            Dune::FieldVector<DF,IG::dimension>
              inside_global = ig.inside().geometry().center();
            Dune::FieldVector<DF,IG::dimension>
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            r_s.accumulate(lfsu_s,0,(g*vn - D_inside*(g-x_s(lfsu_s,0))/distance)*face_volume);
            // r_s.accumulate(lfsu_s,0,( - D_inside*(g-x_s(lfsu_s,0))/(2.*distance))*face_volume);
            // std::cout << "hlocal is " << distance << std::endl;
            celloutflux += D_inside*face_volume/distance;
            return;
          }

      }

      /*
#if COMPTRANST_JACOBIAN
      // jacobian of boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // no fluxes if cell is not active !
        if (!active_cell) return;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        BCType bctype = tp.bctype(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bctype==BCType::Neumann) // Neumann boundary
          {
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        if (bctype==BCType::Outflow) // Outflow boundary
          {
            mat_ss.accumulate(lfsu_s,0,lfsu_s,0,vn*face_volume);
            return;
          }

        if (bctype==BCType::Dirichlet) // Dirichlet boundary
          {
            if (vn>=0)
              mat_ss.accumulate(lfsu_s,0,lfsu_s,0,vn*face_volume);
            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension>
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension>
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            mat_ss.accumulate(lfsu_s,0,lfsu_s,0,D_inside/distance*face_volume);
            return;
          }
      }
#endif // COMPTRANST_JACOBIAN
      */
      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        // typedef typename LFSU::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        if (!first_stage) return; // time step calculation is only done in first stage

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);


        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = tp.c(eg.entity(),inside_local)*eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(celloutflux+1E-40);


#ifdef COMPTRANS_DEBUG
        if (active_cell)
          std::cout << "A: time=" << time
                    << " pos=" << eg.geometry().center()
                    << " snew=" << snew
                    << " capacityvol=" << cellcapacity
                    << " outflux=" << celloutflux
            //    << " residual=" << r[0]
                    << " celldt=" << celldt
                    << " dtmin " <<  std::min(dtmin,celldt)
                    << " x=" << x(lfsu,0)
                    << std::endl;
#endif // COMPTRANS_DEBUG

        dtmin = std::min(dtmin,celldt);
      }


      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        // typedef typename LFSV::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        // no sources if cell is not active !
        if (!active_cell) return;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate source term
        typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local);

        if (q>0)
          {
            //std::cout << "lambda value is " << q << std::endl;
            r.accumulate(lfsv,0,-q*eg.geometry().volume());

          }
          }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        flux.setTime(t);
        ra.setTime(t);
        tp.setTime(t);
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        tp.preStep(time,dt,stages);
        flux.updateReconstruction(time, stages);
      }

      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
        //         std::cout << "preStage on transport operator called" << std::endl;
        if (r==1)
          {
            first_stage = true;
            dtmin = 1E100;
            active_cell_count = 0;
            inactive_cell_count = 0;
            vmax = 0.;
          }
        else first_stage = false;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
#ifdef COMPTRANS_DEBUG
        std::cout << "active cells: " << active_cell_count << " dtmin: " << dtmin << std::endl;
        std::cout << "inactive cells: " << inactive_cell_count << " dtmin: " << dtmin << " suggest " << dtmin * flux.timestepFactor()<<  " vmax " << vmax << std::endl;
#endif // COMPTRANS_DEBUG

        //      std::cout << "inactive cells: " << inactive_cell_count << " dtmin: " << dtmin << " suggest " << dtmin * flux.timestepFactor()<< std::endl;
        return dtmin * flux.timestepFactor();
      }

    private:

      static RA & raDefault()
      {
        static RA ra;
        return ra;
      }

      typename TP::Traits::RangeFieldType regularization (typename TP::Traits::RangeFieldType x)
      {
        const typename TP::Traits::RangeFieldType regeps = 1e-20;
        const typename TP::Traits::RangeFieldType regeps2 = regeps*regeps;
        if (x<-regeps) return -1;
        if (x<0.0) return (x+regeps)*(x+regeps)/regeps2-1.0;
        if (x<regeps) return 1.0-(regeps-x)*(regeps-x)/regeps2;
        return +1.0;
      }

      static FluxReconstruction & fluxDefault()
      {
        static FluxReconstruction f;
        return f;
      }

      TP& tp;
      FluxReconstruction& flux;
      RA& ra;
      bool first_stage;
      mutable bool active_cell;
      mutable typename TP::Traits::RangeFieldType snew;
      mutable int active_cell_count;
      mutable int inactive_cell_count;
      typename TP::Traits::RangeFieldType time;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      mutable typename TP::Traits::RangeFieldType celloutflux;
      typename TP::Traits::RangeFieldType zero;
      const bool neighbor_computation;
      mutable typename TP::Traits::RangeFieldType vmax;
    };




    /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x,t) uv dx
     * \f}
     *
     * version where c(x,t) may become zero.
     */
    template<class TP>
    class ModifiedCCFVTemporalOperator
      : public NumericalJacobianVolume<ModifiedCCFVTemporalOperator<TP> >,
      public NumericalJacobianApplyVolume<ModifiedCCFVTemporalOperator<TP> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      typedef typename TP::BCType BCType;

      ModifiedCCFVTemporalOperator (TP& tp_)
        : tp(tp_), zero(1e-7)
      {
      }

      ModifiedCCFVTemporalOperator (TP& tp_, typename TP::Traits::RangeFieldType z)
        : tp(tp_), zero(z)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local);
        bool active_cell = snew>zero;

        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);

        if (active_cell)
          {
            // residual contribution
            r.accumulate(lfsu,0,c*x(lfsu,0)*eg.geometry().volume());
          }

#ifdef COMPTRANST_DEBUG
        if (active_cell)
          std::cout << "M: time=" << time
            //    << " pos=" << eg.geometry().center()
                    << " snew=" << snew
                    << " capacityvol=" << c*eg.geometry().volume()
            // << " residual=" << r[0]
                    << std::endl;
#endif // COMPTRANS_DEBUG
      }

      //#ifdef COMPTRANST_JACOBIAN
      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate saturation at end of big step for cell activity
        typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local);
        bool active_cell = snew>zero;

        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);

        // residual contribution
        if (active_cell)
          mat.accumulate(lfsu,0,lfsu,0,c*eg.geometry().volume());
        else
          mat.accumulate(lfsu,0,lfsu,0,1.0); // rhs should be zero so we get zero update

#ifdef COMPTRANST_DEBUG
        if (active_cell)
          std::cout << "J: time=" << time
                    << " pos=" << eg.geometry().center()[0] << " " << eg.geometry().center()[1]
                    << " snew=" << snew
            //<< " mat=" << mat
                    << std::endl;
#endif // COMPTRANS_DEBUG
      }
      //#endif // COMPTRANST_JACOBIAN

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
      }

      //suggest time step, asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return std::numeric_limits<typename TP::Traits::RangeFieldType>::max(); //initial value should be big enough
      }



    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
      typename TP::Traits::RangeFieldType zero;
    };



    /** a local operator for a cell-centered finite folume scheme for
        the transport equation

        \nabla \cdot \{v u - D \nabla u \} = q in \Omega
        u = g on \Gamma_D
        \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
        outflow on \Gamma_O

        Implicit version for the case

        d_t (c(x,t)u(x,t)) + \nabla \cdot \{v u - D \nabla u \} = q in \Omega

        where c(x,t) may become zero. We assume that the following holds:

        c(x,t+dt) <= eps  ==>  c(x,t) <= eps
        m
        \tparam TP  parameter class implementing ComponentTransportParameterInterface
    */
    template<typename TP, typename RA = ReactionBaseAdapter >
    class ImplicitCCFVSpatialTransportOperator :
      public NumericalJacobianSkeleton<ImplicitCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplySkeleton<ImplicitCCFVSpatialTransportOperator<TP,RA> >,

      public NumericalJacobianBoundary<ImplicitCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplyBoundary<ImplicitCCFVSpatialTransportOperator<TP,RA> >,

      public NumericalJacobianVolume<ImplicitCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplyVolume<ImplicitCCFVSpatialTransportOperator<TP,RA> >,

      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      typedef typename TP::BCType BCType;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      ImplicitCCFVSpatialTransportOperator (TP& tp_, const bool& central_ = false)
        : tp(tp_), ra(raDefault()), central(central_)
      {
      }


      ImplicitCCFVSpatialTransportOperator (TP& tp_, RA& ra_, const bool& central_ = false)
        : tp(tp_), ra(ra_), central(central_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        //typedef LFSV Space;

        // domain and range field type
        //typedef typename Space::Traits::FiniteElementType::
        // Traits::LocalBasisType::Traits::DomainFieldType DF;
        //typedef typename Space::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeFieldType RF;
        //const int dim = EG::Geometry::dimension;

        // cell center
        //  const Dune::FieldVector<DF,dim>&
        //  inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate capacity
        // typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);

        // residual contribution
        //r.accumulate(lfsu,0,x(lfsu,0)*eg.geometry().volume()*lambda*c);

        ra.evaluate(eg.entity(),lfsu,x,lfsv,r);

      }


      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename LFSU::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeType RangeType;
        // typedef typename LFSU::Traits::GridFunctionSpaceType::Traits::BackendType B;


        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        // convective flux
        RF u_upwind=0.0;
        if (vn>=0) u_upwind = x_s(lfsu_s,0); else u_upwind = x_n(lfsu_n,0);


        // evaluate diffusion coefficients
        typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
        typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local);
        typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();


        if (central==true)
          u_upwind = (x_n(lfsu_n,0)+x_s(lfsu_s,0))/2.0;

        r_s.accumulate(lfsu_s,0,(u_upwind*vn)*face_volume);
        r_n.accumulate(lfsu_n,0,-(u_upwind*vn)*face_volume);


        // diffusive flux
        // note: we do only one-sided evaluation here
        // diffusive flux
        r_s.accumulate(lfsu_s,0,-(D_avg*(x_n(lfsu_n,0)-x_s(lfsu_s,0))/distance)*face_volume);
        r_n.accumulate(lfsu_n,0,(D_avg*(x_n(lfsu_n,0)-x_s(lfsu_s,0))/distance)*face_volume);
      }



      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        // evaluate boundary condition type
        BCType bctype = tp.bctype(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bctype==BCType::Neumann) // Neumann boundary
          {
            typename TP::Traits::RangeFieldType j = tp.j(ig.intersection(),face_local);
            r_s.accumulate(lfsu_s,0,j*face_volume);
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();

        if (bctype==BCType::Outflow) // Outflow boundary
          {
            r_s.accumulate(lfsu_s,0,vn*x_s(lfsu_s,0)*face_volume);
            return;
          }

        if (bctype==BCType::Dirichlet) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            g=tp.g(ig.intersection(),face_local);

            const Dune::FieldVector<DF,IG::dimension>&
              inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local);
            Dune::FieldVector<DF,IG::dimension>
              inside_global = ig.inside()->geometry().center();
            Dune::FieldVector<DF,IG::dimension>
              outside_global = ig.geometry().center();
            inside_global -= outside_global;
            RF distance = inside_global.two_norm();
            r_s.accumulate(lfsu_s,0,(g*vn - D_inside*(g-x_s(lfsu_s,0))/distance)*face_volume);
            return;
          }
      }


      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
        ra.setTime(t);
      }


      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        tp.preStep(time,dt,stages);
      }

      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }


    private:
      static RA & raDefault()
      {
        static RA ra;
        return ra;
      }


      TP& tp;
      RA& ra;
      const bool central;
      typename TP::Traits::RangeFieldType time;

    };



    /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x,t) uv dx
     * \f}
     *
     * version where c(x,t) may become zero.
     */
    template<class TP>
    class ImplicitCCFVTemporalOperator
      : public NumericalJacobianVolume<ImplicitCCFVTemporalOperator<TP> >,
        public NumericalJacobianApplyVolume<ImplicitCCFVTemporalOperator<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      ImplicitCCFVTemporalOperator (TP& tp_)
        : tp(tp_)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        typedef LFSV Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        //typedef typename Space::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);


        // evaluate saturation at end of big step for cell activity
        typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local);
        bool active_cell = snew>0.;

        // evaluate capacity
        typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local);

        if (active_cell)
          {
            // residual contribution
            r.accumulate(lfsu,0,c*x(lfsu,0)*eg.geometry().volume());
          }

      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
      }


      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        tp.preStep(time,dt,stages);
        tp.setTimeTarget(time,dt);
      }

      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //suggest time step, asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return std::numeric_limits<typename TP::Traits::RangeFieldType>::max(); //initial value should be big enough
      }

    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
    };


    template<int b, typename T, typename U>
    struct selecti
    {
      typedef T type;
    };

    template<typename T, typename U>
    struct selecti<1, T, U>
    {
      typedef U type;
    };



    template<typename  TP, typename C>
    class FluxPostprocessing
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename TP::Traits::GridViewType,
                                                                               typename TP::Traits::RangeFieldType,
                                                                               TP::Traits::GridViewType::dimension,
                                                                               Dune::FieldVector<typename TP::Traits::RangeFieldType,TP::Traits::GridViewType::dimension> >,
                                              FluxPostprocessing<TP,C> >
    {
      // extract useful types
      typedef typename TP::Traits::GridViewType GV;
      typedef typename GV::Grid::ctype DF;
      typedef typename TP::Traits::RangeFieldType RF;
      typedef typename TP::Traits::RangeType RangeType;
      enum { dim = TP::Traits::GridViewType::dimension };
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;
      typedef typename TP::BCType BCType;

      typedef typename  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> RT0FEM;
      typedef typename  Dune::Q1LocalFiniteElement<DF,RF,1> P1FEM;

      typedef typename selecti<dim,RT0FEM,P1FEM>::type FEM;

      TP& tp;
      C& c;
      FEM rt0fe;
      const bool higherorder;
      typename TP::Traits::RangeFieldType time;


      typedef typename FEM::Traits::LocalBasisType::Traits::RangeType RT0RangeType;


    public:
      typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

      typedef Dune::PDELab::GridFunctionBase<Traits,FluxPostprocessing<TP,C> > BaseT;

      FluxPostprocessing ( TP& tp_, C& c_, const bool higherorder_=false) :
        tp(tp_), c(c_), higherorder(higherorder_) {}
      // set time where operator is to be evaluated (i.e. end of the time intervall)
      void set_time (typename TP::Traits::RangeFieldType time_)
      {
        time = time_;
      }

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        if (higherorder)
          {
            c.evaluate(e,x,y);
            return;
          }

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
          general(e.type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = e.geometry().global(inside_cell_center_local);

        // for coefficient computation
        RF vn[2*dim];    // normal velocities
        RF coeff[2*dim]; // RT0 coefficient
        Dune::FieldMatrix<typename Traits::DomainFieldType,dim,dim>
          B = e.geometry().jacobianInverseTransposed(x); // the transformation. Assume it is linear
        RF determinant = B.determinant();

        std::vector<typename Dune::FieldVector<DF,dim> > facecenters;
        facecenters.resize(2*dim);

        typename C::Traits::RangeType concentration_in = 0;
        c.evaluate(e,inside_cell_center_local,concentration_in);

        // loop over cell neighbors
        IntersectionIterator endit = c.getGridView().iend(e);
        for (IntersectionIterator iit = c.getGridView().ibegin(e); iit!=endit; ++iit)
          {
            // set to zero for processor boundary
            vn[iit->indexInInside()] = 0.0;

            // face geometry
            const Dune::FieldVector<DF,dim-1>&
              face_local = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);

            // face center in element coordinates
            Dune::FieldVector<DF,dim> face_center_in_element = iit->geometryInInside().global(face_local);

            facecenters[iit->indexInInside()]=face_center_in_element;
            // interior face
            if (iit->neighbor())
              {
                const Dune::FieldVector<DF,dim>&
                  outside_cell_center_local = Dune::ReferenceElements<DF,dim>::
                  general(iit->outside()->type()).position(0,0);
                Dune::FieldVector<DF,dim>
                  outside_cell_center_global = iit->outside()->geometry().global(outside_cell_center_local);

                // distance of cell centers
                Dune::FieldVector<DF,dim> d(outside_cell_center_global);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // convective flux
                RF u_upwind=0.0;
                typename C::Traits::RangeType concentration_out = 0;
                c.evaluate(*(iit->outside()), iit->geometryInOutside().center(),concentration_out);
                vn[iit->indexInInside()] = -(concentration_out-concentration_in)/distance;
              }

            // boundary face
            else if (iit->boundary())
              {
                // distance of cell center to boundary
                Dune::FieldVector<DF,dim> d = iit->geometry().global(face_local);
                d -= inside_cell_center_global;
                RF distance = d.two_norm();

                // evaluate boundary condition type
                BCType bctype = tp.bctype(*iit,face_local);

                if (bctype==BCType::Dirichlet) // Dirichlet boundary
                  {
                    typename TP::Traits::RangeFieldType g;
                    g=tp.g(*iit,face_local);

                    Dune::FieldVector<DF,dim>
                      inside_global = iit->inside()->geometry().center();
                    Dune::FieldVector<DF,dim>
                      outside_global = iit->geometry().center();
                    inside_global -= outside_global;
                    RF distance = inside_global.two_norm();
                    vn[iit->indexInInside()] = - (g-concentration_in)/distance;
                  }
              }

            // compute coefficient
               Dune::FieldVector<DF,dim> vstar=iit->unitOuterNormal(face_local); // normal on tranformef element
            vstar *= vn[iit->indexInInside()];
            Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
            if (iit->indexInInside()%2==0)
              normalhat[iit->indexInInside()/2] = -1.0;
            else
              normalhat[iit->indexInInside()/2] =  1.0;


            Dune::FieldVector<DF,dim> vstarhat(0);
            B.umtv(vstar,vstarhat); // Piola backward transformation
            vstarhat *= determinant;
            coeff[iit->indexInInside()] = vstarhat*normalhat;

          }

         RF coefficients[2*dim+1];    // coefficients for a+bx+cy
         std::vector<typename Traits::RangeType> fluxes; // fluxes which are computed before
         fluxes.resize(2*dim);

            B.invert();
            for (int i=0;i<2*dim;i++)
              {
                // compute velocity on reference element
                std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
                rt0fe.localBasis().evaluateFunction(facecenters[i],rt0vectors);
                typename Traits::RangeType yhat(0);
                for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
                  yhat.axpy(coeff[i],rt0vectors[i]);

                // apply Piola transformation

                fluxes[i]=0;
                B.umtv(yhat,fluxes[i]);
                fluxes[i] /= determinant;
              }

            RF f1 = vn[0]/determinant;
            RF f2 = -vn[1]/determinant;


            coefficients[2*dim]=concentration_in;
            for (int i=0;i<dim;i++)
              {
                coefficients[2*i]=fluxes[2*i];
                coefficients[2*i+1]=(coefficients[2*i]-fluxes[2*i+1])/2.0;

                if (dim==1)
                  {
                    coefficients[2*i]/=determinant;
                    coefficients[2*i+1]=(-fluxes[2*i+1]/determinant-coefficients[2*i])/2.0;
                  }
                coefficients[2*dim]-=(coefficients[2*i]/2.0+coefficients[2*i+1]/3.0);
              }


            y = coefficients[2*dim];
            for (int i=0;i<dim;i++)
              y+=coefficients[2*i]*x[i]+coefficients[2*i+1]*x[i]*x[i];

      }



      inline const typename Traits::GridViewType& getGridView () const
      {
        return c.getGridView();
      }

    private:

      template<typename T>
      T aavg (T a, T b) const
      {
        return 0.5*(a+b);
      }

      template<typename T>
      T havg (T a, T b) const
      {
        T eps = 1E-30;
        return 2.0/(1.0/(a+eps) + 1.0/(b+eps));
      }
    };


  }
}

#endif
