// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_COMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_COMPONENTTRANSPORTOP_HH

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include "fluxreconstruction.hh"

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
    template<typename TP, typename FR=DefaultFluxReconstruction<typename TP::Traits> >
    class ModifiedCCFVSpatialTransportOperator :
      //  public NumericalJacobianSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      //  public NumericalJacobianBoundary<ModifiedCCFVSpatialTransportOperator<TP> >,
      //  public NumericalJacobianApplySkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      //  public NumericalJacobianApplyBoundary<ModifiedCCFVSpatialTransportOperator<TP> >,
      //  public NumericalJacobianApplyVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      public NumericalJacobianVolumePostSkeleton<ModifiedCCFVSpatialTransportOperator<TP> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      typedef FR FluxReconstruction;

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
        : tp(tp_), flux(fluxDefault()), zero(z), neighbor_computation(neighbor_computation_)
      {
      }

      ModifiedCCFVSpatialTransportOperator (TP& tp_, FluxReconstruction & f, typename TP::Traits::RangeFieldType z=1e-7, bool neighbor_computation_=false)
        : tp(tp_), flux(f), zero(z), neighbor_computation(neighbor_computation)
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

        // skip if the neighbor is not active
        if (neighbor_computation)
          {
            bool active_neighbor = tp.snew(*(ig.outside()),outside_local)>zero;

            if (!active_neighbor)
              return;
          }

        // convective flux
        RF u_upwind=0.0;
        if (vn>=0)
          {
            // evaluate face flux on inside
            flux.evaluate(*ig.inside(), ig.geometryInInside().center(), x_s(lfsu_s,0), u_upwind);
          }
        else
          {
            // evaluate face flux on outside
            flux.evaluate(*ig.outside(), ig.geometryInOutside().center(), x_n(lfsu_n,0), u_upwind);
          }

        r_s.accumulate(lfsu_s,0,(u_upwind*vn)*face_volume);

        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

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
        r_s.accumulate(lfsu_s,0,-(D_avg*(x_n(lfsu_n,0)-x_s(lfsu_s,0))/distance)*face_volume);

        celloutflux += D_avg*face_volume/distance;
      }



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
        int bc = tp.bc(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bc==0) // Neumann boundary
          {
            typename TP::Traits::RangeFieldType j = tp.j(ig.intersection(),face_local);
            r_s.accumulate(lfsu_s,0,j*face_volume);
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        vmax = std::max(vmax, std::abs(vn));

        if (vn>=0)
          celloutflux += vn*face_volume; // dt computation

        if (bc==2) // Outflow boundary
          {
            r_s.accumulate(lfsu_s,0,vn*x_s(lfsu_s,0)*face_volume);
            return;
          }

        if (bc==1) // Dirichlet boundary
          {
            typename TP::Traits::RangeFieldType g;
            //  if (vn>0)
            //     g=x_s(lfsu_s,0);
            //< else
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
            celloutflux += D_inside*face_volume/distance;
            return;
          }
      }

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
        int bc = tp.bc(ig.intersection(),face_local);

        // do things depending on boundary condition type
        if (bc==0) // Neumann boundary
          {
            return;
          }

        // evaluate velocity
        typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element));

        // the normal velocity
        RF vn = v*ig.centerUnitOuterNormal();
        if (bc==2) // Outflow boundary
          {
            mat_ss.accumulate(lfsu_s,0,lfsu_s,0,vn*face_volume);
            return;
          }

        if (bc==1) // Dirichlet boundary
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

      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSU::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        const int dim = EG::Geometry::dimension;

        if (!first_stage) return; // time step calculation is only done in first stage

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);


        // compute optimal dt for this cell
        typename TP::Traits::RangeFieldType cellcapacity = tp.c(eg.entity(),inside_local)*eg.geometry().volume();
        typename TP::Traits::RangeFieldType celldt = cellcapacity/(celloutflux+1E-40);
        //std::cout << cellcapacity << " " << celloutflux << " " << celldt << std::endl;

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

        Dune::FieldVector<DF, dim>
          global = eg.geometry().global(inside_local);

        //  std::cout << global[0] << " " << global[1] << " " << inside_local[0] << std::endl;
        //  if (active_cell && global[1]<0.2)
        dtmin = std::min(dtmin,celldt);
      }


      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef typename LFSV::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        const int dim = EG::Geometry::dimension;

        // no sources if cell is not active !
        if (!active_cell) return;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        // evaluate source term
        typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local);

        r.accumulate(lfsv,0,-q*eg.geometry().volume());
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        flux.setTime(t);
        tp.setTime(t);
        //std::cout << "set time spatial " << t << std::endl;
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
      : // public NumericalJacobianVolume<ModifiedCCFVTemporalOperator<TP> >,
      // public NumericalJacobianApplyVolume<ModifiedCCFVTemporalOperator<TP> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      ModifiedCCFVTemporalOperator (TP& tp_)
        : tp(tp_), zero(1e-7)
      {
      }

      ModifiedCCFVTemporalOperator (TP& tp_, typename TP::Traits::RangeFieldType z)
        : tp(tp_), zero(z)
      {
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
        // std::cout << "set time temporal" << t << std::endl;
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

#ifdef COMPTRANS_DEBUG
        if (active_cell)
          std::cout << "M: time=" << time
                    << " pos=" << eg.geometry().center()
                    << " snew=" << snew
                    << " capacityvol=" << c*eg.geometry().volume()
                    << " residual=" << r[0]
                    << std::endl;
#endif // COMPTRANS_DEBUG
      }

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

#ifdef COMPTRANS_DEBUG
        if (active_cell)
          std::cout << "J: time=" << time
                    << " pos=" << eg.geometry().center()
                    << " snew=" << snew
                    << " mat=" << mat(0,0)
                    << std::endl;
#endif // COMPTRANS_DEBUG
      }

      //suggest time step, asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return 1e100; //initial value should be big enough
      }

    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
      typename TP::Traits::RangeFieldType zero;
    };

  }
}

#endif
