// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_MULTICOMPONENTTRANSPORTOP_HH
#define DUNE_PDELAB_MULTICOMPONENTTRANSPORTOP_HH

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include"transporttraits.hh"
#include "fluxreconstruction.hh"

namespace Dune {
  namespace PDELab {


    /** a local operator for a cell-centered finite folume scheme for
        the diffusion-reaction equation

        \nabla \cdot \{v u- D \nabla u \} = q in \Omega
        u = g on \Gamma_D
        \{v u - D \nabla u \} \cdot \nu = j on \Gamma_N
        outflow on \Gamma_O

        Modified version for the case

        d_t (c(x,t)u(x,t)) + \nabla \cdot \{v u - D \nabla u \} = q in \Omega

        where c(x,t) may become zero. We assume that the following holds:

        c(x,t+dt) <= eps  ==>  c(x,t) <= eps



        \tparam TP  parameter class implementing ComponentTransportParameterInterface
    */
    template<typename TP, typename RA = ReactionBaseAdapter>
    class MulticomponentCCFVSpatialTransportOperator :
      public NumericalJacobianSkeleton<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplySkeleton<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,

      public NumericalJacobianBoundary<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplyBoundary<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,

      public NumericalJacobianVolume<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,
      public NumericalJacobianApplyVolume<MulticomponentCCFVSpatialTransportOperator<TP,RA> >,

      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };

      enum { dim = TP::Traits::GridViewType::dimension };
      typedef typename TP::BCType BCType;

      MulticomponentCCFVSpatialTransportOperator (TP& tp_, typename TP::Traits::RangeFieldType z=1e-7)
        : tp(tp_), zero(z), ra(raDefault())
      {
      }

      MulticomponentCCFVSpatialTransportOperator (TP& tp_, RA& ra_, typename TP::Traits::RangeFieldType z=1e-7)
        : tp(tp_), zero(z), ra(ra_)
      {
      }


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {

        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);


        // here could described the source (right hand side of the sytem of equations)
        // for each equation (big ammount of code)
        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {
            // evaluate source term
            typename TP::Traits::RangeFieldType q = tp.q(eg.entity(),inside_local,k);
            r.accumulate(lfsv,k,-q*eg.geometry().volume());
          }

        // evaluate the reaction adapter which is done on the local basis
        // if we are using CCFV, then x contains DOF's corresponding to the entity eg
        // sorting of x depends on the powergridfunction structure and ordering
        // in this case it is sorted lexikographically
        ra.evaluate(eg.entity(),lfsu,x,lfsv,r);
      }



      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);



        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            // evaluate boundary condition type
            BCType bctype = tp.bctype(ig.intersection(),face_local,k);

            // do things depending on boundary condition type
            if (bctype==BCType::Neumann) // Neumann boundary
              {
                typename TP::Traits::RangeFieldType j = tp.j(ig.intersection(),face_local,k);
                r_s.accumulate(lfsu_s,k,j*face_volume);
                continue;
              }

            // evaluate velocity
            typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element,k));

            // the normal velocity
            RF vn = v*ig.centerUnitOuterNormal();

            if (bctype==BCType::Outflow) // Outflow boundary
              {
                r_s.accumulate(lfsu_s,k,vn*x_s(lfsu_s,k)*face_volume);
                continue;
              }

             if (bctype==BCType::Dirichlet) // Dirichlet boundary
              {
                typename TP::Traits::RangeFieldType g;
                g=tp.g(ig.intersection(),face_local,k);
                r_s.accumulate(lfsu_s,k,(g*vn - D_inside*(g-x_s(lfsu_s,k))/distance)*face_volume);
                continue;
              }
          }

      }


      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {

        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

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

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {
            // evaluate velocity
            typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element,k));

            // the normal velocity
            RF vn = v*ig.centerUnitOuterNormal();

            // convective flux
            RF u_upwind=0.0;
            if (vn>=0) u_upwind = x_s(lfsu_s,k); else u_upwind = x_n(lfsu_n,k);
            // central differences for concentration
            if (std::abs(vn)<1.e-8)
              u_upwind = (x_n(lfsu_n,k)+x_s(lfsu_s,k))/2.0;

            r_s.accumulate(lfsu_s,k,(u_upwind*vn)*face_volume);
            r_n.accumulate(lfsu_n,k,-(u_upwind*vn)*face_volume);

            // evaluate diffusion coefficients
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local,k);
            typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));


            // diffusive flux
            r_s.accumulate(lfsu_s,k,-(D_avg*(x_n(lfsu_n,k)-x_s(lfsu_s,k))/distance)*face_volume);
            r_n.accumulate(lfsu_n,k,(D_avg*(x_n(lfsu_n,k)-x_s(lfsu_s,k))/distance)*face_volume);
          }

      }
      /*
      // jacobian of skeleton term
      // each face is only visited TWICE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_skeleton (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                              M& mat_ss, M& mat_sn,
                              M& mat_ns, M& mat_nn) const
      {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;


        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        const Dune::FieldVector<DF,IG::dimension>&
          inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,IG::dimension>&
          outside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.outside()->type()).position(0,0);

        // distance between cell centers in global coordinates
        Dune::FieldVector<DF,IG::dimension>
          inside_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_global = ig.outside()->geometry().center();
        inside_global -= outside_global;
        RF distance = inside_global.two_norm();

        // face center in element coordinates
        Dune::FieldVector<DF,IG::dimension> face_center_in_element = ig.geometryInInside().global(face_local);

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {

            // evaluate saturation at end of big step for cell activity
            //typename TP::Traits::RangeFieldType snew = tp.snew(*(ig.inside()),face_center_in_element,k);
            //bool active_cell = snew>zero;
            //if (!active_cell) continue;
            // evaluate velocity
            typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element,k));

            // the normal velocity
            RF vn = v*ig.centerUnitOuterNormal();


            // activity of neighbor
            bool active_neighbor = tp.snew(*(ig.outside()),outside_local,k)>zero;
            // skip if the neighbor is not active
            //if(!active_neighbor) continue;

            // convective flux
            if (vn>=0)
              {
                mat_ss.accumulate(lfsu_s,k,lfsu_s,k,vn*face_volume);
                if (active_neighbor)
                  mat_ns.accumulate(lfsu_n,k,lfsu_s,k,-vn*face_volume);
              }
            else
              {
                if(active_neighbor)
                  {
                    mat_nn.accumulate(lfsu_n,k,lfsu_n,k,vn*face_volume);
                    mat_sn.accumulate(lfsu_s,k,lfsu_n,k,-vn*face_volume);
                  }
              }

            // evaluate diffusion coefficients
            typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
            typename TP::Traits::RangeFieldType D_outside = tp.D(*(ig.outside()),outside_local,k);
            typename TP::Traits::RangeFieldType D_avg = 2.0/(1.0/(D_inside+1E-40) + 1.0/(D_outside+1E-40));



            // diffusive flux, in self
            mat_ss.accumulate(lfsu_s,k,lfsu_s,k,D_avg/distance*face_volume);

            if(active_neighbor)
              {
                // in self
                mat_ns.accumulate(lfsu_n,k,lfsu_s,k,D_avg/distance*face_volume);

                // jiggle in neighbor
                mat_sn.accumulate(lfsu_s,k,lfsu_n,k,D_avg/distance*face_volume);
                mat_nn.accumulate(lfsu_n,k,lfsu_n,k,D_avg/distance*face_volume);
              }
          }

      }


      // jacobian of boundary term
      template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_boundary (const IG& ig,
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              M& mat_ss) const
      {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;



        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();
        Dune::FieldVector<DF,dim> face_center_in_element = ig.geometryInInside().global(face_local);

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {

            // evaluate saturation at end of big step for cell activity
            typename TP::Traits::RangeFieldType snew = tp.snew(*(ig.inside()),face_center_in_element,k);
            bool active_cell = snew>zero;
            if (!active_cell) continue;
            // evaluate boundary condition type
            BCType bctype = tp.bctype(ig.intersection(),face_local,k);

            // do things depending on boundary condition type
             if (bctype==BCType::Neumann) // Neumann boundary
              {
                return;
              }

            // evaluate velocity
            typename TP::Traits::RangeType v(tp.v(*(ig.inside()),face_center_in_element,k));

            // the normal velocity
            RF vn = v*ig.centerUnitOuterNormal();
             if (bctype==BCType::Outflow) // Outflow boundary
              {
                mat_ss.accumulate(lfsu_s,k,lfsu_s,k,vn*face_volume);
                return;
              }

             if (bctype==BCType::Dirichlet) // Dirichlet boundary
              {
                if (vn>=0)
                  mat_ss.accumulate(lfsu_s,k,lfsu_s,k,vn*face_volume);
                const Dune::FieldVector<DF,IG::dimension>&
                  inside_local = Dune::ReferenceElements<DF,IG::dimension>::general(ig.inside()->type()).position(0,0);
                typename TP::Traits::RangeFieldType D_inside = tp.D(*(ig.inside()),inside_local,k);
                Dune::FieldVector<DF,IG::dimension>
                  inside_global = ig.inside()->geometry().center();
                Dune::FieldVector<DF,IG::dimension>
                  outside_global = ig.geometry().center();
                inside_global -= outside_global;
                RF distance = inside_global.two_norm();
                mat_ss.accumulate(lfsu_s,k,lfsu_s,k,D_inside/distance*face_volume);
                return;
              }
          }
      }
      */

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
      typename TP::Traits::RangeFieldType zero;
      typename TP::Traits::RangeFieldType time;
      RA& ra;
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
    class MulticomponentCCFVTemporalOperator
      : public NumericalJacobianVolume<MulticomponentCCFVTemporalOperator<TP> >,
      public NumericalJacobianApplyVolume<MulticomponentCCFVTemporalOperator<TP> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      MulticomponentCCFVTemporalOperator (TP& tp_, typename TP::Traits::RangeFieldType z = 1.e-7)
        : tp(tp_), zero(z)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        // typedef typename Space::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeFieldType RF;
        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {

            // evaluate saturation at end of big step for cell activity
            typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local,k);
            bool active_cell = snew>zero;

            // evaluate capacity
            typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local,k);

            if (active_cell)
              {
                // residual contribution
                r.accumulate(lfsu,k,c*x(lfsu,k)*eg.geometry().volume());
              }
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M& mat) const
      {
        typedef typename LFSV::template Child<0>::Type Space;

        // domain and range field type
        typedef typename Space::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;

        const int dim = EG::Geometry::dimension;

        // cell center
        const Dune::FieldVector<DF,dim>&
          inside_local = Dune::ReferenceElements<DF,dim>::general(eg.entity().type()).position(0,0);

        for (std::size_t k = 0;k<TP::COMPONENTS;k++)
          {

            // evaluate saturation at end of big step for cell activity
            typename TP::Traits::RangeFieldType snew = tp.snew(eg.entity(),inside_local,k);
            bool active_cell = snew>zero;

            // evaluate capacity
            typename TP::Traits::RangeFieldType c = tp.c(eg.entity(),inside_local,k);

            if (active_cell)
              {
                mat.accumulate(lfsu,k,lfsu,k,c*eg.geometry().volume());
              }
            else
              mat.accumulate(lfsu,k,lfsu,k,1.0);
          }
      }

      /*
      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        //tp.preStep(time,dt,stages);
        //tp.setTimeTarget(time,dt);
      }

      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }
      */
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
      typename TP::Traits::RangeFieldType zero;
      typename TP::Traits::RangeFieldType time;
    };

  }
}

#endif
