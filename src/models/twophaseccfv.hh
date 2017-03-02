// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TWOPHASEOP_HH
#define DUNE_PDELAB_TWOPHASEOP_HH

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
//#include<dune/pdelab/localoperator/defaultimp.hh>
#include"twosidedimp.hh"
#include<dune/pdelab/localoperator/idefault.hh>

#include"twophasetraits.hh"

namespace Dune {
  namespace PDELab {

    //! enum to choose between gas and liquid phase
    enum TwoPhasePhaseName
      {
        TwoPhaseLiquid = 0,
        TwoPhaseGas = 1
      };


    // a local operator for solving the two-phase flow in pressure-pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    // V  : Vector holding last time step
    template<typename TP>
    class TwoPhaseTwoPointFluxOperator
      : public NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseTwoPointFluxOperator<TP> >,

        public NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseTwoPointFluxOperator<TP> >,

        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { gas = TwoPhaseGas };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (const TP& tp_, Real scale_l_, Real scale_g_)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_)
      {}

      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (const TP& tp_)
        : tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g())
      {}

      // volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<liquid>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().volume();
        RF q_l = tp.q_l(eg.entity(),cell_center_local,time);
        RF q_g = tp.q_g(eg.entity(),cell_center_local,time);

        // contribution from source term
        r.accumulate(lfsv,liquid,-scale_l * q_l * cell_volume);
        r.accumulate(lfsv,gas,-scale_g * q_g * cell_volume);
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // cell geometries
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.outside()->type()).position(0,0);
        Dune::FieldVector<DF,IG::dimension>
          inside_cell_center_global = ig.inside()->geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_cell_center_global = ig.outside()->geometry().center();

        // distance of cell centers
        Dune::FieldVector<DF,dim> d(outside_cell_center_global);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);
        RF k_abs_outside = tp.k_abs(*(ig.outside()),outside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF rho_l_outside = tp.rho_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        // determines direction of liquid phase flow
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;
        RF pc_upwind, s_l_upwind_s, s_l_upwind_n, s_g_upwind_s, s_g_upwind_n;
        RF nu_l = aavg(tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid)),
                       tp.nu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid)));

        if (w_l>=0) // upwind capillary pressure on face
          pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
        else
          pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);

        s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
        s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);

        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * sigma_l
                       * w_l * face_volume);
        r_n.accumulate(lfsu_n,liquid,-scale_l * nu_l * sigma_l
                       * w_l * face_volume);

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas));
        // determines direction of gase phase flow
        RF w_g = (x_s(lfsu_s,gas)-x_n(lfsu_n,gas))/distance + aavg(rho_g_inside,rho_g_outside)*gn;
        RF nu_g = aavg(tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)),
                       tp.nu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas)));

        if (w_l*w_g<0) // new evaluation necessary only if signs differ
          {
            if (w_g>=0) // upwind capillary pressure on face
              pc_upwind = x_s(lfsu_s,gas)-x_s(lfsu_s,liquid);
            else
              pc_upwind = x_n(lfsu_n,gas)-x_n(lfsu_n,liquid);

            s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
            s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
          }

        // compute upwind value of gas saturation
        s_g_upwind_s = 1-s_l_upwind_s;
        s_g_upwind_n = 1-s_l_upwind_n;

        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_upwind_s)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_upwind_n)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,gas));
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

        r_s.accumulate(lfsu_s,gas,scale_g * nu_g * sigma_g * w_g * face_volume);
        r_n.accumulate(lfsu_n,gas,-scale_g * nu_g * sigma_g * w_g * face_volume);

#ifdef TEST_DEBUG_BODY
        std::cout << " nu_l " << nu_l
                  << "\n w_l " << w_l
                  << "\nsigma_l " << sigma_l
                  << "\nnu_l " << nu_l
                  << "\nw_g " << w_g
                  << "\nsigma_l " << sigma_g << std::endl;
        //std::cout << scale_l * nu_l * sigma_l  * w_l * face_volume << " " <<  scale_g * nu_g * sigma_g  * w_g * face_volume << std::endl;
#endif
      }

      // skeleton integral depending on test and ansatz functions
      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // face geometry
        const Dune::FieldVector<DF,dim-1>&
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        // if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = ig.inside()->geometry().global(inside_cell_center_local);

        // distance of cell center to boundary
        Dune::FieldVector<DF,dim> d = ig.geometry().global(face_local);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(*(ig.inside()),inside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * sigma_l
                           * w_l * face_volume);

#ifdef TEST_DEBUG_AC
            std::cout << "bc_l " << scale_l * (nu_l) * sigma_l * w_l * face_volume << std::endl;
#endif

#ifdef TEST_DEBUG
            std::cout << " pl" <<  x_s(lfsu_s,liquid)
                      << "\nrho " << rho_l_inside
                      << "\n g_l " << g_l
                      << "\n w_l " << w_l
                      << "\n s " << s_l
                      << "\n sigma " << sigma_l
                      <<"\n nu_l " << nu_l << std::endl;
#endif

          }

        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s(lfsu_s,gas)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas)-x_s(lfsu_s,liquid));
            RF s_g = 1-s_l;
            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            RF sigma_g = lambda_g_inside*k_abs_inside;
            RF nu_g = tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            r_s.accumulate(lfsu_s,gas,scale_g * nu_g * sigma_g
                           * w_g * face_volume);

#ifdef TEST_DEBUG_AC
            std::cout << "bc_g " << scale_g * nu_g * sigma_g * w_g * face_volume << std::endl;
#endif

#ifdef TEST_DEBUG
            std::cout << " pg " << x_s(lfsu_s,gas)
                      << "\nrho " << rho_g_inside
                      << "\n g_g " << g_g
                      << "\n w_g " << w_g
                      << "\n s " << s_l
                      << "\n sigma " << sigma_g
                      <<"\n nu_g " << nu_g << std::endl;
#endif
          }

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * j_l * face_volume);

#ifdef TEST_DEBUG_AC
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            std::cout << "bc_l_n " << scale_l* j_l * face_volume << std::endl;
#endif
          }

        // gas phase Neumann boundary
        if (bc_g==0)
          {
            RF nu_g = tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,gas));
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,gas,scale_g * nu_g * j_g * face_volume);

#ifdef TEST_DEBUG_AC
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            std::cout << "bc_g_n " << scale_g* j_g * face_volume << std::endl;
#endif
          }

      }

      // boundary integral independent of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
      {
        /*    // select the two components
              typedef typename LFSV::template Child<0>::Type PLSpace;

              // domain and range field type
              typedef typename PLSpace::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::DomainFieldType DF;
              typedef typename PLSpace::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeFieldType RF;
              typedef typename PLSpace::Traits::FiniteElementType::
              Traits::LocalBasisType::Traits::RangeType RangeType;

              // face geometry
              const Dune::FieldVector<DF,dim-1>&
              face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
              RF face_volume = ig.geometry().integrationElement(face_local)*
              Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

              // evaluate boundary condition type
              int bc_l = tp.bc_l(ig.intersection(),face_local,time);
              int bc_g = tp.bc_g(ig.intersection(),face_local,time);
              if (bc_l!=0 && bc_g!=0) return; // no Neumann boundary conditions

              // liquid phase Neumann boundary
              if (bc_l==0)
              {
              RF j_l = tp.j_l(ig.intersection(),face_local,time);
              r_s.accumulate(lfsv,liquid,scale_l * j_l * face_volume);

              }

              // gas phase Neumann boundary
              if (bc_g==0)
              {
              RF j_g = tp.j_g(ig.intersection(),face_local,time);
              r_s.accumulate(lfsv,gas,scale_g * j_g * face_volume);
              }
        */
      }

      //! set time for subsequent evaluation
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        tp.preStep(time,dt,stages);
      }

    private:
      const TP& tp;  // two phase parameter class
      typename TP::Traits::RangeFieldType time;
      Real scale_l, scale_g;

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


    /** a local operator for the storage operator
     *
     * \f{align*}{
     \int_\Omega c(x) uv dx
     * \f}
     */
    template<class TP>
    class TwoPhaseOnePointTemporalOperator
      : public NumericalJacobianVolume<TwoPhaseOnePointTemporalOperator<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseOnePointTemporalOperator<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { gas = TwoPhaseGas };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseOnePointTemporalOperator (const TP& tp_, Real scale_l_, Real scale_g_)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_)
      {
      }

      //! constructor: pass parameter object
      TwoPhaseOnePointTemporalOperator (const TP& tp_)
        : tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g())
      {}



      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<0>::Type PLSpace;

        // domain and range field type
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PLSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);
        RF cell_volume = eg.geometry().volume();

        RF phi = tp.phi(eg.entity(),cell_center_local);
        RF s_l = tp.s_l(eg.entity(),cell_center_local,
                        x(lfsu,gas)-x(lfsu,liquid));
        RF nu_l = tp.nu_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        RF nu_g = tp.nu_g(eg.entity(),cell_center_local,x(lfsu,gas));

        r.accumulate(lfsu,liquid,scale_l * phi * s_l * nu_l * cell_volume);
        r.accumulate(lfsu,gas,scale_g * phi * (1-s_l) * nu_g * cell_volume);

#ifdef TEST_DEBUG
        //std::cout <<scale_l * phi * s_l * nu_l * cell_volume << " " << scale_g * phi * (1-s_l) * nu_g * cell_volume << std::endl;
#endif
      }

      //! set time for subsequent evaluation
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
      }

      //! to be called once before each time step
      void preStep (typename TP::Traits::RangeFieldType time, typename TP::Traits::RangeFieldType dt,
                    int stages)
      {
        tp.preStep(time,dt,stages);
      }

    private:
      const TP& tp;
      typename TP::Traits::RangeFieldType time;
      Real scale_l, scale_g;
    };
  }
}

#endif
