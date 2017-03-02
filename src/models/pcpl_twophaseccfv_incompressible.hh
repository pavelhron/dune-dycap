
// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_PCPLTWOPHASEOP_INCOMPRESSIBLE_HH
#define DUNE_PDELAB_PCPLTWOPHASEOP_INCOMPRESSIBLE_HH

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include"defaultimp.hh"
#include<dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include"twophasetraits.hh"


namespace Dune {
  namespace PDELab {


    //! enum to choose between gas and liquid phase
    enum TwoPhasePhaseName
      {
        TwoPhaseLiquid = 0,
        TwoPhaseCapillary = 1
      };

    /*! \brief a local operator for solving the two-phase flow problem: \n

      \f{align*}{
      \partial_t(\phi S_l \nu_l) + \nabla\cdot\{ \nu_l u_l \}  = q_l, \\
      \partial_t(\phi S_g \nu_g) + \nabla\cdot\{ \nu_g u_g \} = q_g,  \\
      \f}

      The phase velocity is related to the phase pressure via the extended Darcy law
      \f{align*}
      u_l = -\frac{k_{rl}(s_l)}{\mu_l} K \left(\nabla p_l - \rho_l g\right)
      u_\g = -\frac{k_{rg}(s_g)}{\mu_g} K \left(\nabla p_g - \rho_g g\right)
      \f}

    */



    // a local operator for solving the two-phase flow in liquid pressure - capillary pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
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
      enum { capillary = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;
      typedef typename TP::Traits::GridViewType GV;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume   = true };
      //  enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (TP& tp_, Real scale_l_, Real scale_g_, bool blend_, Real blendmin_, Real blendmax_)
        : NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), tp(tp_), scale_l(scale_l_), scale_g(scale_g_), time(0.), verbosity(0)
      {}


      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (TP& tp_, Real scale_l_, Real scale_g_)
        : NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), tp(tp_), scale_l(scale_l_), scale_g(scale_g_),time(0.), verbosity(0)
      {}

      //! constructor: pass parameter object
      TwoPhaseTwoPointFluxOperator (TP& tp_)
        : NumericalJacobianSkeleton<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTwoPointFluxOperator<TP> >(eps_numdiff), tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g()),time(0.), verbosity(0)
      {
      }


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

        r.accumulate(lfsv,liquid,-scale_l * q_l * cell_volume);
        r.accumulate(lfsv,capillary,-scale_g * q_g * cell_volume);

        if (verbosity>1)
          std::cout << -scale_l * q_l * cell_volume << " " << -scale_g * q_g * cell_volume << std::endl;

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
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside().type()).position(0,0);
        const Dune::FieldVector<DF,dim>&
          outside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.outside().type()).position(0,0);
        Dune::FieldVector<DF,IG::dimension>
          inside_cell_center_global = ig.inside().geometry().center();
        Dune::FieldVector<DF,IG::dimension>
          outside_cell_center_global = ig.outside().geometry().center();

        // distance of cell centers
        Dune::FieldVector<DF,dim> d(outside_cell_center_global);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();


        int verbosity_level=0;
        ElementMapper<GV> cell_mapper(tp.getGridView());
        //typename GV::IndexSet::IndexType ids = cell_mapper.map(*(ig.inside()));
        //typename GV::IndexSet::IndexType idn = cell_mapper.map(*(ig.outside()));

        //        if (ids<100 || idn < 100)
        //  verbosity_level=2;


        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<DF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(ig.inside(),inside_cell_center_local);
        RF k_abs_outside = tp.k_abs(ig.outside(),outside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF rho_l_outside = tp.rho_l(ig.outside(),outside_cell_center_local,x_n(lfsu_n,liquid));
        // determines liquid flow direction
        const RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF rho_g_outside = tp.rho_g(ig.outside(),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        // determines gas flow direction
        RF w_g = ((x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid))-(x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid)))/distance + aavg(rho_g_inside,rho_g_outside)*gn;

        // first equation
        RF pc_upwind, s_l_upwind_s, s_g_upwind_s, s_l_upwind_n, s_g_upwind_n;
        RF nu_l = aavg(tp.nu_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid)),
                       tp.nu_l(ig.outside(),outside_cell_center_local,x_n(lfsu_n,liquid)));

        // upwind capillary pressure on face
        if (w_l>=0)
          // flow from the current element
          pc_upwind = x_s(lfsu_s,capillary);
        else
          // flow to the current element
          pc_upwind = x_n(lfsu_n,capillary);

        RF seps = 0.;

        // compute upwind value of liquid saturation
        s_l_upwind_s = tp.s_l(ig.inside(),inside_cell_center_local,pc_upwind);
        s_l_upwind_n = tp.s_l(ig.outside(),outside_cell_center_local,pc_upwind);

        s_l_upwind_s = std::min(1.0,s_l_upwind_s+seps);
        s_l_upwind_n = std::min(1.0,s_l_upwind_n+seps);

        // harmonic average of k_r
        RF lambda_l_inside = tp.kr_l(ig.inside(),inside_cell_center_local,s_l_upwind_s)/
          tp.mu_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(ig.outside(),outside_cell_center_local,s_l_upwind_n)/
          tp.mu_l(ig.outside(),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        // residual accumulation
        r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * sigma_l
                       * w_l * face_volume);
        r_n.accumulate(lfsu_n,liquid,-scale_l * nu_l * sigma_l
                       * w_l * face_volume);


        int precision = 15;
        // remember old flags
        std::ios_base::fmtflags oldflags = std::cout.flags();

        // set the output format
        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        int oldprec = std::cout.precision();
        std::cout.precision(precision);


        if (verbosity_level)
          std::cout << "l " << scale_l * nu_l * sigma_l
            * w_l * face_volume<< " " <<  w_l << " " << distance << " " <<x_s(lfsu_s,liquid) << " " << x_n(lfsu_n,liquid)<< " "<< (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance << " " <<  aavg(rho_l_inside,rho_l_outside)*gn<< std::endl;

        // second equation
        RF nu_g = aavg(tp.nu_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid)),
                       tp.nu_g(ig.outside(),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid)));

        // new evaluation necessary only if signs differ
        if (w_l*w_g<=0)
          {
            // upwind capillary pressure on face
            if (w_g>=0)
              pc_upwind = x_s(lfsu_s,capillary);
            else
              pc_upwind = x_n(lfsu_n,capillary);

            //RF bl = Blend(std::abs(w_g), 0.0, 1.e-2);
            // pc_upwind = bl * pc_upwind + (1.- bl) * aavg (x_s(lfsu_s,capillary),x_n(lfsu_n,capillary));

            // compute upwind value of liquid saturation if necessary
            s_l_upwind_s = tp.s_l((ig.inside()),inside_cell_center_local,pc_upwind);
            s_l_upwind_n = tp.s_l((ig.outside()),outside_cell_center_local,pc_upwind);

            s_l_upwind_s = std::min(1.0,s_l_upwind_s+seps);
            s_l_upwind_n = std::min(1.0,s_l_upwind_n+seps);
          }
        // compute upwind value of gas saturation
        // here could be problem if there is no air phas
        s_g_upwind_s = 1.-s_l_upwind_s;
        s_g_upwind_n = 1.-s_l_upwind_n;

        RF lambda_g_inside = tp.kr_g((ig.inside()),inside_cell_center_local,s_g_upwind_s)/
          tp.mu_g((ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF lambda_g_outside = tp.kr_g((ig.outside()),outside_cell_center_local,s_g_upwind_n)/
          tp.mu_g((ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

        //      if ((scale_g * nu_g * sigma_g * w_g * face_volume)>1.e20)
        //        std::cout << nu_g << " " << sigma_g << " " << w_g << std::endl;

        r_s.accumulate(lfsu_s,capillary,scale_g * nu_g * sigma_g * w_g * face_volume);
        r_n.accumulate(lfsu_n,capillary,-scale_g * nu_g * sigma_g * w_g * face_volume);


        if (verbosity_level)
          std::cout <<"g " << scale_g * nu_g * sigma_g * w_g * face_volume << " " <<  sigma_g  << " " << s_g_upwind_s << " " << s_g_upwind_n << " " <<  lambda_g_inside << " " << lambda_g_outside<< std::endl;

        std::cout.precision(oldprec);

        // reset the output format
        std::cout.flags(oldflags);

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
        //typedef typename PLSpace::Traits::FiniteElementType::
        //Traits::LocalBasisType::Traits::RangeType RangeType;

        // face geometry
        const Dune::FieldVector<DF,dim-1>&
          face_local = Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        //  if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside().type()).position(0,0);
        Dune::FieldVector<DF,dim>
          inside_cell_center_global = ig.inside().geometry().global(inside_cell_center_local);

        // distance of cell center to boundary
        Dune::FieldVector<DF,dim> d = ig.geometry().global(face_local);
        d -= inside_cell_center_global;
        RF distance = d.two_norm();

        // absolute permeability
        RF k_abs_inside = tp.k_abs(ig.inside(),inside_cell_center_local);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);


        // liquid phase Dirichlet boundary
        //
        if (bc_l==1)
          {
            RF rho_l_inside = tp.rho_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF lambda_l_inside = tp.kr_l(ig.inside(),inside_cell_center_local,s_l)/
              tp.mu_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            RF nu_l = tp.nu_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
            r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * sigma_l
                           * w_l * face_volume);

            if (verbosity>1)
              std::cout <<scale_l * nu_l * sigma_l * w_l * face_volume << std::endl;
          }

        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF rho_g_inside = tp.rho_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF s_g = 1-s_l;

            //  s_g= std::max(0.0,s_g-0.01);

            RF lambda_g_inside = tp.kr_g(ig.inside(),inside_cell_center_local,s_g)/
              tp.mu_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF sigma_g = lambda_g_inside*k_abs_inside;
            RF nu_g =tp.nu_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            r_s.accumulate(lfsu_s,capillary,scale_g * nu_g * sigma_g
                           * w_g * face_volume);

            if (verbosity>1)
              std::cout << scale_g * nu_g * sigma_g * w_g * face_volume<< std::endl;
          }

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF nu_l = tp.nu_l(ig.inside(),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,liquid,scale_l * nu_l * j_l * face_volume);
            if (verbosity>1)
              std::cout << scale_l * nu_l * j_l * face_volume << std::endl;
          }

        // gas phase Neumann boundary
        if (bc_g==0)
          {
            RF nu_g =tp.nu_g(ig.inside(),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,capillary,scale_g * nu_g * j_g * face_volume);
            if (verbosity>1)
              std::cout << scale_g * nu_g * j_g * face_volume << std::endl;
          }
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
      TP& tp;  // two phase parameter class
      Real scale_l, scale_g;
      bool blend;
      Real blendmin,blendmax;
      typename TP::Traits::RangeFieldType time;
      const int verbosity;

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


    /*! \brief a local operator for solving the temporal term of the two-phase flow problem: \n

      \f{align*}{
      \partial_t(\phi S_l \nu_l), \     \
      \partial_t(\phi S_g \nu_g),  \    \
      \f}

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
      enum { capillary = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseOnePointTemporalOperator (TP& tp_, Real scale_l_, Real scale_g_)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_), time(0.), verbosity(0)
      {
      }

      //! constructor: pass parameter object
      TwoPhaseOnePointTemporalOperator (TP& tp_)
        : tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g()), time(0.), verbosity(1)
      {
      }

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
        RF s_l = tp.s_l(eg.entity(),cell_center_local,x(lfsu,capillary));
        RF s_g = 1. -s_l;
        RF nu_l = tp.nu_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        RF nu_g = tp.nu_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        r.accumulate(lfsu,liquid,scale_l * phi * s_l * nu_l * cell_volume);
        r.accumulate(lfsu,capillary,scale_g * phi * s_g * nu_g * cell_volume);
        if (verbosity>1)
          std::cout << "t " << scale_l * phi * s_l * nu_l * cell_volume << " " << scale_g * phi * s_g * nu_g * cell_volume << std::endl;
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
      TP& tp;
      Real scale_l, scale_g;
      typename TP::Traits::RangeFieldType time;
      const int verbosity;

    };






  }
}

#endif
