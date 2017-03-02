
// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_PCPLTWOPHASEOP_TOTALFLUX_HH
#define DUNE_PDELAB_PCPLTWOPHASEOP_TOTALFLUX_HH

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/idefault.hh>

#include"twophasetraits.hh"


namespace Dune {
  namespace PDELab {


    //! enum to choose between gas and liquid phase
    enum TwoPhasePhaseName
      {
        TwoPhaseLiquid = 0,
        TwoPhaseCapillary = 1
      };


    // a local operator for solving the two-phase flow in liquid pressure - capillary pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    template<typename TP>
    class TwoPhaseTotalFlux
      : public NumericalJacobianVolume<TwoPhaseTotalFlux<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseTotalFlux<TP> >,
        public NumericalJacobianSkeleton<TwoPhaseTotalFlux<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseTotalFlux<TP> >,
        public NumericalJacobianBoundary<TwoPhaseTotalFlux<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseTotalFlux<TP> >,


        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doAlphaVolume   = true };
      //  enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseTotalFlux (const TP& tp_, Real scale_l_, Real scale_g_)
        : NumericalJacobianSkeleton<TwoPhaseTotalFlux<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTotalFlux<TP> >(eps_numdiff), tp(tp_), scale_l(scale_l_), scale_g(scale_g_),time(0.)
      {}

      //! constructor: pass parameter object
      TwoPhaseTotalFlux (const TP& tp_)
        : NumericalJacobianSkeleton<TwoPhaseTotalFlux<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTotalFlux<TP> >(eps_numdiff), tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g()),time(0.)
      {
      }


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
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
        // liquid phase calculation
        RF rho_l = tp.rho_l(eg.entity(),cell_center_local,x(lfsu,liquid));

        // gas phase calculation
        RF rho_g = tp.rho_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        RF q_l = tp.q_l(eg.entity(),cell_center_local,time);
        RF q_g = tp.q_g(eg.entity(),cell_center_local,time);

        // contribution from source term
        r.accumulate(lfsv,flux,  -(q_l/rho_l+q_g/rho_g) * cell_volume);
        r.accumulate(lfsv,transport, - q_l/rho_l * cell_volume);
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
        RF k_abs = aavg(k_abs_inside, k_abs_outside);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF rho_l_outside = tp.rho_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF rho_l =  aavg(rho_l_inside,rho_l_outside);
        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_s(lfsu_s,liquid));
        RF rho_g = aavg(rho_g_inside,rho_g_outside);

        // compute upwind value of liquid saturation
        RF s_l_s = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
        RF s_l_n = tp.s_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_s,capillary));

        RF s_g_s = 1. - s_l_s;
        RF s_g_n = 1. - s_l_n;


        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_s)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_n)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        RF lambda_g =  havg(lambda_g_inside,lambda_g_outside);


        // harmonic average of k_r
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF lambda_l =  havg(lambda_l_inside,lambda_l_outside);

        //  - grad p_l + rho_l * g * n,e = - (grad p_l - rho_l * g * n,e)
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;
        //  - grad p_l - grad p_c + rho_l * g * n,e = - (grad p_l + grad p_c - rho_g * g * n,e)
        RF w_g = ((x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid))-(x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid)))/distance + aavg(rho_g_inside,rho_g_outside)*gn;

        RF nu_g = aavg(tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid)), tp.nu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid)));

        RF nu_l = aavg(tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid)), tp.nu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid)));

        //  RF sum = rho_l*lambda_l + rho_g*lambda_g;

        //  RF G = (rho_l*rho_l*lambda_l + rho_g*rho_g*lambda_g)/sum;

        RF qt = 0.;

        RF pc_upwind, s_l_upwind_s, s_g_upwind_s, s_l_upwind_n, s_g_upwind_n;
        // upwind capillary pressure on face
        if (w_l>=0)
          // flow from the current element
          pc_upwind = x_s(lfsu_s,capillary);
        else
          // flow to the current element
          pc_upwind = x_n(lfsu_n,capillary);

        // compute upwind value of liquid saturation
        s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
        s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);

        // harmonic average of k_r
        RF sigma_l = havg(tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind_s)/
                          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid))*k_abs_inside,tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind_n)/
                          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid))*k_abs_outside);

        qt-=sigma_l * w_l;

        // new evaluation necessary only if signs differ
        if (w_l*w_g<=0)
          {
            // upwind capillary pressure on face
            if (w_g>=0)
              pc_upwind = x_s(lfsu_s,capillary);
            else
              pc_upwind = x_n(lfsu_n,capillary);

            // compute upwind value of liquid saturation if necessary
            s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
            s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);
          }
        // compute upwind value of gas saturation
        s_g_upwind_s = 1-s_l_upwind_s;
        s_g_upwind_n = 1-s_l_upwind_n;

        RF sigma_g = havg(tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_upwind_s)/
                          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid))*k_abs_inside,tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_upwind_n)/
                          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid))*k_abs_outside);

        qt-=sigma_g * w_g;


        //      qt = sum * k_abs *( (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid) )/distance + rho_g*lambda_g/sum*(x_s(lfsu_s,capillary)-x_n(lfsu_n,capillary) )/distance + G*gn );

        // - (grad p_c + (rho_l- rho_g) * g * n,e)
        RF v = w_g - w_l;

        // Gl = - rho_n * lambda_n * v
        RF Gl = -rho_g * v;
        // Gg = rho_l * lambda_l * v
        RF Gg = rho_l* v;

        RF lambda_l_tilde = 0.0;
        RF lambda_g_tilde = 0.0;
        //  RF fluxl(0), fluxg(0);

        if ((Gl*qt)>=0) // alpa = liquid
          {
            if (qt>=0)
              lambda_l_tilde = lambda_l_inside;
            else
              lambda_l_tilde = lambda_l_outside;

            Gg = rho_l*lambda_l_tilde*v*k_abs;

            if ((qt+Gg)>=0)
              lambda_g_tilde = lambda_g_inside;
            else
              lambda_g_tilde = lambda_g_outside;

            Gl = -rho_g*lambda_g_tilde*v*k_abs;
          }
        else
          {
            if (qt>=0)
              lambda_g_tilde = lambda_g_inside;
            else
              lambda_g_tilde = lambda_g_outside;

            Gl = - rho_g*lambda_g_tilde*v*k_abs;

            if ((qt+Gl)>=0)
              lambda_l_tilde = lambda_l_inside;
            else
              lambda_l_tilde = lambda_l_outside;

            Gg = rho_l*lambda_l_tilde*v*k_abs;
          }

        // second equation

        RF fluxl=scale_l * nu_l*lambda_l_tilde/(rho_l*lambda_l_tilde + rho_g * lambda_g_tilde)*(qt + Gl);
        RF fluxg=scale_g * nu_g*lambda_g_tilde/(rho_l*lambda_l_tilde + rho_g * lambda_g_tilde)*(qt + Gg);

        //RF fluxl = scale_l * nu_l * sigma_l * w_l;
        //RF fluxg = scale_g * nu_g * sigma_g * w_g;


        //  // residual accumulation
        r_s.accumulate(lfsu_s,transport,fluxl * face_volume);
        r_n.accumulate(lfsu_n,transport,-fluxl * face_volume);
        r_s.accumulate(lfsu_s,flux, (fluxl+fluxg)* face_volume);
        r_n.accumulate(lfsu_n,flux, -(fluxl+fluxg)* face_volume);
        //r_s.accumulate(lfsu_s,flux,fluxg * face_volume);
        //r_n.accumulate(lfsu_n,flux,-fluxg * face_volume);

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


        //  RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid)+x_s(lfsv_s,capillary));

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);


        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid));
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            r_s.accumulate(lfsu_s,transport,scale_l * nu_l * sigma_l * w_l * face_volume);
            r_s.accumulate(lfsu_s,flux,scale_l * nu_l * sigma_l * w_l * face_volume);
          }


        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF s_g = 1-s_l;

            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF sigma_g = lambda_g_inside*k_abs_inside;
            RF nu_g =tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            r_s.accumulate(lfsu_s,flux,scale_g * nu_g * sigma_g
                           * w_g * face_volume);

          }

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF nu_l = tp.nu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,transport,scale_l * nu_l * j_l * face_volume);
          }


        // gas phase Neumann boundary
        if (bc_g==0)
          {
            if (bc_l!=0)
              DUNE_THROW(Dune::Exception, "mixed boundary conditions (dirichlet for liquid phase and neumann for gas phase NOT implemented");
            RF nu_g =tp.nu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF j_g = tp.j_g(ig.intersection(),face_local,time);

            r_s.accumulate(lfsu_s,flux, scale_g * nu_g * j_g * face_volume);
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
      const TP& tp;  // two phase parameter clas
      Real  scale_l, scale_g;
      typename TP::Traits::RangeFieldType time;

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


    template<class TP>
    class TwoPhaseTotalFluxTemporal
      : public NumericalJacobianVolume<TwoPhaseTotalFluxTemporal<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseTotalFluxTemporal<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };


      TwoPhaseTotalFluxTemporal (const TP& tp_, Real scale_l_, Real scale_g_)
        : tp(tp_), scale_l(scale_l_), scale_g(scale_g_), time(0.)
      {
      }

      //! constructor: pass parameter object
      TwoPhaseTotalFluxTemporal(const TP& tp_)
        : tp(tp_), scale_l(tp.scale_l()), scale_g(tp.scale_g()), time(0.)
      {
      }



      //! set time for subsequent evaluation
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
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
        //RF rho_l = tp.rho_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        //RF rho_g = tp.rho_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));
        RF nu_l = tp.nu_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        RF nu_g = tp.nu_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        r.accumulate(lfsu,flux, phi * ( scale_g * (1. - s_l)  * nu_g *+ scale_l * s_l * nu_l)* cell_volume);
        r.accumulate(lfsu,transport, scale_l * phi * s_l* nu_l * cell_volume);
      }


    private:
      const TP& tp;
      Real scale_l, scale_g;
      typename TP::Traits::RangeFieldType time;
    };








    // a local operator for solving the two-phase flow in liquid pressure - capillary pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    template<typename TP>
    class TwoPhaseLiquidPresentSumSpatial
      : public NumericalJacobianVolume<TwoPhaseLiquidPresentSumSpatial<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseLiquidPresentSumSpatial<TP> >,
        public NumericalJacobianSkeleton<TwoPhaseLiquidPresentSumSpatial<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseLiquidPresentSumSpatial<TP> >,
        public NumericalJacobianBoundary<TwoPhaseLiquidPresentSumSpatial<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseLiquidPresentSumSpatial<TP> >,


        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doAlphaVolume   = true };
      //  enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseLiquidPresentSumSpatial (const TP& tp_)
        : NumericalJacobianSkeleton<TwoPhaseLiquidPresentSumSpatial<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseLiquidPresentSumSpatial<TP> >(eps_numdiff), tp(tp_),time(0.)
      {}


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
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
        r.accumulate(lfsv,flux, -(q_l+q_g) * cell_volume);
        r.accumulate(lfsv,transport, - q_l * cell_volume);
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

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_s(lfsu_s,liquid));

        // compute upwind value of liquid saturation
        RF s_l_s = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
        RF s_l_n = tp.s_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_s,capillary));

        RF s_g_s = 1. - s_l_s;
        RF s_g_n = 1. - s_l_n;


        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_s)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_n)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

        // gas phase velocity
        RF qn = sigma_g *( (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid) )/distance + (x_s(lfsu_s,capillary)-x_n(lfsu_n,capillary) )/distance + aavg(rho_g_inside,rho_g_outside)*gn);

        // harmonic average of k_r
        RF clambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF clambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF csigma_l = havg(clambda_l_inside*k_abs_inside,clambda_l_outside*k_abs_outside);


        r_s.accumulate(lfsu_s,flux,  qn *  aavg(rho_g_inside,rho_g_outside) * face_volume);
        r_n.accumulate(lfsu_n,flux, - qn *  aavg(rho_g_inside,rho_g_outside) *face_volume);



        // determines liquid flow direction
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;

        // first equation
        RF pc_upwind, s_l_upwind_s, s_l_upwind_n;

        // upwind capillary pressure on face
        if (w_l>=0)
          // flow from the current element
          pc_upwind = x_s(lfsu_s,capillary);
        else
          // flow to the current element
          pc_upwind = x_n(lfsu_n,capillary);

        // compute upwind value of liquid saturation
        s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
        s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);

        // harmonic average of k_r
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        RF ql = sigma_l * w_l;

        // residual accumulation
        r_s.accumulate(lfsu_s,transport,sigma_l
                       * w_l * aavg(rho_l_inside,rho_l_outside) * face_volume);
        r_n.accumulate(lfsu_n,transport,-sigma_l
                       * w_l * aavg(rho_l_inside,rho_l_outside) * face_volume);
        r_s.accumulate(lfsu_s,flux, ql * aavg(rho_l_inside,rho_l_outside)* face_volume);
        r_n.accumulate(lfsu_n,flux, -ql * aavg(rho_l_inside,rho_l_outside)* face_volume);

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

        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid));
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid)+x_s(lfsv_s,capillary));

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);


        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            r_s.accumulate(lfsu_s,transport,sigma_l
                           * w_l * rho_l_inside * face_volume);
            r_s.accumulate(lfsu_s,flux,sigma_l
                           * w_l * rho_l_inside * face_volume);
          }


        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)+x_s(lfsu_s,capillary)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF s_g = 1. - s_l;
            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid) + x_s(lfsu_s,capillary));
            RF sigma_g = lambda_g_inside*k_abs_inside;

            r_s.accumulate(lfsu_s,flux, sigma_g * w_l *rho_l_inside * face_volume);

          }

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,transport, j_l *  rho_l_inside * face_volume);
            r_s.accumulate(lfsu_s,flux, j_l * rho_l_inside * face_volume);
          }


        // gas phase Neumann boundary
        if (bc_g==0)
          {
            if (bc_l!=0)
              DUNE_THROW(Dune::Exception, "mixed boundary conditions (dirichlet for liquid phase and neumann for gas phase NOT implemented");
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,flux, j_g * rho_g_inside * face_volume);
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
      const TP& tp;  // two phase parameter class
      typename TP::Traits::RangeFieldType time;

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


    template<class TP>
    class TwoPhaseLiquidPresentSumTemporal
      : public NumericalJacobianVolume<TwoPhaseLiquidPresentSumTemporal<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseLiquidPresentSumTemporal<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseLiquidPresentSumTemporal (const TP& tp_)
        : tp(tp_),time(0.)
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
        RF rho_l = tp.rho_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        RF rho_g = tp.rho_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        r.accumulate(lfsu,flux,phi * (s_l*rho_l + (1. - s_l) * rho_g) * cell_volume);
        r.accumulate(lfsu,transport, phi * s_l * rho_l * cell_volume);
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
    };



    // a local operator for solving the two-phase flow in liquid pressure - capillary pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    template<typename TP>
    class TwoPhaseLiquidPresentSpatial
      : public NumericalJacobianVolume<TwoPhaseLiquidPresentSpatial<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseLiquidPresentSpatial<TP> >,
        public NumericalJacobianSkeleton<TwoPhaseLiquidPresentSpatial<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseLiquidPresentSpatial<TP> >,
        public NumericalJacobianBoundary<TwoPhaseLiquidPresentSpatial<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseLiquidPresentSpatial<TP> >,


        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doAlphaVolume   = true };
      enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseLiquidPresentSpatial (const TP& tp_)
        : NumericalJacobianSkeleton<TwoPhaseLiquidPresentSpatial<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseLiquidPresentSpatial<TP> >(eps_numdiff), tp(tp_),time(0.)
      {}


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
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

        // liquid phase calculation
        RF rho_l = tp.rho_l(eg.entity(),cell_center_local,x(lfsu,liquid));

        // gas phase calculation
        RF rho_g = tp.rho_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        RF q_l = tp.q_l(eg.entity(),cell_center_local,time);
        RF q_g = tp.q_g(eg.entity(),cell_center_local,time);

        // contribution from source term
        r.accumulate(lfsv,flux, -(q_l/rho_l+q_g/rho_g) * cell_volume);
        r.accumulate(lfsv,transport, - q_l/rho_l * cell_volume);
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

        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_s(lfsu_s,liquid));

        // compute upwind value of liquid saturation
        RF s_l_s = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
        RF s_l_n = tp.s_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_s,capillary));

        RF s_g_s = 1. - s_l_s;
        RF s_g_n = 1. - s_l_n;


        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_s)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_n)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        RF sigma_g = havg(lambda_g_inside*k_abs_inside,lambda_g_outside*k_abs_outside);

        // gas phase velocity
        RF qn = sigma_g *( (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid) )/distance + (x_s(lfsu_s,capillary)-x_n(lfsu_n,capillary) )/distance + aavg(rho_g_inside,rho_g_outside)*gn);

        // harmonic average of kr_l
        RF clambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF clambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF csigma_l = havg(clambda_l_inside*k_abs_inside,clambda_l_outside*k_abs_outside);

        r_s.accumulate(lfsu_s,flux,  qn * face_volume);
        r_n.accumulate(lfsu_n,flux, - qn * face_volume);

        // determines liquid flow direction
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;

        // first equation
        RF pc_upwind, s_l_upwind_s, s_l_upwind_n;

        // upwind capillary pressure on face
        if (w_l>=0)
          // flow from the current element
          pc_upwind = x_s(lfsu_s,capillary);
        else
          // flow to the current element
          pc_upwind = x_n(lfsu_n,capillary);

        // compute upwind value of liquid saturation
        s_l_upwind_s = tp.s_l(*(ig.inside()),inside_cell_center_local,pc_upwind);
        s_l_upwind_n = tp.s_l(*(ig.outside()),outside_cell_center_local,pc_upwind);

        // harmonic average of k_r
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_upwind_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_upwind_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF sigma_l = havg(lambda_l_inside*k_abs_inside,lambda_l_outside*k_abs_outside);

        RF ql = sigma_l * w_l;

        // residual accumulation
        r_s.accumulate(lfsu_s,transport,sigma_l
                       * w_l * face_volume);
        r_n.accumulate(lfsu_n,transport,-sigma_l
                       * w_l * face_volume);
        r_s.accumulate(lfsu_s,flux, ql * face_volume);
        r_n.accumulate(lfsu_n,flux, -ql * face_volume);
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
        if (bc_l!=1 && bc_g!=1) return; // no Dirichlet boundary conditions

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

        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid));
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid)+x_s(lfsv_s,capillary));

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);


        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            r_s.accumulate(lfsu_s,transport,sigma_l
                           * w_l * face_volume);
            r_s.accumulate(lfsu_s,flux,sigma_l
                           * w_l * face_volume);
          }


        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)+x_s(lfsu_s,capillary)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF s_g = 1. - s_l;
            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid) + x_s(lfsu_s,capillary));
            RF sigma_g = lambda_g_inside*k_abs_inside;

            r_s.accumulate(lfsu_s,flux, sigma_g * w_l * face_volume);

          }

      }

      // boundary integral independent of ansatz functions
      template<typename IG, typename LFSV, typename R>
      void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r_s) const
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
        RF face_volume = ig.geometry().integrationElement(face_local)*
          Dune::ReferenceElements<DF,dim-1>::general(ig.geometry().type()).volume();

        // evaluate boundary condition type
        int bc_l = tp.bc_l(ig.intersection(),face_local,time);
        int bc_g = tp.bc_g(ig.intersection(),face_local,time);
        if (bc_l!=0 && bc_g!=0) return; // no Neumann boundary conditions

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          inside_cell_center_local = Dune::ReferenceElements<DF,dim>::general(ig.inside()->type()).position(0,0);

        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,0.0);

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsv,transport, j_l * face_volume);
            r_s.accumulate(lfsv,flux, j_l * face_volume);
          }


        // gas phase Neumann boundary
        if (bc_g==0)
          {
            if (bc_l!=0)
              DUNE_THROW(Dune::Exception, "mixed boundary conditions (dirichlet for liquid phase and neumann for gas phase NOT implemented");
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsv,flux, j_g * face_volume);
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
      const TP& tp;  // two phase parameter class
      typename TP::Traits::RangeFieldType time;

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


    template<class TP>
    class TwoPhaseLiquidPresentTemporal
      : public NumericalJacobianVolume<TwoPhaseLiquidPresentTemporal<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseLiquidPresentTemporal<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseLiquidPresentTemporal (const TP& tp_)
        : tp(tp_),time(0.)
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

        r.accumulate(lfsu,transport,phi * s_l * cell_volume);
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
    };




    // a local operator for solving the two-phase flow in liquid pressure - capillary pressure formulation
    // with two-point flux approximation
    // TP : parameter class, see above
    template<typename TP>
    class TwoPhaseTotalFluxOriginal
      : public NumericalJacobianVolume<TwoPhaseTotalFluxOriginal<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseTotalFluxOriginal<TP> >,
        public NumericalJacobianSkeleton<TwoPhaseTotalFluxOriginal<TP> >,
        public NumericalJacobianApplySkeleton<TwoPhaseTotalFluxOriginal<TP> >,
        public NumericalJacobianBoundary<TwoPhaseTotalFluxOriginal<TP> >,
        public NumericalJacobianApplyBoundary<TwoPhaseTotalFluxOriginal<TP> >,


        public FullSkeletonPattern,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,

        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doAlphaVolume   = true };
      //  enum { doLambdaBoundary = true };

      //! constructor: pass parameter object
      TwoPhaseTotalFluxOriginal (const TP& tp_)
        : NumericalJacobianSkeleton<TwoPhaseTotalFluxOriginal<TP> >(eps_numdiff), NumericalJacobianBoundary<TwoPhaseTotalFluxOriginal<TP> >(eps_numdiff), tp(tp_),time(0.)
      {}


      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
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
        r.accumulate(lfsv,flux, -(q_l+q_g) * cell_volume);
        r.accumulate(lfsv,transport, - q_l * cell_volume);
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
        RF k_abs = aavg(k_abs_inside, k_abs_outside);

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);

        // liquid phase calculation
        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF rho_l_outside = tp.rho_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));
        RF rho_l =  aavg(rho_l_inside,rho_l_outside);
        // gas phase calculation
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF rho_g_outside = tp.rho_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_s(lfsu_s,liquid));
        RF rho_g = aavg(rho_g_inside,rho_g_outside);

        // compute upwind value of liquid saturation
        RF s_l_s = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
        RF s_l_n = tp.s_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_s,capillary));

        RF s_g_s = 1. - s_l_s;
        RF s_g_n = 1. - s_l_n;


        RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g_s)/
          tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
        RF lambda_g_outside = tp.kr_g(*(ig.outside()),outside_cell_center_local,s_g_n)/
          tp.mu_g(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid));
        RF lambda_g =  havg(lambda_g_inside,lambda_g_outside);


        // harmonic average of k_r
        RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l_s)/
          tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
        RF lambda_l_outside = tp.kr_l(*(ig.outside()),outside_cell_center_local,s_l_n)/
          tp.mu_l(*(ig.outside()),outside_cell_center_local,x_n(lfsu_n,liquid));

        //  - grad p_l + rho_l * g * n,e = - (grad p_l - rho_l * g * n,e)
        RF w_l = (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid))/distance + aavg(rho_l_inside,rho_l_outside)*gn;
        //  - grad p_l - grad p_c + rho_l * g * n,e = - (grad p_l + grad p_c - rho_g * g * n,e)
        RF w_g = ((x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid))-(x_n(lfsu_n,capillary)+x_n(lfsu_n,liquid)))/distance + aavg(rho_g_inside,rho_g_outside)*gn;



        RF lambda_l =  havg(lambda_l_inside,lambda_l_outside);

        RF sum = rho_l*lambda_l + rho_g*lambda_g;

        RF G = (rho_l*rho_l*lambda_l + rho_g*rho_g*lambda_g)/sum;

        RF qt = sum * k_abs *( (x_s(lfsu_s,liquid)-x_n(lfsu_n,liquid) )/distance + rho_g*lambda_g/sum*(x_s(lfsu_s,capillary)-x_n(lfsu_n,capillary) )/distance + G*gn );

        // - (grad p_c + (rho_l- rho_g) * g * n,e)
        RF v = w_g - w_l;

        // Gl = - rho_n * lambda_n * v
        RF Gl = -rho_g * v;
        // Gg = rho_l * lambda_l * v
        RF Gg = rho_l* v;

        RF lambda_l_tilde = 0.0;
        RF lambda_g_tilde = 0.0;
        RF fluxl(0), fluxg(0);

        if ((Gl*qt)>=0) // alpa = liquid
          {
            if (qt>=0)
              lambda_l_tilde = lambda_l_inside;
            else
              lambda_l_tilde = lambda_l_outside;

            Gg = rho_l*lambda_l_tilde*v*k_abs;

            if ((qt+Gg)>=0)
              lambda_g_tilde = lambda_g_inside;
            else
              lambda_g_tilde = lambda_g_outside;

            Gl = -rho_g*lambda_g_tilde*v*k_abs;
          }
        else
          {
            if (qt>=0)
              lambda_g_tilde = lambda_g_inside;
            else
              lambda_g_tilde = lambda_g_outside;

            Gl = - rho_g*lambda_g_tilde*v*k_abs;

            if ((qt+Gl)>=0)
              lambda_l_tilde = lambda_l_inside;
            else
              lambda_l_tilde = lambda_l_outside;

            Gg = rho_l*lambda_l_tilde*v*k_abs;
          }

        fluxl=rho_l*lambda_l_tilde/(rho_l*lambda_l_tilde + rho_g * lambda_g_tilde)*(qt + Gl);
        fluxg=rho_g*lambda_g_tilde/(rho_l*lambda_l_tilde + rho_g * lambda_g_tilde)*(qt + Gg);


        // residual accumulation
        r_s.accumulate(lfsu_s,transport,fluxl * face_volume);
        r_n.accumulate(lfsu_n,transport,-fluxl * face_volume);
        r_s.accumulate(lfsu_s,flux, (fluxg + fluxl)* face_volume);
        r_n.accumulate(lfsu_n,flux, -(fluxg + fluxl)* face_volume);

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

        RF rho_l_inside = tp.rho_l(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid));
        RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsv_s,liquid)+x_s(lfsv_s,capillary));

        // gravity times normal
        RF gn = tp.gravity()*ig.unitOuterNormal(face_local);


        // liquid phase Dirichlet boundary
        if (bc_l==1)
          {
            RF g_l = tp.g_l(ig.intersection(),face_local,time);
            RF w_l = (x_s(lfsu_s,liquid)-g_l)/distance + rho_l_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF lambda_l_inside = tp.kr_l(*(ig.inside()),inside_cell_center_local,s_l)/
              tp.mu_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,liquid));
            RF sigma_l = lambda_l_inside*k_abs_inside;
            r_s.accumulate(lfsu_s,transport,sigma_l
                           * w_l * rho_l_inside * face_volume);
            r_s.accumulate(lfsu_s,flux,sigma_l
                           * w_l * rho_l_inside * face_volume);
          }


        // gas phase Dirichlet boundary
        if (bc_g==1)
          {
            RF rho_g_inside = tp.rho_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF g_g = tp.g_g(ig.intersection(),face_local,time);
            RF w_g = (x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid)-g_g)/distance + rho_g_inside*gn;
            RF s_l = tp.s_l(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary));
            RF s_g = 1-s_l;

            RF lambda_g_inside = tp.kr_g(*(ig.inside()),inside_cell_center_local,s_g)/
              tp.mu_g(*(ig.inside()),inside_cell_center_local,x_s(lfsu_s,capillary)+x_s(lfsu_s,liquid));
            RF sigma_g = lambda_g_inside*k_abs_inside;

            r_s.accumulate(lfsu_s,flux, sigma_g * w_g *rho_g_inside * face_volume);

          }

        // liquid phase Neumann boundary
        if (bc_l==0)
          {
            RF j_l = tp.j_l(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,transport, j_l *  rho_l_inside * face_volume);
            r_s.accumulate(lfsv_s,flux, j_l *  rho_l_inside * face_volume);
          }


        // gas phase Neumann boundary
        if (bc_g==0)
          {
            if (bc_l!=0)
              DUNE_THROW(Dune::Exception, "mixed boundary conditions (dirichlet for liquid phase and neumann for gas phase NOT implemented");
            RF j_g = tp.j_g(ig.intersection(),face_local,time);
            r_s.accumulate(lfsu_s,flux, j_g * rho_g_inside * face_volume);
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
      const TP& tp;  // two phase parameter class
      typename TP::Traits::RangeFieldType time;

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


    template<class TP>
    class TwoPhaseTotalFluxOriginalTemporal
      : public NumericalJacobianVolume<TwoPhaseTotalFluxOriginalTemporal<TP> >,
        public NumericalJacobianApplyVolume<TwoPhaseTotalFluxOriginalTemporal<TP> >,
        public FullVolumePattern,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };
      enum { liquid = TwoPhaseLiquid };
      enum { capillary = TwoPhaseCapillary };

      enum { flux = TwoPhaseLiquid };
      enum { transport = TwoPhaseCapillary };

      typedef typename TP::Traits::RangeFieldType Real;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TwoPhaseTotalFluxOriginalTemporal (const TP& tp_)
        : tp(tp_),time(0.)
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
        RF rho_l = tp.rho_l(eg.entity(),cell_center_local,x(lfsu,liquid));
        RF rho_g = tp.rho_g(eg.entity(),cell_center_local,x(lfsu,capillary)+x(lfsu,liquid));

        r.accumulate(lfsu,flux,phi * ( (1. - s_l) * rho_g + s_l * rho_l)* cell_volume);
        r.accumulate(lfsu,transport, phi * s_l * rho_l * cell_volume);
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
    };



  }
}

#endif
