
// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PM_REACTIONOP_HH
#define DUNE_PM_REACTIONOP_HH
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include"defaultimp.hh"
#include<dune/pdelab/localoperator/idefault.hh>



namespace Dune {
  namespace PDELab {

    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct ReactionParameterTraits
    {
      //! \brief the grid view
      typedef GV GridViewType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = GV::dimension
      };

      //! \brief Export type for domain field
      typedef typename GV::Grid::ctype DomainFieldType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

      //! \brief domain type
      typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief range type
      typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class ReactionParameterInterface
    {
    public:
      typedef T Traits;

      //! porosity
      typename Traits::RangeFieldType
      phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().phi(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      sl(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().sl(e,x);
      }


      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Mumax (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Mumax(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Ks (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Ks(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Ko (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Ko(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Rd (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Rd(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Ys (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Ys(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Yo (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Yo(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Exp (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Exp(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Mumax_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Mumax_an(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Ks_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Ks_an(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Rd_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Rd_an(e,x);
      }

      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Ys_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Ys_an(e,x);
      }


      //! maximum specific growth rate
      typename Traits::RangeFieldType
      Exp_an (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().Exp_an(e,x);
      }

      //! gas-liquid mass transfer coefficient * volumetric gas-liquid interfacial area (33.4 h^-1)
      typename Traits::RangeFieldType
      k_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().k_l(e,x);
      }

      //! Henry's Constant
      typename Traits::RangeFieldType
      K_H(const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().K_H(e,x);
      }

      //! oxygen concentration in water
      typename Traits::RangeFieldType
      c_O2_liquid (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
                   typename Traits::RangeFieldType c) const
      {
        return asImp().c_O2_liquid(e,x,c);
      }

      //! oxygen concentration in air
      typename Traits::RangeFieldType
      c_O2_gas (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
                typename Traits::RangeFieldType c) const
      {
        return asImp().c_O2_gas(e,x,c);
      }

      //! DOC concentration
      typename Traits::RangeFieldType
      c_DOC (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
             typename Traits::RangeFieldType c) const
      {
        return asImp().c_DOC(e,x,c);
      }

      //! Ecoli concentration
      typename Traits::RangeFieldType
      c_Ecoli (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
               typename Traits::RangeFieldType c) const
      {
        return asImp().c_Ecoli(e,x,c);
      }

      //! set time for subsequent evaluation
      void setTime (typename Traits::RangeFieldType t)
      {
        return asImp().setTime(t);
      }

      //! to be called once before each time step
      void preStep (typename Traits::RangeFieldType time_, typename Traits::RangeFieldType dt_, int stages)
      {
        return asImp().preStep(time_,dt_,stages);
      }


    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    /*! \brief a local operator for a Dopel Monod Model, oxygen dissolution is fast enough \n

      \f$C_{0_2} \dots \quad\f$ dissolved oxygen concentration \n
      \f$C_{0_2,g} \dots \quad\f$ air oxygen concentration \n
      \f$C_{0_2}^* \dots \quad\f$ dissolved oxygen concentration at equilibrium    \n
      \f$k_l \, \alpha \dots \quad\f$  gas/liquid mass-transfer coefficient, oxygen mass transfer coefficient based on liquid volume and concetration driving force  \n
      \f$K_o\dots \quad\f$ oxygen half saturation constant  \n
      \f$K_s \dots \quad\f$  DOC half saturation constant \n
      \f$K_{s,ae} \dots \quad\f$  DOC half saturation constant, anaerobic growth \n
      \f$ C_{DOC} \dots \quad\f$  organic carbon (DOC) concentration  \n
      \f$X \dots \quad \f$  cell mass concentration \n
      \f$Y_{x,o} \dots \quad\f$ growth yield coefficient based on oxygen utilized  \n
      \f$Y_{x,s} \dots \quad\f$ growth yield coefficient based on DOC utilized  \n
      \f$\mu_{max} \dots \quad\f$ maximum specific growth rate  \n
      \f$\mu_{max,ae} \dots \quad\f$ maximum specific growth rate, aerobic growth  \n
      \f$exp \dots \dots \quad\f$ exponent  \n

      Let denote:
      \f{equation*}{
      \lambda = \left( \mu_{max}  \left(\frac{C_{DOC}}{K_s + C_{DOC}}\right)^{exp} \frac{C_{0_2}}{K_o + C_{0_2}} +  \mu_{max,ae}  \frac{C_{DOC}}{K_{s,ae} + C_{DOC}} \right)\, X,
      \f}

      Dissolved oxygen concentration:
      \f{equation*}{
      \frac{d \, C_{0_2} + C_{0_2,g}}{d \, t} =  -\frac{\lambda}{Y_{x,o}}
      \f}

      Equilibrium concentration \f$C_{0_2}^*\f$ for gases of low solubility is supplied by Henry's law in the form
      \f{equation*}{
      C_{0_2}^*= k_h \, C_{0_2,g} \,R \,T = k_H \, C_{0_2,g},
      \f}
      where \f$k_H \f$ is a Henry's law constant.

      The cell growth:
      \f{equation*}{
      \frac{d \, X}{d \, t} = \lambda,
      \f}

      DOC concentration change:
      \f{equation*}{
      \frac{d \, C_{DOC}}{d \, t} = -\frac{\lambda}{Y_{x,s}},
      \f}

    */


    template<typename TP>
    class SpatialDMMStaticOperator:
      public NumericalJacobianVolume<SpatialDMMStaticOperator<TP> >,
      public NumericalJacobianApplyVolume<SpatialDMMStaticOperator<TP> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      enum { oliquid = 0};
      enum { ogas = 1};
      enum { doc = 2};
      enum { ecoli = 3};

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };

      SpatialDMMStaticOperator (TP& tp_, typename TP::Traits::RangeFieldType z=1e-7) DUNE_DEPRECATED
        : tp(tp_), zero(z)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {

        // select the two components
        typedef typename LFSV::template Child<oliquid>::Type O2LiquidSpace;

        // domain and range field type
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);


        RF mumax = tp.Mumax(eg.entity(),cell_center_local);
        RF mumax_ae = tp.Mumax_an(eg.entity(),cell_center_local);
        RF Ks =  tp.Ks(eg.entity(),cell_center_local);
        RF Ks_ae =  tp.Ks_an(eg.entity(),cell_center_local);
        RF Ko =  tp.Ko(eg.entity(),cell_center_local);
        RF Kd =  tp.Rd(eg.entity(),cell_center_local);
        RF Ys =  tp.Ys(eg.entity(),cell_center_local);
        RF Yo =  tp.Yo(eg.entity(),cell_center_local);
        //  RF exponent =  tp.exponent(eg.entity(),cell_center_local);

        RF c_O2_liquid = tp.c_O2_liquid(eg.entity(),cell_center_local,x(lfsu,oliquid));
        RF c_DOC = tp.c_DOC(eg.entity(),cell_center_local,x(lfsu,doc));
        RF X = tp.c_Ecoli(eg.entity(),cell_center_local,x(lfsu,ecoli));
        RF K_H =  tp.K_H(eg.entity(),cell_center_local);

        RF ft1 = mumax * c_DOC/ (c_DOC + Ks) * c_O2_liquid/ (c_O2_liquid + Ko)*X;
        RF ft2 = mumax_ae * c_DOC/ (c_DOC + Ks_ae)*X;

        // contribution from source term
        r.accumulate(lfsv,oliquid,ft1/Yo);
        r.accumulate(lfsv,ogas,x(lfsu,ogas)*K_H-x(lfsu,oliquid));
        r.accumulate(lfsv,doc,(ft1+ft2)/Ys);
        r.accumulate(lfsv,ecoli,-(ft1+ft2)+Kd*X);
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
      }


      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
        dtmin = 1E100;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return dtmin;
      }

    private:

      TP& tp;
      typename TP::Traits::RangeFieldType time;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      typename TP::Traits::RangeFieldType zero;
    };




    /** a local operator for the storage operator
     */
    template<class TP>
    class TemporalDMMStaticOperator
      :   public NumericalJacobianVolume<TemporalDMMStaticOperator<TP> >,
          public NumericalJacobianApplyVolume<TemporalDMMStaticOperator<TP> >,
          public FullVolumePattern,
          public LocalOperatorDefaultFlags,
          public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      enum { oliquid = 0};
      enum { ogas = 1};
      enum { doc = 2};
      enum { ecoli = 3};

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TemporalDMMStaticOperator (TP& tp_)
        : tp(tp_), zero(1e-7)
      {
      }

      TemporalDMMStaticOperator (TP& tp_, typename TP::Traits::RangeFieldType z) DUNE_DEPRECATED
        : tp(tp_), zero(z)
      {
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<oliquid>::Type O2LiquidSpace;

        // domain and range field type
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        //typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);


        RF c_O2_liquid = tp.c_O2_liquid(eg.entity(),cell_center_local,x(lfsu,oliquid));
        RF c_O2_gas = tp.c_O2_gas(eg.entity(),cell_center_local,x(lfsu,ogas));
        RF c_DOC = tp.c_DOC(eg.entity(),cell_center_local,x(lfsu,doc));
        RF X = tp.c_Ecoli(eg.entity(),cell_center_local,x(lfsu,ecoli));

        // contribution from source term
        r.accumulate(lfsv,oliquid,c_O2_liquid+c_O2_gas); //?? 2*c_O2_liquid
        r.accumulate(lfsv,ogas,0.);
        r.accumulate(lfsv,doc,c_DOC);
        r.accumulate(lfsv,ecoli,X);
      }

      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return 1e100;
      }

    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
      typename TP::Traits::RangeFieldType zero;
    };


    /*! \brief a local operator for a Dopel Monod Model, dynamic oxygen dissolution \n

      \f$C_{0_2} \dots \quad\f$ dissolved oxygen concentration \n
      \f$C_{0_2,g} \dots \quad\f$ air oxygen concentration \n
      \f$C_{0_2}^* \dots \quad\f$ dissolved oxygen concentration at equilibrium    \n
      \f$k_l \, \alpha \dots \quad\f$  gas/liquid mass-transfer coefficient, oxygen mass transfer coefficient based on liquid volume and concetration driving force  \n
      \f$K_o\dots \quad\f$ oxygen half saturation constant  \n
      \f$K_s \dots \quad\f$  DOC half saturation constant \n
      \f$ C_{DOC} \dots \quad\f$  organic carbon (DOC) concentration  \n
      \f$X \dots \quad \f$  cell mass concentration \n
      \f$Y_{x,o} \dots \quad\f$ growth yield coefficient based on oxygen utilized  \n
      \f$Y_{x,s} \dots \quad\f$ growth yield coefficient based on DOC utilized  \n
      \f$\mu_{max} \dots \quad\f$ maximum specific growth rate  \n
      \f$exp \dots \dots \quad\f$ exponent  \n

      Equilibrium concentration \f$C_{0_2}^*\f$ for gases of low solubility is supplied by Henry's law in the form
      \f{equation*}{
      C_{0_2}^*= k_h \, C_{0_2,g} \,R \,T = k_H \, C_{0_2,g},
      \f}
      where \f$k_H \f$ is a Henry's law constant.

      Dissolved oxygen concentration:
      \f{equation*}{
      \frac{d \, C_{0_2}}{d \, t} = k_l \left(C_{0_2}^* -C_{0_2}\right)  -\frac{\mu_{max}}{Y_{x,o}}  \frac{C_{DOC}}{K_s + C_{DOC}} \frac{C_{0_2}}{K_o + C_{0_2}}X.
      \f}

      Oxygen transfer from gas phase
      \f{equation*}{
      \frac{d \, C_{0_2,g}}{d \, t} = - k_l \left(C_{0_2}^* -C_{0_2}\right) \frac{s_l}{s_g}.
      \f}


      The cell growth:
      \f{equation*}{
      \frac{d \, X}{d \, t} = \mu_{max}  \frac{C_{DOC}}{K_s + C_{DOC}} \frac{C_{0_2}}{K_o + C_{0_2}}X,
      \f}

      Organic carbon concentration change:
      \f{equation*}{
      \frac{d \, C_{DOC}}{d \, t} = -\frac{\mu_{max}}{Y_{x,s}}  \frac{C_{DOC}}{K_s + C_{DOC}} \frac{C_{0_2}}{K_o + C_{0_2}}X,
      \f}

      Time step restriction:
      \f{equation*}{
      k_l k_H \, C_{0_2,g}  -\frac{\mu_{max}}{Y_{x,o}}  \frac{C_{DOC}}{K_s + C_{DOC}} \frac{1}{K_o}X< \frac{1}{dt}
      \f}

      \f{equation*}{
      k_l k_H \frac{s_l}{s_g}< \frac{1}{dt}
      \f}

      \f{equation*}{
      \mu_{max}  \frac{C_{DOC}}{K_s + C_{DOC}} \frac{C_{0_2}}{K_o + C_{0_2}}<\frac{1}{dt}
      \f}

      \f{equation*}{
      \frac{\mu_{max}}{Y_{x,s}}  \frac{1}{K_s} \frac{C_{0_2}}{K_o + C_{0_2}}X <\frac{1}{dt}
      \f}


    */


    template<typename TP>
    class SpatialDMMDynamicOperator:
      public NumericalJacobianVolume<SpatialDMMDynamicOperator<TP> >,
      public NumericalJacobianApplyVolume<SpatialDMMDynamicOperator<TP> >,
    // public NumericalJacobianVolumePostSkeleton<SpatialDMMDynamicOperator<TP> >,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      enum { oliquid = 0};
      enum { ogas = 1};
      enum { doc = 2};
      enum { ecoli = 3};

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaVolumePostSkeleton = true };

      SpatialDMMDynamicOperator (TP& tp_, typename TP::Traits::RangeFieldType z=1e-7) DUNE_DEPRECATED
        : tp(tp_), zero(z)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {

        // select the two components
        typedef typename LFSV::template Child<oliquid>::Type O2LiquidSpace;

        // domain and range field type
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //  Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

        celldt = 0.0;

        RF mumax = tp.mumax(eg.entity(),cell_center_local);
        RF k_l = tp.k_l(eg.entity(),cell_center_local);
        RF Ks =  tp.Ks(eg.entity(),cell_center_local);
        RF Ko =  tp.Ko(eg.entity(),cell_center_local);
        RF Ys =  tp.Ys(eg.entity(),cell_center_local);
        RF Yo =  tp.Yo(eg.entity(),cell_center_local);
        RF exponent =  tp.exponent(eg.entity(),cell_center_local);

        RF sl = tp.sl(eg.entity(),cell_center_local);

        RF c_O2_liquid = tp.c_O2_liquid(eg.entity(),cell_center_local,x(lfsu,oliquid));
        RF c_DOC = tp.c_DOC(eg.entity(),cell_center_local,x(lfsu,doc));
        RF X = tp.c_Ecoli(eg.entity(),cell_center_local,x(lfsu,ecoli));
        RF K_H =  tp.K_H(eg.entity(),cell_center_local);

        RF ft = mumax * std::pow(c_DOC/ (c_DOC + Ks),exponent) * c_O2_liquid/ (c_O2_liquid + Ko)*X;


        celldt = std::abs(X/ft);
        celldt = std::min(celldt,std::abs(c_DOC*Ys/ft));
        celldt = std::min(celldt,std::abs(c_O2_liquid / (k_l*(K_H*x(lfsu,ogas)-x(lfsu,oliquid))+ft/Yo)));

        //  std::cout << celldt << std::endl;
        // contribution from source term
        r.accumulate(lfsv,oliquid,-k_l*(K_H*x(lfsu,ogas)-x(lfsu,oliquid))+ft/Yo);
        r.accumulate(lfsv,ogas,-k_l*(K_H*x(lfsu,ogas)-x(lfsu,oliquid))*sl/(1.-sl));
        r.accumulate(lfsv,doc,ft/Ys);
        r.accumulate(lfsv,ecoli,-ft);

      }

      // post skeleton: compute time step allowable for cell; to be done later
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                      const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<oliquid>::Type O2LiquidSpace;

        // domain and range field type
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        //  typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;

        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);

        dtmin = std::min(dtmin,celldt);
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
      }


      //! to be called once before each stage
      void preStage (typename TP::Traits::RangeFieldType time, int r)
      {
        dtmin = 1E100;
      }

      //! to be called once at the end of each stage
      void postStage ()
      {
      }

      //! to asked after first stage
      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        //std::cout << "dtmin " << dtmin << std::endl;
        return dtmin*10;
      }

    private:

      TP& tp;
      typename TP::Traits::RangeFieldType time;
      mutable typename TP::Traits::RangeFieldType dtmin; // accumulate minimum dt here
      mutable typename TP::Traits::RangeFieldType celldt;
      typename TP::Traits::RangeFieldType zero;
    };




    /** a local operator for the storage operator
     */
    template<class TP>
    class TemporalDMMDynamicOperator
      :   public NumericalJacobianVolume<TemporalDMMDynamicOperator<TP> >,
          public NumericalJacobianApplyVolume<TemporalDMMDynamicOperator<TP> >,
          public FullVolumePattern,
          public LocalOperatorDefaultFlags,
          public InstationaryLocalOperatorDefaultMethods<typename TP::Traits::RangeFieldType>
    {
      enum { dim = TP::Traits::GridViewType::dimension };

      enum { oliquid = 0};
      enum { ogas = 1};
      enum { doc = 2};
      enum { ecoli = 3};

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      TemporalDMMDynamicOperator (TP& tp_)
        : tp(tp_), zero(1e-7)
      {
      }

      TemporalDMMDynamicOperator (TP& tp_, typename TP::Traits::RangeFieldType z)
        : tp(tp_), zero(z)
      {
      }

      //! set time in parameter class
      void setTime (typename TP::Traits::RangeFieldType t)
      {
        time = t;
        tp.setTime(t);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // select the two components
        typedef typename LFSV::template Child<oliquid>::Type O2LiquidSpace;

        // domain and range field type
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename O2LiquidSpace::Traits::FiniteElementType::
          Traits::LocalBasisType::Traits::RangeFieldType RF;
        // typedef typename O2LiquidSpace::Traits::FiniteElementType::
        //   Traits::LocalBasisType::Traits::RangeType RangeType;

        // dimensions
        const int dim = EG::Geometry::dimension;
        // cell geometry
        const Dune::FieldVector<DF,dim>&
          cell_center_local = Dune::ReferenceElements<DF,dim>::general(eg.geometry().type()).position(0,0);


        RF c_O2_liquid = tp.c_O2_liquid(eg.entity(),cell_center_local,x(lfsu,oliquid));
        RF c_O2_gas = tp.c_O2_gas(eg.entity(),cell_center_local,x(lfsu,ogas));
        RF c_DOC = tp.c_DOC(eg.entity(),cell_center_local,x(lfsu,doc));
        RF X = tp.c_Ecoli(eg.entity(),cell_center_local,x(lfsu,ecoli));

        // contribution from source term
        r.accumulate(lfsv,oliquid,c_O2_liquid);
        r.accumulate(lfsv,ogas,c_O2_gas);
        r.accumulate(lfsv,doc,c_DOC);
        r.accumulate(lfsv,ecoli,X);
      }

      typename TP::Traits::RangeFieldType suggestTimestep (typename TP::Traits::RangeFieldType dt) const
      {
        return 1e100;
      }

    private:
      TP& tp;
      typename TP::Traits::RangeFieldType time;
      typename TP::Traits::RangeFieldType zero;
    };


  }
}

#endif
