// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TWOPHASETRAITS_HH
#define DUNE_PDELAB_TWOPHASETRAITS_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>

namespace Dune {
  namespace PDELab {
    const double eps_numdiff = 1.e-10;

    //! /brief discontinuity is replaced with smooth function
    template<class T>
    T Blend (T x, T x0, T x1)
    {
      if (x>=x1) return 1.0;
      if (x<=x0) return 0.0;
      x = (x-x0)/(x1-x0);
      if (x<=0.5)
        return(2.*x*x);
      else
        return(1.-2.*(1.-x)*(1.-x));
    }

    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct TwoPhaseParameterTraits
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

      //! \brief permeability tensor type
      typedef RangeFieldType PermTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;
    };

    template<typename GV, typename RF>
    struct TwoPhaseFullTensorParameterTraits : TwoPhaseParameterTraits<GV, RF>
    {
      typedef TwoPhaseParameterTraits<GV, RF> Base;
      typedef typename Base::RangeFieldType RangeFieldType;

      //! \brief permeability tensor type
      typedef Dune::FieldMatrix<RangeFieldType,Base::dimDomain,Base::dimDomain> PermTensorType;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class TwoPhaseParameterInterface
    {
    public:
      typedef T Traits;

      //! porosity
      typename Traits::RangeFieldType
      phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().phi(e,x);
      }

      //! entry pressure to detect media discontinuities
      typename Traits::RangeFieldType
      pe (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().pe(e,x);
      }

      //! capillary pressure function
      typename Traits::RangeFieldType
      pc (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
          typename Traits::RangeFieldType s_l) const
      {
        return asImp().pc(e,x,s_l);
      }

      //! inverse capillary pressure function
      typename Traits::RangeFieldType
      s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType pc) const
      {
        return asImp().s_l(e,x,pc);
      }

      //! liquid phase relative permeability
      typename Traits::RangeFieldType
      kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType s_l) const
      {
        return asImp().kr_l(e,x,s_l);
      }

      //! gas phase relative permeability
      typename Traits::RangeFieldType
      kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType s_g) const
      {
        return asImp().kr_g(e,x,s_g);
      }

      //! liquid phase viscosity
      typename Traits::RangeFieldType
      mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_l) const
      {
        return asImp().mu_l(e,x,p_l);
      }

      //! gas phase viscosity
      typename Traits::RangeFieldType
      mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_g) const
      {
        return asImp().mu_g(e,x,p_g);
      }

      //! absolute permeability (scalar!)
      typename Traits::PermTensorType
      k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().k_abs(e,x);
      }

      //! gravity vector
      const typename Traits::RangeType& gravity () const
      {
        return asImp().gravity();
      }

      //! liquid phase molar density
      template<typename E>
      typename Traits::RangeFieldType
      nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_l) const
      {
        return asImp().nu_l(e,x,p_l);
      }

      //! gas phase molar density
      typename Traits::RangeFieldType
      nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
            typename Traits::RangeFieldType p_g) const
      {
        return asImp().nu_g(e,x,p_g);
      }

      //! liquid phase mass density
      typename Traits::RangeFieldType
      rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
             typename Traits::RangeFieldType p_l) const
      {
        return asImp().rho_l(e,x,p_l);
      }

      //! gas phase mass density
      typename Traits::RangeFieldType
      rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
             typename Traits::RangeFieldType p_g) const
      {
        return asImp().rho_g(e,x,p_g);
      }

      //! liquid phase boundary condition type
      int
      bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().bc_l(is,x,time);
      }

      //! gas phase boundary condition type
      int
      bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().bc_g(is,x,time);
      }

      //! liquid phase Dirichlet boundary condition
      typename Traits::RangeFieldType
      g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().g_l(is,x,time);
      }

      //! gas phase Dirichlet boundary condition
      typename Traits::RangeFieldType
      g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().g_g(is,x,time);
      }

      //! liquid phase Neumann boundary condition
      typename Traits::RangeFieldType
      j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().j_l(is,x,time);
      }

      //! gas phase Neumann boundary condition
      typename Traits::RangeFieldType
      j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
      {
        return asImp().j_g(is,x,time);
      }

      //! liquid phase source term
      typename Traits::RangeFieldType
      q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType time) const
      {
        return asImp().q_l(e,x,time);
      }

      //! gas phase source term
      typename Traits::RangeFieldType
      q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
           typename Traits::RangeFieldType time) const
      {
        return asImp().q_g(e,x,time);
      }

      //! scale function for liquid phase
      typename Traits::RangeFieldType
      scale_l () const
      {
        return 1.0;
      }

      //! scale function for gas phase
      typename Traits::RangeFieldType
      scale_g () const
      {
        return 1.0;
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

  }
}

#endif
