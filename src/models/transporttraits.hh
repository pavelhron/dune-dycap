// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_TRANSPORTTRAITS_HH
#define DUNE_PDELAB_TRANSPORTTRAITS_HH

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>

namespace Dune {
  namespace PDELab {

    // reaction base class
    // do not compute anything
    class ReactionBaseAdapter
    {
    public:
      //! constructor stores reference to the model
      ReactionBaseAdapter()
      {}

      //! do one step
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void evaluate (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
      {

      }

      //! set time step for subsequent steps
      template<typename RF>
      void setTime (RF time_)
      {
      }
    };


    //! traits class for two phase parameter class
    template<typename GV, typename RF>
    struct TransportParameterTraits
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
      typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> PermTensorType;

      //! grid types
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Intersection IntersectionType;

      /** \brief Class to define the boundary condition types
       */
      struct ConvectionDiffusionBoundaryConditions
      {
        enum Type { Dirichlet=1, Neumann=-1, Outflow=-2, None=-3, Robin=2 }; // BC requiring constraints must be >0 if
        // constraints assembler coming with PDELab is used
      };

      typedef typename ConvectionDiffusionBoundaryConditions::Type BCType;
    };

    //! base class for parameter class
    template<class T, class Imp>
    class TransportSpatialParameterInterface
    {
    public:
      typedef T Traits;
      typedef typename Traits::BCType BCType;

      //! velocity field
      typename Traits::RangeType
      v (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().v(e,x);
      }

      //! scalar diffusion coefficient
      typename Traits::RangeFieldType
      D (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().D(e,x);
      }

      //! source term
      typename Traits::RangeFieldType
      q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().q(e,x);
      }

      //! boundary condition type function
      /**
       * 0 means Neumann
       * 1 means Dirichlet
       * 2 means Outflow (zero diffusive flux, velocity points outward)
       */
      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return asImp().bctype(is,x);
      }

      //! Dirichlet boundary condition on inflow
      typename Traits::RangeFieldType
      g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return asImp().g(is,x);
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j  (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
      {
        return asImp().j(is,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };




    //! base class for parameter class
    template<class T, class Imp>
    class ModifiedTransportSpatialParameterInterface
      : public TransportSpatialParameterInterface<T,Imp>
    {
    public:
      typedef T Traits;
      typedef typename Traits::BCType BCType;

      //! saturation at end of big step
      typename Traits::RangeFieldType
      snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().snew(e,x);
      }

      //! saturation at end of big step
      typename Traits::RangeFieldType
      sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().sold(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };



    //! base class for parameter class
    template<class T, class Imp>
    class TransportTemporalParameterInterface
    {
    public:
      typedef T Traits;
      typedef typename Traits::BCType BCType;

      //! source term
      typename Traits::RangeFieldType
      c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().c(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    //! base class for parameter class
    template<class T, class Imp>
    class ModifiedTransportTemporalParameterInterface
      : public TransportTemporalParameterInterface<T,Imp>
    {
    public:
      typedef T Traits;
      typedef typename Traits::BCType BCType;

      //! saturation at end of big step
      typename Traits::RangeFieldType
      snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
      {
        return asImp().snew(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };



    //! base class for parameter class
    template<class T, class Imp>
    class TransportMulticomponentInterface
    {
    public:
      typedef T Traits;
      typedef typename Traits::BCType BCType;

      //! velocity field
      typename Traits::RangeType
      v (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).v(e,x);
      }

      //! scalar diffusion coefficient
      typename Traits::RangeFieldType
      D (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).D(e,x);
      }

      //! source term
      typename Traits::RangeFieldType
      q (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).q(e,x);
      }

      BCType
      bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i)
      {
        return asImp().component(i).bctype(is,x);
      }

      //! Dirichlet boundary condition on inflow
      typename Traits::RangeFieldType
      g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i) const
      {
        return asImp().component(i).g(is,x);
      }

      //! Neumann boundary condition
      typename Traits::RangeFieldType
      j  (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, std::size_t i) const
      {
        return asImp().component(i).j(is,x);
      }

      //! saturation at end of big step
      typename Traits::RangeFieldType
      snew (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).snew(e,x);
      }

      //! saturation at end of big step
      typename Traits::RangeFieldType
      sold (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).sold(e,x);
      }

      void setTime(typename Traits::RangeFieldType t)
      {
        asImp().setTime(t);
      }

      void preStep(typename Traits::RangeFieldType t)
      {
        asImp().preStep(t);
      }

      void setTimeTarget(typename Traits::RangeFieldType t)
      {
        asImp().setTimeTarget(t);
      }



    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


       //! base class for parameter class
    template<class T, class Imp>
    class TransportMulticomponentTemporalInterface
    {
    public:
      typedef T Traits;

      //! source term
      typename Traits::RangeFieldType
      c (const typename Traits::ElementType& e, const typename Traits::DomainType& x, std::size_t i) const
      {
        return asImp().component(i).c(e,x);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


  }
}

#endif
