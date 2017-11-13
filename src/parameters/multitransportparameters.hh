// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_MULTITRANSPORTPARAMETERS_HH
#define DUNE_DYCAP_MULTITRANSPORTPARAMETERS_HH

#include <assert.h>
#include <iostream>
#include <type_traits>
#include <memory>
#include <cstdint>
#include <dune/common/array.hh>
#include"../models/multicomponenttransportop.hh"
#include"transportparameters.hh"

namespace Dune {
  namespace PDELab {


    namespace {

      // prototype and end of recursion
      template<typename T, typename It, typename... Args>
      void assign_reference_pack_to_shared_ptr_array_unpack(It it, Args&&... args) {}

      template<typename T, typename It, typename Arg, typename... Args>
      void assign_reference_pack_to_shared_ptr_array_unpack(It it, Arg&& arg, Args&&... args)
      {
        static_assert(std::is_same<T,typename std::remove_const<typename std::remove_reference<Arg>::type>::type>::value,"type mismatch during array conversion");
        *it = Dune::TypeTree::convert_arg(std::forward<Arg>(arg));
        assign_reference_pack_to_shared_ptr_array_unpack<T>(++it,std::forward<Args>(args)...);
      }

      template<typename T, std::size_t n, typename... Args>
      void assign_reference_pack_to_shared_ptr_array(array<std::shared_ptr<T>,n>& res, Args&&... args)
      {
        static_assert(sizeof...(Args) == n, "invalid number of arguments");
        return assign_reference_pack_to_shared_ptr_array_unpack<T>(res.begin(),std::forward<Args>(args)...);
      }

      // prototype and end of recursion
      template<typename T, typename It, typename... Args>
      void assign_pointer_pack_to_shared_ptr_array_unpack(It it, Args*... args) {}

      template<typename T, typename It, typename Arg, typename... Args>
      void assign_pointer_pack_to_shared_ptr_array_unpack(It it, Arg* arg, Args*... args)
      {
        static_assert(std::is_same<T,typename std::remove_const<typename std::remove_reference<Arg>::type>::type>::value,"type mismatch during array conversion");
        *it = std::shared_ptr<Arg>(arg);
        assign_pointer_pack_to_shared_ptr_array_unpack<T>(++it,args...);
      }

      template<typename T, std::size_t n, typename... Args>
      void assign_pointer_pack_to_shared_ptr_array(array<std::shared_ptr<T>,n>& res, Args*... args)
      {
        static_assert(sizeof...(Args) == n, "invalid number of arguments");
        return assign_pointer_pack_to_shared_ptr_array_unpack<T>(res.begin(),args...);
      }

    } // anonymous namespace


        template<class T, class Imp>
        class TransportMulticomponentInterface;

        template<class T, class Imp>
        class TransportMulticomponentTemporalInterface;


        template<typename GV, typename RF>
        struct TransportParameterTraits;


      /** \brief Collect k instances of type T within a \ref TypeTree.
       *
       *  \tparam T The base type
       *  \tparam k The number of instances this node should collect
       */

    //! Base class for composite nodes based on variadic templates.
    template<typename GV, typename RF, typename T, std::size_t k>
    class MulticomponentTransport:
      public Dune::PDELab::TransportMulticomponentInterface<Dune::PDELab::TransportParameterTraits<GV,RF>, MulticomponentTransport<GV,RF,T,k> >,
      public Dune::PDELab::TransportMulticomponentTemporalInterface<Dune::PDELab::TransportParameterTraits<GV,RF>, MulticomponentTransport<GV,RF,T,k> >
    {

    public:
      typedef Dune::PDELab::TransportParameterTraits<GV,RF> Traits;

      //! The number of Components.
      static const std::size_t COMPONENTS = k;

      //! The type of each component.
      typedef T ComponentType;

      //! The storage type of each component.
      typedef std::shared_ptr<T> ComponentStorageType;

      //! The const version of the storage type of each component.
      typedef std::shared_ptr<const T> ComponentConstStorageType;

      //! The type used for storing the Components.
      typedef array<ComponentStorageType,k> NodeStorage;


      //! Access to the type and storage type of the i-th component.
      template<std::size_t i>
      struct Component
      {

        static_assert((i < COMPONENTS), "component index out of range");

        //! The type of the component.
        typedef T Type;

        //! The storage type of the component.
        typedef ComponentStorageType Storage;

        //! The const storage type of the component.
        typedef ComponentConstStorageType ConstStorage;
      };



      T& component (std::size_t i)
      {
        assert(i < COMPONENTS && "component index out of range");
        return *_Components[i];
      }

      //! Returns the i-th component (const version).
      /**
       * \returns a const reference to the i-th component.
       */
      const T& component (std::size_t i) const
      {
        assert(i < COMPONENTS && "component index out of range");
        return *_Components[i];
      }


      //! Returns the storage of the i-th component.
      /**
       * \returns a copy of the object storing the i-th component.
       */
      ComponentStorageType componentStorage(std::size_t i)
      {
        assert(i < COMPONENTS && "component index out of range");
        return _Components[i];
      }

      ComponentConstStorageType componentStorage (std::size_t i) const
      {
        assert(i < COMPONENTS && "component index out of range");
        return (_Components[i]);
      }

      //! Sets the i-th component to the passed-in value.
      void setComponent (std::size_t i, T& t)
      {
        assert(i < COMPONENTS && "component index out of range");
        _Components[i] = stackobject_to_shared_ptr(t);
      }

      //! Sets the stored value representing the i-th component to the passed-in value.
      void setComponent (std::size_t i, ComponentStorageType st)
      {
        assert(i < COMPONENTS && "component index out of range");
        _Components[i] = st;
      }

      //! Sets the stored value representing the i-th component to the passed-in value.
      std::string getName (std::size_t i)
      {
        assert(i < COMPONENTS && "component index out of range");
        return _Components[i]->getName();
      }

      void setTime(RF time)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->setTime(time);
      }

      void setTimeTarget(RF time, RF dt)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->setTimeTarget(time,dt);
      }

      void preStep(RF time, RF dt, int stages)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->preStep(time,dt,stages);
      }

      // constructor for objects
      template<typename... Components>
      MulticomponentTransport (Components&&... components)
      {
        assign_reference_pack_to_shared_ptr_array(_Components,std::forward<Components>(components)...);
      }

      // constructor for pointers!!
      template<typename... Components>
      MulticomponentTransport (Components*... components)
      {
        assign_pointer_pack_to_shared_ptr_array(_Components,std::forward<Components*>(components)...);
      }

      //! @}

    private:
      NodeStorage _Components;
    };



    template<class T, class Imp>
    class DiffusionMulticomponentInterface;

    template<typename GV, typename RF>
    struct DiffusionParameterTraits;


    /** \brief Collect k instances of type T within a \ref TypeTree.
     *
     *  \tparam T The base type
     *  \tparam k The number of instances this node should collect
     */

    //! Base class for composite nodes based on variadic templates.
    template<typename GV, typename RF, typename T, std::size_t k>
    class MulticomponentDiffusion:
      public Dune::PDELab::DiffusionMulticomponentInterface<Dune::PDELab::DiffusionParameterTraits<GV,RF>, MulticomponentDiffusion<GV,RF,T,k> >
    {

    public:
      typedef typename T::Traits Traits;

      //! The number of Components.
      static const std::size_t COMPONENTS = k;

      //! The type of each component.
      typedef T ComponentType;

      //! The storage type of each component.
      typedef shared_ptr<T> ComponentStorageType;

      //! The const version of the storage type of each component.
      typedef shared_ptr<const T> ComponentConstStorageType;

      //! The type used for storing the Components.
      typedef array<ComponentStorageType,k> NodeStorage;


      //! Access to the type and storage type of the i-th component.
      template<std::size_t i>
      struct Component
      {

        static_assert((i < COMPONENTS), "component index out of range");

        //! The type of the component.
        typedef T Type;

        //! The storage type of the component.
        typedef ComponentStorageType Storage;

        //! The const storage type of the component.
        typedef ComponentConstStorageType ConstStorage;
      };



      T& component (std::size_t i)
      {
        assert(i < COMPONENTS && "component index out of range");
        return *_Components[i];
      }

      //! Returns the i-th component (const version).
      /**
       * \returns a const reference to the i-th component.
       */
      const T& component (std::size_t i) const
      {
        assert(i < COMPONENTS && "component index out of range");
        return *_Components[i];
      }


      //! Returns the storage of the i-th component.
      /**
       * \returns a copy of the object storing the i-th component.
       */
      ComponentStorageType componentStorage(std::size_t i)
      {
        assert(i < COMPONENTS && "component index out of range");
        return _Components[i];
      }

      ComponentConstStorageType componentStorage (std::size_t i) const
      {
        assert(i < COMPONENTS && "component index out of range");
        return (_Components[i]);
      }

      //! Sets the i-th component to the passed-in value.
      void setComponent (std::size_t i, T& t)
      {
        assert(i < COMPONENTS && "component index out of range");
        _Components[i] = stackobject_to_shared_ptr(t);
      }

      //! Sets the stored value representing the i-th component to the passed-in value.
      void setComponent (std::size_t i, ComponentStorageType st)
      {
        assert(i < COMPONENTS && "component index out of range");
        _Components[i] = st;
      }

      void setTime(RF time)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->setTime(time);
      }

      void setTimeTarget(RF time, RF dt)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->setTimeTarget(time,dt);
      }

      void preStep(RF time, RF dt, int stages)
      {
        for (std::size_t i = 0; i<COMPONENTS; i++)
          _Components[i]->preStep(time,dt,stages);
      }

      // constructor for objects
      template<typename... Components>
      MulticomponentDiffusion (Components&&... components)
      {
        assign_reference_pack_to_shared_ptr_array(_Components,std::forward<Components>(components)...);
      }

      // constructor for pointers!!
      template<typename... Components>
      MulticomponentDiffusion (Components*... components)
      {
        assign_pointer_pack_to_shared_ptr_array(_Components,std::forward<Components*>(components)...);
      }

      //! @}

    private:
      NodeStorage _Components;
    };



  } // namespace PDELab
} //namespace Dune

#endif
