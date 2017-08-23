// -*- tab-width: 4; c-basic-offset: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FLUXRECONSTRUCTION_HH
#define DUNE_PDELAB_FLUXRECONSTRUCTION_HH

#include <limits>
#include<dune/pdelab/finiteelementmap/monomfem.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

namespace Dune {
  namespace PDELab {

    //! don't do any kind of reconstruciton, just assume piecewise constant data
    template <class Traits>
    class DefaultFluxReconstruction
    {
    public:

      //! set time of current time step
      inline void setTime(typename Traits::RangeFieldType t)
      {}

      //! update the construction values for a given timestep
      inline void updateReconstruction(typename Traits::RangeFieldType time, int stage)
      {}

      /*! \brief evaulate the reconstructed flux with in a certain element

        \param e entity to evaluate in
        \param x local position in the entity
        \param c center value in this cell
        \param y result
      */

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            const typename Traits::RangeFieldType& c,
                            typename Traits::RangeFieldType& y) const
      {
        // default is the center value
        y = c;
      }

      //! factor for the timestep in order to fulfill the CFL condition
      inline typename Traits::RangeFieldType timestepFactor() const
      {
        return 1.0;
      }

      template<typename X>
      void prestage(X& x)
      {}
      template<typename X>
      void poststage(X& x)
      {}
    };


       //! don't do any kind of reconstruciton, just assume piecewise constant data
    template <class Traits>
    class DGFluxReconstruction
    {
    public:

      //! set time of current time step
      inline void setTime(typename Traits::RangeFieldType t)
      {}

      //! update the construction values for a given timestep
      inline void updateReconstruction(typename Traits::RangeFieldType time, int stage)
      {}

      /*! \brief evaulate the reconstructed flux with in a certain element

        \param e entity to evaluate in
        \param x local position in the entity
        \param c center value in this cell
        \param y result
      */

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            const typename Traits::RangeFieldType& c,
                            typename Traits::RangeFieldType& y) const
      {
        // default is the center value
        y = c;
      }

      //! factor for the timestep in order to fulfill the CFL condition
      inline typename Traits::RangeFieldType timestepFactor() const
      {
        return 0.5;
      }
      template<typename X>
      void prestage(X& x)
      {}
      template<typename X>
      void poststage(X& x)
      {}
    };

    //! second order flux reconstruction with minmod slope limiter in cartesian grids
    template <class Traits, typename F>
    class SecondOrderFluxReconstruction
    {
    private:
      typedef F FunctionU;
      enum { dim = Traits::dimDomain };
      typedef typename Traits::GridViewType GridView;
      typedef typename Traits::RangeFieldType RangeFieldType;
      typedef typename Traits::DomainFieldType DomainFieldType;
      typedef Dune::PDELab::MonomLocalFiniteElementMap<DomainFieldType, RangeFieldType, dim, 1> rFEM;
      typedef Dune::PDELab::GridFunctionSpace<GridView, rFEM> rGFS;
      typedef typename Dune::PDELab::Backend::Vector<rGFS,RangeFieldType>::Type rVEC;
      typedef Dune::PDELab::LocalFunctionSpace<rGFS> rLFS;
      typedef Dune::PDELab::LFSIndexCache<rLFS> rLFSCache;
      typedef Dune::PDELab::DiscreteGridFunction<rGFS,rVEC> rDGF;

      static GeometryType getCube() {
        GeometryType gt;
        gt.makeCube(dim);
        return gt;
      }
    public:
      SecondOrderFluxReconstruction(const FunctionU & _u, RangeFieldType _t = 1.0, const std::string model_="minmod") :
        u(_u),
        fem(getCube()),
        gfs(u.getGridView(),fem),
        lfs(gfs),
        lfs_cache(lfs),
        vec(gfs),
        vec_view(vec),
        r(gfs,vec),
        time1(0.0),
        time2(0.0),
        theta(_t),
        model(model_)
      {
        std::cout << "second order reconstruction with theta " << _t << std::endl;
      }

      //! set time of current time step
      void setTime(RangeFieldType t)
      {
        time2 = t;
      }

      template<typename X>
      void prestage(X& x)
      {}
      template<typename X>
      void poststage(X& x)
      {}




      //! update the construction values for a given timestep
      void updateReconstruction(RangeFieldType t, int stages)
      {
        // current time
        time1 = t;
        assert(stages == 1);

        // store slopes
        Dune::FieldVector< Dune::FieldVector<RangeFieldType, 1>, 2*dim+1> values;
        Dune::FieldVector<RangeFieldType, dim> slopes(0.0);

        // for each cell
        typedef typename GridView::template Codim <0>::Iterator ElementIterator;
        ElementIterator it = gfs.gridView().template begin<0>();
        ElementIterator endit = gfs.gridView().template end<0>();
        for (; it != endit; ++it)
          {
            /*
              as we have structured grid, the direction and the axis are computed as follows:

              axis = face/2
              dir = face%2

              where
              face \in [0, 2*dim[
              axis \in [0, dim[
              dir \in {0, 1}.

              the cell value is stored in values[0]
              the neighboring values are stored as values[1+axis*2+dir]
            */

            // compute cell value
            values = std::numeric_limits<RangeFieldType>::quiet_NaN();
            u.evaluate(*it, it->geometry().center(), values[0]);
            // compute neighbor values
            typedef typename GridView::IntersectionIterator IntersectionIterator;
            IntersectionIterator endit = gfs.gridView().iend(*it);
            IntersectionIterator iit = gfs.gridView().ibegin(*it);
            for (; iit!=endit; ++iit)
              {
                int face = iit->indexInInside();
                int axis = face/2;
                int dir = face%2;

                if (iit->neighbor())
                  u.evaluate(*(iit->outside()), iit->outside()->geometry().center(), values[1+axis*2+dir]);
                else
                  values[1+axis*2+dir] = values[0];
              }

            // compute slopes for each neighbor
            // and apply slope limiter
            for (int a = 0; a<dim; a++)
              {
                if (model=="minmod")
                  slopes[a] = minmod(values[0] - values[1+2*a], values[1+2*a+1] - values[0]);
                else if (model=="minmodchanged")
                  slopes[a] = minmod(values[0] - values[1+2*a],(values[1+2*a+1]-values[1+2*a])/2.0, values[1+2*a+1] - values[0]);
                else if (model=="superbee")
                  {
                    RangeFieldType s1 = minmod(2*(values[0] - values[1+2*a]), values[1+2*a+1] - values[0]);
                    RangeFieldType s2 = minmod((values[0] - values[1+2*a]), 2*(values[1+2*a+1] - values[0]));
                    slopes[a]=maxmod(s1,s2);
                  }
                else
                  DUNE_THROW(Dune::Exception, "limiter model " << model << " is not known.");


              }

            // project onto p1-dg
            lfs.bind(*it);
            lfs_cache.update();
            assert (lfs.size() == dim+1);
            std::vector<RangeFieldType> x(dim+1);
            x[0] = values[0];
            for (int a = 0; a < dim; a++)
              {
                x[a+1] = slopes[a];
                x[0] -= 0.5 * slopes[a];
              }
            //  vec.write_sub_container(lfs,x);

            vec_view.bind(lfs_cache);
            vec_view.write(x);
            vec_view.unbind();
            //lfs.unbind(*it);
            //lfs.vwrite(x, vec);
          }

        // communicate overlap
        if (gfs.gridView().comm().size() != 1 &&
            gfs.gridView().overlapSize(0) < 2)
          DUNE_THROW(Dune::Exception,
                     "SecondOrderFluxReconstruction in parallel "
                     "requires at least overlap 2.");
      }

      /*! \brief evaulate the reconstructed flux with in a certain element

        \param e entity to evaluate in
        \param x local position in the entity
        \param c center value in this cell
        \param y result
      */
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            const typename Traits::RangeFieldType& c,
                            typename Traits::RangeFieldType& y) const
      {
        // assert we are in explicit mode
#ifndef NDEBUG
        if (time1 != time2)
          {
            DUNE_THROW(Dune::Exception,
                       "SecondOrderFluxReconstruction does only work "
                       "for the explicit Euler time stepping scheme,\n"
                       "you seem to be using a different scheme, as the "
                       "evaluation time\ndoes not correspond to the "
                       "starting time of the current time step.");
          }
#endif
        // default is the center value
        FieldVector<typename Traits::RangeFieldType, 1> vn;
        r.evaluate(e,x,vn);
        y = vn[0];
      }
      /*
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeFieldType& y) const
      {
        // default is the center value
        FieldVector<typename Traits::RangeFieldType, 1> vn;
        r.evaluate(e,x,vn);
        y = vn[0];
      }
      */

      //! factor for the timestep in order to fulfill the CFL condition
      inline RangeFieldType timestepFactor() const
      {
        return 1.0 / 2.;//Traits::dimDomain;
      }

    private:

      // minmod slope limiter
      RangeFieldType minmod(RangeFieldType a, RangeFieldType b) const
      {
        if (a*b < 0.0)
          return 0.0;
        if (a > 0.0)
          return theta * std::min(a,b);
        if (a < 0.0)
          return theta * std::max(a,b);
        return 0;
      }

      // minmod slope limiter
      RangeFieldType maxmod(RangeFieldType a, RangeFieldType b) const
      {
        if (std::abs(a)>std::abs(b))
          return a;
        else
          return b;
      }


      // minmod slope limiter
      RangeFieldType minmod(RangeFieldType a, RangeFieldType b, RangeFieldType c) const
      {
        return minmod(a,minmod(b,c));
      }

      // function which is to be reconstructed
      const FunctionU & u;
      // functionspace related variables
      rFEM fem;
      rGFS gfs;
      rLFS lfs;
      rLFSCache lfs_cache;
      // data vector and function
      rVEC vec;
      typename rVEC::template LocalView<rLFSCache> vec_view;
      rDGF r;

      // just needed for assertions
      RangeFieldType time1;
      RangeFieldType time2;

      // limiter parameters
      RangeFieldType theta;
      std::string model;
    };




  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_FLUXRECONSTRUCTION_HH
