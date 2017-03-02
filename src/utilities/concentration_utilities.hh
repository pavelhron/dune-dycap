// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_CONCENTRATION_UTILITIES_HH
#define DUNE_DYCAP_CONCENTRATION_UTILITIES_HH

#include <dune/common/parametertreeparser.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>


//! class to compute concentration outflow for each inlet
/**
 * \brief concentration outflow, stores in text file

 * \tparam Domain     domain
 * \tparam C          discrete grid function representing concentration
 */
template <class Domain, typename C>
class ConcentrationOutflow
{
private:

  typedef typename C::Traits Traits;

  //! \brief Enum for domain dimension
  enum { dim = Traits::dimDomain };

  //! \brief the grid view
  typedef typename Traits::GridViewType GridView;

  //! \brief export type for range field
  typedef typename Traits::RangeFieldType RF;

  //! \brief export type for doman field
  typedef typename Traits::DomainFieldType DF;

  typedef typename GridView::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GridView::Traits::template Codim<0>::Entity Element;

  typedef typename Dune::PDELab::ElementGeometry<Element> EG;

  typedef typename GridView::IntersectionIterator IntersectionIterator;
  typedef typename Dune::PDELab::IntersectionGeometry<typename IntersectionIterator::Intersection> IG;

public:
  ConcentrationOutflow(const Domain & dom_, const C & c_, const Dune::ParameterTree & param_, const std::string cname = "TransportOxygenWater", const std::string name = "xx") :
    dom(dom_),
    c(c_),
    inletsize(param_.sub(cname).template get<RF>("outputinletsize")),
    inlets(param_.sub(cname).template get<std::vector<RF> >("outputinlets", std::vector<RF>(1, 0.0))),
    background_conc(param_.sub(cname).template get<RF>("input")),
    localinletsize(inlets.size(),0.0),
    igstorage()
  {
    igstorage.resize(inlets.size());
    // loop once over the grid
    for (ElementIterator it = c.getGridView().template begin<0>();
         it!=c.getGridView().template end<0>(); ++it)
      {
        // Traverse intersections and store intersections and local size of outlets
        unsigned int intersection_index = 0;
        IntersectionIterator endit = c.getGridView().iend(*it);
        for (IntersectionIterator iit = c.getGridView().ibegin(*it); iit!=endit; ++iit, ++intersection_index)
          {
            if (!iit->boundary())
              continue;
            else
              {
                // I don't know
                IG ig(*iit,intersection_index);


                Dune::FieldVector<DF,IG::dimension>
                  global = ig.geometry().center();

                for (std::size_t i=0; i<inlets.size();++i)
                  {
                    if (global[0]>dom.width-1e-6 && global[dim-1]>(inlets[i]-inletsize/2.) && global[dim-1]<(inlets[i]+inletsize/2.) )
                      {
                        igstorage[i].push_back(iit);
                        RF face_volume = ig.geometry().volume();
                        localinletsize[i]+=face_volume;
                      }
                  }

              }
          }
      }

    if (c.getGridView().comm().size()>1)
      for (std::size_t i=0; i<inlets.size();++i)
        localinletsize[i]=c.getGridView().comm().sum(localinletsize[i]);

    if (c.getGridView().comm().rank()==0)
      for (std::size_t i=0; i<inlets.size();++i)
        std::cout << "outflow inlet " << i << " with volume " << localinletsize[i] << std::endl;

    if (c.getGridView().comm().rank() == 0)
      {
        std::ostringstream sname;
        if (name == "xx")
          sname << cname << "outflow.dat";
        else
          sname << name << "outflow.dat";
        s.exceptions(std::ios_base::badbit | std::ios_base::eofbit
                     | std::ios_base::failbit);
        s.open(sname.str());
        s << std::setprecision(14) << std::scientific;
        s << "#results, oxygen outflow\n";
        s << "#outflow consists of " << inlets.size() << "inlets\n";
      }
  }

  /**
   * \brief stores concentrations at each time t in text file
   */
  template<class Time>
  void output(Time t)
  {

    std::cout << "concentration outflow" << std::endl;
    if (c.getGridView().comm().rank() == 0)
      s << t << " ";

    for (std::size_t i=0; i<inlets.size();++i)
      {
        RF norm_conc(0.);
        for(auto pit=igstorage[i].begin();pit!=igstorage[i].end();++pit)
          {
            //     IG ig = (*pit)->geometry();
            //   EG eg(*((*pit)->inside()));
            //std::cout << "same types " << std::is_same<typename EG::Entity,typename C::Traits::ElementType>::value << '\n';

            //std::cout << "concentration outflow" << i << std::endl;
            // cell geometry
            //const Dune::FieldVector<RF,dim>&
            //  cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
            //std::cout << "cell center local"  << std::endl;
            //  Dune::FieldVector<RF, dim>
            //  cell_center_global = eg.geometry().global(cell_center_local);

            typename C::Traits::RangeType value;


            //std::cout << "ccl " << cell_center_local << " " << (*pit)->geometryInInside().center() << std::endl;
            // corresponds to local coordinates in element
            c.evaluate(*(*pit)->inside(),(*pit)->geometryInInside().center(),value);
            RF face_volume = (*pit)->geometry().volume();
            norm_conc+=value*face_volume;
          }

        RF global_norm_conc = c.getGridView().comm().sum(norm_conc);
        // RF inlet_volume = 2.*inletsize;//static_cast<RF>(igstorage[i].size());
        //std::cout << "sum concentration outflow" << i << std::endl;
        if (c.getGridView().comm().rank() == 0)
          s << global_norm_conc/(localinletsize[i]*background_conc) << " ";
        std::cout << "volume " << localinletsize[i] << " conc " << global_norm_conc << std::endl;
        //std::cout << norm_conc << " " << background_conc << " " << global_inlet_volume*background_conc << std::endl;
      }
    std::cout << "concentration outflow" << std::endl;
    if (c.getGridView().comm().rank() == 0)
      s << "\n";
  }

  template<class Time, typename... DGF>
  void output(Time t, const DGF&... dgf) {
    if (c.getGridView().comm().rank() == 0)
      std::cout << "time " << t;

    static const unsigned short int argsize = sizeof...(DGF);
    std::vector<RF> value(argsize);
    std::vector<RF> value_global(argsize);

    for (std::size_t i=0; i<inlets.size();++i)
      {
        for(auto pit=igstorage[i].begin();pit!=igstorage[i].end();++pit)
          {
            EG eg(*(*pit)->inside());
            const Dune::FieldVector<RF,dim>&
              cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);

            Dune::FieldVector<RF, dim>
              cell_center_global = eg.geometry().global(cell_center_local);


            auto it = value.begin();
            DGFoutput(it,eg,dgf...);
          }
        for(auto it = value.begin(); it < value.end(); ++it)
          {
            *it/=background_conc;
          }
        for (size_t j = 0; j<value.size();++j)
          {
            value_global[j]=c.getGridView().comm().sum(value[j]);
          }
        if (c.getGridView().comm().rank() == 0)
          {
            std::cout << " inlet " << i << " values ";
            for(auto it = value.begin(); it < value.end(); ++it)
              std::cout << *it << " ";
          }
      }
    if (c.getGridView().comm().rank() == 0) std::cout << "\n";
  }

  //! output from DGF (only rank 0), ellipses
  template<typename IT, typename EG, typename HDGF, typename... DGF>
  void DGFoutput(IT& it, const EG& eg, const HDGF& hdgf, const DGF&... dgf) {

    typename HDGF::Traits::RangeType value, global_value;
    // cell geometry
    const Dune::FieldVector<RF,dim>&
      cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);
    (*it)+= value;
    ++it;
    DGFoutput(it,eg,dgf...);
  }

  //! output from DGF
  template<typename IT, typename EG, typename HDGF>
  void DGFoutput(IT &it, const EG& eg, const HDGF& hdgf) {
    typename HDGF::Traits::RangeType value, global_value;
    // cell geometry
    const Dune::FieldVector<RF,dim>&
      cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
    hdgf.evaluate(eg.entity(),cell_center_local,value);
    (*it)+= value;
  }



private:
  // functionspace related variables
  const Domain& dom;
  const C& c;   //!< grid function space
  RF inletsize;
  const std::vector<RF> inlets;
  RF background_conc;
  std::vector<RF> localinletsize;
  std::vector<std::vector<IntersectionIterator> > igstorage;
  std::ofstream s;

};




//! \brief parameter class for absolute concentration output
/**
 * as a result is concentration multiplied by saturation and porosity,
 * output is component'th concentration
 * \tparam RP ReactionParameters
 * \tparam C  Concentration Discrete Grid Function
 */

template<typename  RP, typename C>
class AbsoluteConcentration
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                   typename C::Traits::RangeFieldType,
                                   C::Traits::dimRange,
                                   typename C::Traits::RangeType>, AbsoluteConcentration<RP,C> >
{
  const RP& rp;
  const C& c;
  const int component;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                           typename C::Traits::RangeFieldType,
                                           C::Traits::dimRange,
                                           typename C::Traits::RangeType> Traits;


  AbsoluteConcentration (const RP& rp_, const C& c_, const int component_) : rp(rp_), c(c_), component(component_)  {}


  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename C::Traits::RangeType Conc, c_abs(0.0);

    c.evaluate(e,x,Conc);

    if (component == 1) // O2 concentration in liquid
      c_abs=rp.c_O2_liquid(e,x,Conc);
    if (component ==2) // O2 concentration in gas
      c_abs=rp.c_O2_gas(e,x,Conc);
    if (component ==3) // DOC concentration in liquid
      c_abs=rp.c_DOC(e,x,Conc);
    if (component ==4) // Ecoli concentration in liquid
      c_abs=rp.c_Ecoli(e,x,Conc)/(6e-10);

    y = c_abs;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return c.getGridView();
  }
};





template<typename  C>
class NormConcentration
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                                                           typename C::Traits::RangeFieldType,C::Traits::dimRange,
                                                                           typename C::Traits::RangeType>,
                                          NormConcentration<C> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                           typename C::Traits::RangeFieldType,C::Traits::dimRange,
                                           typename C::Traits::RangeType> Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,NormConcentration<C> > BaseT;
  typedef typename C::Traits::RangeType RF;

  NormConcentration (const C& c_, Dune::ParameterTree param) :
    c(c_),
    height_watertable(param.sub("Setup").get<RF>("waterheight")),
    c_equilibrium(param.sub("TransportOxygenWater").get<RF>("equil")),
    c_background(param.sub("TransportOxygenWater").get<RF>("bg"))

  {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    RF normalconcentration, concentration;
    c.evaluate(e,x,concentration);
    //  typename Traits::RangeType global = e.geometry().global(x);
    //  if (global[dim-1] < height_watertable)
    normalconcentration = (concentration - c_background)/(c_equilibrium-c_background);
    y = normalconcentration;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return c.getGridView();
  }

private:
  const C& c;
  const RF height_watertable;
  const RF c_equilibrium;
  const RF c_background;
};


template<typename  RP,typename  C>
class UnitConcentration
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                                                           typename C::Traits::RangeFieldType,C::Traits::dimRange,
                                                                           typename C::Traits::RangeType>,
                                          UnitConcentration<RP,C> >
{
public:
  typedef Dune::PDELab::GridFunctionTraits<typename C::Traits::GridViewType,
                                           typename C::Traits::RangeFieldType,C::Traits::dimRange,
                                           typename C::Traits::RangeType> Traits;

  typedef typename C::Traits::RangeType RF;
  typedef typename Traits::GridViewType GV;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  enum {dim = GV::dimension};

  UnitConcentration (const RP& rp_,const C& c_) :
    rp(rp_), c(c_)
  {
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    RF concentration, c_abs;
    c.evaluate(e,x,concentration);
    c_abs=rp.c_Ecoli(e,x,concentration);
    y = c_abs/maxConc();
  }

  RF maxConc() const
  {
    RF maxc = 0.;
    for (ElementIterator it = c.getGridView().template begin<0>();
         it!=c.getGridView().template end<0>(); ++it)
      {
        EG eg(*it);
        const Dune::FieldVector<RF,dim>&
          cell_center_local = Dune::ReferenceElements<RF,dim>::general(eg.geometry().type()).position(0,0);
        typename C::Traits::RangeType cvalue, c_abs;

        c.evaluate(eg.entity(),cell_center_local,cvalue);
        c_abs=rp.c_Ecoli(eg.entity(),cell_center_local,cvalue);

        if (c_abs>maxc)
          maxc = c_abs;
      }
    return maxc;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return c.getGridView();
  }

private:
  const RP& rp;
  const C& c;
};


#endif
