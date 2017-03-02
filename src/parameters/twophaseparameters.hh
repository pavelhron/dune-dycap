// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_TWOPHASEPARAMETERS_HH
#define DUNE_DYCAP_TWOPHASEPARAMETERS_HH

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>
#include <dune/pdelab/common/function.hh>
#include<dune/geometry/referenceelements.hh>

#include<src/physics/hydraulicparameters.hh>
#include<src/physics/physical_chemistry.hh>

#include"heleshawdomain.hh"

namespace Dune {
  namespace PDELab {
    template<class T, class Imp>
    class TwoPhaseParameterInterface;

    template<typename GV, typename RF>
    struct TwoPhaseParameterTraits;
  }
}

namespace Dune {
  namespace Dycap {
    template<class Traits, typename RF>
    class TwophaseInlets;
  }
}

template<typename GV, typename RF>
class TwoPhaseParameter
  : public Dune::PDELab::TwoPhaseParameterInterface<Dune::PDELab::TwoPhaseParameterTraits<GV,RF>,
                                                    TwoPhaseParameter<GV,RF> >
{
  static constexpr RF eps = 1E-5; // eps in geometry calculations

public:
  typedef Dune::PDELab::TwoPhaseParameterTraits<GV,RF> Traits;
  enum {dim=GV::Grid::dimension};


  typedef Dune::Dycap::TwophaseInlets<Traits,RF> Inlets;


  //! constructor
  TwoPhaseParameter(const GV& gv_, const HeleshawDomain<dim>& dom_, const Dune::ParameterTree & param_) :

    param(param_),  // parameter class
    material(param.sub("Setup").get<std::string>("material")), // material to use
    incompressible(param.sub("Setup").get<bool>("incompressible", true)), // compressible / uncompressible air
    densitydependent(param.sub("Setup").get<bool>("densitydependent", false)), // compressible / uncompressible air
    gv(gv_),
    height_watertable(param.sub("Setup").get<RF>("waterheight")), // initial hight of water
    initialsaturation(param.sub("Setup").get<RF>("initialsaturation",0.0)), // initial hight of water
    scalel(param.sub("Setup").get<RF>("scalel",1.)), // initial hight of water
    scaleg(param.sub("Setup").get<RF>("scaleg",1.)), // initial hight of water
    domain(dom_), is(gv.indexSet()), perm(is.size(0)),density(is.size(0)),
    porosity(param.sub(material).get<RF>("porosity")), // medium porosity
    temperature(param.sub("Setup").get<RF>("temperature",293.16)), // temperature
    water_conservative(param.sub("Setup").get<bool>("conservative",true)), // water conservative
    random_permeability(param.sub("Setup").get<bool>("random_permeability",false)),
    reduce_permeability(param.sub(material).get<bool>("reduce_permeability",false)),
    reduce_height(param.sub(material).get<RF>("reduce_height",0)),
    reduce_factor(param.sub(material).get<RF>("reduce_factor",1.)),
    inlets(gv_,param_)

  {

    gvector=0; gvector[dim-1]=-Dune::PM::PhysicalChemistry<RF>::GRAVITY; // gravity vector

    if (!(initialsaturation<1.))
      DUNE_THROW(Dune::Exception, "initialsaturation in Setup.initialsaturation should be 0<Sl<1");

    if (gv.comm().rank()==0 && densitydependent)
      std::cout << "Twophase parameters: gas density is component dependent!" << std::endl;

    typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
    //typedef typename Traits::DomainFieldType DF;

    RF permbar = param.sub(material).get<RF>("permbar"); // permeability of material
    RF permbarlense = param.sub(material).get<RF>("permbarlense",permbar); // permeability of material

    //if (domain.lense_width_min>0)
      std::cout << "lense: " << " xmin " << domain.lense_width_min << " xmax " << domain.lense_width_max << std::endl;

    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        int id = is.index(*it);
        Dune::FieldVector<RF, dim> global = it->geometry().center();
        perm[id]=permbarlense;
        //initial density
        density[id]=patm/(R_S*temperature);
        //for (int i=0; i<dim-1; i++)
          if (global[0]<domain.lense_width_min || global[0]>domain.lense_width_max)
            perm[id]=permbar;
        if (global[dim-1]<domain.lense_height_min || global[dim-1]>domain.lense_height_max)
          perm[id]=permbar;

        if (random_permeability)
          {
            int v  = std::rand() % 101 - 50;
            perm[id]+=permbar*0.000001*v;
          }

        if (reduce_permeability)
          {
            const RF distance = domain.height-reduce_height;
            if (global[dim-1]>distance)
              perm[id]=permbar/(reduce_height*reduce_factor)*(1-reduce_factor)*(global[dim-1])+permbar-distance*permbar/(reduce_height*reduce_factor)*(1-reduce_factor);
          }

      }

    // hydraulic parameters
    std::string model = param.sub(material).get<std::string>("model");
    bool linear_interpolation = param.sub(material).get<bool>("linearinterpolation",false);
    bool cubic_interpolation = param.sub(material).get<bool>("cubicinterpolation",false);

    if (linear_interpolation && cubic_interpolation)
      DUNE_THROW(Dune::Exception, "Only one interpolation can be used");
    if (model == "regBrooksCorey")
      {
        Dune::PM::RegBrooksCoreyParam<RF>* brooksCorey = new Dune::PM::RegBrooksCoreyParam<RF>;
        hydrParam = brooksCorey;
        brooksCorey->Setdelta(param.sub(material).get<RF>("delta",0.) );
        brooksCorey->SetLambda( param.sub(material).get<RF>("lambda",0.) );
        brooksCorey->SetTau( param.sub(material).get<RF>("tau", 0.5) );
        brooksCorey->SetEntryPressure( param.sub(material).get<RF>("pentry") );
        brooksCorey->InitRegularisation(param.sub(material).get<RF>("deltareg", 1e-3) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "BrooksCorey")
      {
        Dune::PM::BrooksCoreyParam<RF>* brooksCorey = new Dune::PM::BrooksCoreyParam<RF>;
        hydrParam = brooksCorey;
        brooksCorey->Setdelta(param.sub(material).get<RF>("delta",0.) );
        brooksCorey->SetLambda( param.sub(material).get<RF>("lambda",0.) );
        brooksCorey->SetTau( param.sub(material).get<RF>("tau", 0.5) );
        brooksCorey->SetEntryPressure( param.sub(material).get<RF>("pentry") );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(brooksCorey);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "vanGenuchten")
      {
        Dune::PM::VanGenuchtenParam<RF>* vanGenuchten = new Dune::PM::VanGenuchtenParam<RF>;
        hydrParam = vanGenuchten;
        vanGenuchten->Setdelta(param.sub(material).get<RF>("delta",0.) );
        vanGenuchten->SetAlpha( param.sub(material).get<RF>("alpha") );
        vanGenuchten->SetN( param.sub(material).get<RF>("n") );
        if (param.sub(material).hasKey("m"))
          vanGenuchten->SetM( param.sub(material).get<RF>("m") );
        vanGenuchten->SetTau( param.sub(material).get<RF>("tau", 0.5) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }
      }
    else if (model == "modVanGenuchten")
      {
        Dune::PM::ModVanGenuchtenParam<RF>* vanGenuchten = new Dune::PM::ModVanGenuchtenParam<RF>;
        hydrParam = vanGenuchten;
        vanGenuchten->Setdelta(param.sub(material).get<RF>("delta",0.) );
        vanGenuchten->SetAlpha( param.sub(material).get<RF>("alpha") );
        vanGenuchten->SetN( param.sub(material).get<RF>("n") );
        if (param.sub(material).hasKey("m"))
          vanGenuchten->SetM( param.sub(material).get<RF>("m") );
        vanGenuchten->SetTau( param.sub(material).get<RF>("tau", 0.5) );

        if (linear_interpolation)
          {
            Dune::PM::HydrLinearInterpolation<RF>* li =
              new Dune::PM::HydrLinearInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }

        if (cubic_interpolation)
          {
            Dune::PM::HydrCubicInterpolation<RF>* li =
              new Dune::PM::HydrCubicInterpolation<RF>(vanGenuchten);
            hydrParam = li;
            li->Init();
          }
      }
    else
      DUNE_THROW(Dune::Exception, "Unknown hydraulic model " + model);
    assert(hydrParam);
  }

  //! porosity
  typename Traits::RangeFieldType
  phi (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return porosity;
  }

  //! porosity
  typename Traits::RangeFieldType
  phi () const
  {
    return porosity;
  }

  //! capillary pressure function
  typename Traits::RangeFieldType
  pe (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! capillary pressure function
  typename Traits::RangeFieldType
  pc (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType s_l) const
  {
    return pc(s_l);
  }

  typename Traits::RangeFieldType
  pc (typename Traits::RangeFieldType s_l) const
  {
    return hydrParam->Pc(s_l);
  }

  //! inverse capillary pressure function
  typename Traits::RangeFieldType
  s_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType pc) const
  {
    return s_l(pc);
  }

  //! inverse capillary pressure function
  typename Traits::RangeFieldType
  s_l (typename Traits::RangeFieldType pc) const
  {
    return hydrParam->Sl(pc);
  }

  //! liquid phase relative permeability
  typename Traits::RangeFieldType
  kr_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType s_l) const
  {
    return hydrParam->KRelL(s_l);
  }

  //! gas phase relative permeability
  typename Traits::RangeFieldType
  kr_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType s_g) const
  {
    return hydrParam->KRelG(s_g);
  }

  //! liquid phase viscosity
  typename Traits::RangeFieldType
  mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_l) const
  {
    return Dune::PM::PhysicalChemistry<RF>::ViscosityWater(temperature);
  }

  //! liquid phase viscosity
  typename Traits::RangeFieldType
  mu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return Dune::PM::PhysicalChemistry<RF>::ViscosityWater(temperature);
  }


  //! liquid surface tension
  typename Traits::RangeFieldType
  sigma_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return Dune::PM::PhysicalChemistry<RF>::RelSurfaceTensionWaterAir(temperature);
  }

  //! liquid surface tension
  typename Traits::RangeFieldType
  sigma_l () const
  {
    return Dune::PM::PhysicalChemistry<RF>::RelSurfaceTensionWaterAir(temperature);
  }

  //! gas phase viscosity
  typename Traits::RangeFieldType
  mu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_g) const
  {
    return Dune::PM::PhysicalChemistry<RF>::ViscosityAir(temperature);
  }

  //! absolute permeability (scalar!)
  typename Traits::RangeFieldType
  k_abs (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    int id = is.index(e);
    return perm[id];
  }

  //! gravity vector
  const typename Traits::RangeType& gravity () const
  {
    return gvector;
  }

  //! liquid phase molar density [mol/m^3]
  typename Traits::RangeFieldType
  nu_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_l) const
  {
    return 1.;//rho_l(e,x,p_l)/M_l;
  }

  //! gas phase molar density [mol/m^3]
  typename Traits::RangeFieldType
  nu_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
        typename Traits::RangeFieldType p_g) const
  {
    // incompressible, pressure of air is atmospheric pressure
    if (incompressible)
      return patm/(R_S*temperature);
    return (p_g+patm)/(R_S*temperature);
  }

  //! liquid phase mass density [kg/m^3]
  typename Traits::RangeFieldType
  rho_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType p_l) const
  {
    return Rho_l;
  }

  //! gas phase mass density [kg/m^3]
  typename Traits::RangeFieldType
  rho_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
         typename Traits::RangeFieldType p_g) const
  {
    // incompressible, pressure of air is atmospheric pressure
    if (densitydependent)
      {
        int id = is.index(e);
        return density[id];
      }
    if (incompressible)
      return patm/(R_S*temperature);
    return (p_g+patm)/(R_S*temperature);
  }

  //! liquid scale
  typename Traits::RangeFieldType
  scale_l () const
  {
    return   scalel;//(perm[0]*Dune::PM::PhysicalChemistry<RF>::ViscosityWater(temperature));
  }

  //! gas scale
  typename Traits::RangeFieldType
  scale_g () const
  {
    return scaleg;
  }


  //! liquid phase boundary condition type
  int
  bc_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);

    if (water_conservative)
      return 0;

    // top dirichlet bc
    if (global[dim-1]>domain.height-eps /*&& global[0] < 0.05*/)
      {
        return 1;
      }

    // front / back
    if (dim==3)
      {
        if (global[1]<eps || global[1]>domain.depth-eps)
          return 0; // left & right boundary: Neumann
      }

    return 0;
  }

  //! gas phase boundary condition type
  int
  bc_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);

    // top / bottom
    if (global[dim-1]>domain.height-eps)
      return 1; // top boundary Dirichlet

    // front / back
    if (dim==3)
      {
        if (global[1]<eps || global[1]>domain.depth-eps)
          return 0; // left & right boundary: Neumann
      }

    // else Neumann
    return 0;
  }

  //! liquid phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);
    // atm pressure
    return g_l_global(global,time);
  }

  //! liquid phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_l_global (const typename Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>& global,
              typename Traits::RangeFieldType time) const
  {
    if (initialsaturation>0)
      return g_g_global(global,time)-pc(initialsaturation);

    else
      return (height_watertable - global[dim-1])*(-gvector[dim-1]*Rho_l)/*+patm*/;  // hydrostatic pressure
    //      //std::max((height_watertable - global[dim-1])*(-gvector[dim-1]*Rho_l),0.)/*+patm*/;  // hydrostatic pressure
  }
  //
  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = is.geometry().global(x);
    return g_g_global(global, time);
  }

  //! gas phase Dirichlet boundary condition
  typename Traits::RangeFieldType
  g_g_global (const typename Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>& global,
              typename Traits::RangeFieldType time) const
  {
    return 0.0; //I would say the hydraulic pressure of gas is in this are negligible
    // this was variant used before
    //return /*patm +*/ (height_watertable - global[dim-1])*(-gvector[dim-1]*patm)/(R_S*temperature);
  }

  //! liquid phase Neumann boundary condition
  typename Traits::RangeFieldType
  j_l (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    // cell geometries
    //  Typename Traits::DomainType cell_center_local = Dune::GenericReferenceElements< typename Traits::RangeFieldType,dim>::general(is.inside()->type()).position(0,0);

    // typename Traits::RangeFieldType nu = nu_l(*(is.inside()),cell_center_local,0.0);

    return inlets.evaluate(is,x,time);
  }

  //! gas phase Neumann boundary condition
  typename Traits::RangeFieldType
  j_g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  //! liquid phase source term
  typename Traits::RangeFieldType
  q_l (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  //! gas phase source term
  typename Traits::RangeFieldType
  q_g (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType time) const
  {
    return 0.0;
  }

  // set timetarget for inlets
  void preStep (RF time_, RF dt_, int stages)
  {
    inlets.setTime(time_+dt_);
  }

  void setScaleL(RF scalel_)
  {
    scalel = scalel_;
  }

  void setScaleG(RF scaleg_)
  {
    scaleg = scaleg_;
  }

  inline const typename Traits::GridViewType& getGridView () const
  {
    return gv;
  }

  inline const Dune::ParameterTree & getParam() const
  {
    return param;
  }

  const Dune::PM::HydrParamBase<RF>* getHydrParam() const
  {
    return hydrParam;
  }

  const RF getWaterTable() const
  {
    return height_watertable;
  }

  const RF getInitialSaturation() const
  {
    return initialsaturation;
  }

  std::vector<RF>& getDensity(){return density;}

  const Dune::ParameterTree & param;
  const std::string material;


private:
  bool incompressible;
  bool densitydependent;
  const GV& gv;
  RF height_watertable;
  const RF initialsaturation;
  RF scalel;
  RF scaleg;
  typename Traits::RangeType gvector;
  const HeleshawDomain<dim>& domain;
  const typename GV::IndexSet& is;
  std::vector<RF> perm;
  std::vector<RF> density;
  const RF porosity;
  const RF temperature;
  Dune::PM::HydrParamBase<RF>* hydrParam;
  const bool water_conservative;
  const bool random_permeability;
  const bool reduce_permeability;
  const RF reduce_height;
  const RF reduce_factor;
  Inlets inlets;

};

/********************************************************************/
/************** saturation and gas pressure output ******************/
/********************************************************************/

template<typename  TP>
class PermeabilityField
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename TP::Traits::GridViewType,
                                   typename TP::Traits::RangeFieldType,
                                   1,
                                   typename TP::Traits::RangeType>,
  PermeabilityField<TP> >
{
  const TP& tp;

public:
  typedef Dune::PDELab::GridFunctionTraits<typename TP::Traits::GridViewType,
                                           typename TP::Traits::RangeFieldType,
                                           1,
                                           typename TP::Traits::RangeType> Traits;


  PermeabilityField (const TP& tp_) : tp(tp_) {}

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    //typename TP::Traits::RangeType pl_value,pc_value;
    y = tp.k_abs(e,x);
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return tp.getGridView();
  }
};


#endif
