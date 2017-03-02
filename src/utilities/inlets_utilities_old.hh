// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifdef DUNE_DYCAP_INLETSUTILITIES_HH
#warning ***** WARNING ***** inlets_utilities.hh was already included ******
#endif

#ifndef DUNE_DYCAP_INLETSUTILITIES_HH
#define DUNE_DYCAP_INLETSUTILITIES_HH

#include <sys/types.h>
#include <sys/stat.h>

#include<iostream>
#include<vector>
#include<map>
#include<dune/common/exceptions.hh>
#include<dune/common/parametertreeparser.hh>
#include <dune/common/timer.hh>
#include<dune/grid/yaspgrid.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/universalmapper.hh>
#include <dune/grid/io/file/vtk/boundarywriter.hh>
#include <dune/pdelab/experimental/common/vtkexport.hh>

#include"../parameters/heleshawdomain.hh"
#include"gridfunction_utilities.hh"
#include <dune/pdelab/common/geometrywrapper.hh>
#include"inlets_utilities_old.hh"

enum Face
  {
    left=0,
    right=1,
    top = 2,
    bottom = 3,
    front = 4,
    back = 5
  };

template <class RF, int dim>
struct NewSpanInletsParameter
{
  const RF start;
  const RF end;
  const std::vector<RF> inletsx;
  const std::vector<RF> inletsy;
  const std::vector<RF> inletsz;
  std::vector<RF> inletsize;
  std::vector<RF> initial;
  std::vector<RF> input;
  const std::vector<RF> flux;

  NewSpanInletsParameter(RF start_, RF end_, std::vector<RF> inletsx_,std::vector<RF> inletsy_,std::vector<RF> inletsz_,std::vector<RF> inletsize_, std::vector<RF> initial_, std::vector<RF> input_, const std::vector<RF> flux_) DUNE_DEPRECATED
    : start(start_), end(end_),
      inletsx(inletsx_), inletsy(inletsy_), inletsz(inletsz_), inletsize(inletsize_), initial(initial_), input(input_),
      flux(flux_)
  {
    const size_t vsize = inletsx.size();
    if (vsize!=inletsy.size())
      {DUNE_THROW(Dune::Exception,"Size of inlets is not the same: inletsx "
                  << inletsx.size() << " inletsy "
                  << inletsy.size());}
    if (dim>2 && vsize!=inletsz.size())
      {DUNE_THROW(Dune::Exception,"Size of inlets is not the same: inletsx "
                  << inletsx.size() << " inletsy "
                  << inletsy.size() << " inletsz "
                  << inletsz.size());}

    if (vsize!=inletsize.size() && inletsize.size()!=1)
      {DUNE_THROW(Dune::Exception,"Size of inletsize vector is not 1 or is not same as dim of inlets ");}

    if (vsize!=initial.size() && initial.size()!=1)
      {DUNE_THROW(Dune::Exception,"Size of initial vector is not 1 or is not same as dim of inlets ");}

    if (vsize!=input.size() && input.size()!=1)
      {DUNE_THROW(Dune::Exception,"Size of initial vector is not 1 or is not same as dim of inlets ");}


    if (inletsize.size()==1){
      inletsize.resize(vsize);
      std::fill (inletsize.begin()+1,inletsize.end(),inletsize[0]);
    }
    if (initial.size()==1){
      initial.resize(vsize);
      std::fill (initial.begin()+1,initial.end(),initial[0]);
    }
    if (input.size()==1){
      input.resize(vsize);
      std::fill (input.begin()+1,input.end(),input[0]);
    }

  }

};

//! original inlets class
template <class Traits, class RF>
class InletsBase
{
public:
  typedef typename Traits::GridViewType GV;
  typedef typename GV::Grid::ctype DF;
  typedef typename Traits::DomainType DomainType;
  enum {dim=Traits::GridViewType::Grid::dimension};

  typedef Dune::LocalUniversalMapper<typename GV::Traits::Grid> BoundaryMapper;
  typedef NewSpanInletsParameter<RF,dim> SIParameter;
  typedef std::list<SIParameter> SpanInletsList;

  InletsBase(const GV& gv_,  const Dune::ParameterTree & param_, std::string cname_="InletsDefault") DUNE_DEPRECATED :
    param(param_),
    gv(gv_),
    domain(param.sub("Domain")),
    cname(cname_),
    fluxvector(0.),
    eps(1.e-7),
    updated(false),
    verbosity_level(param.get<int>("InletsDefault.verbosity",0.)),
    verbosity_visualize(param.get<bool>("InletsDefault.visualize",false))

  {
    if (gv.comm().rank()>0 && verbosity_level<2)
      verbosity_level = 0;

    if (verbosity_level)
      std::cout << "creating InletsBase class for " << cname << std::endl;

    boundaryMapper();
    fluxvector.resize(boundaryMapper_->size());
    if (verbosity_level)
      std::cout <<"boundary mapper created, size of fluxvector is " << fluxvector.size() << std::endl;

    typedef typename Dune::ParameterTree::KeyVector KeyVector;
    const Dune::ParameterTree intervals = param.Dune::ParameterTree::sub("intervals");
    const KeyVector keyvector = intervals.getSubKeys();
    KeyVector::const_iterator it = keyvector.begin();
    KeyVector::const_iterator eit = keyvector.end();

    std::vector<RF> inletsx_default =  param.get<std::vector<RF> >(cname+".inletsx", std::vector<RF>(1, 0.0));
    std::vector<RF> inletsy_default =  param.get<std::vector<RF> >(cname+".inletsy", std::vector<RF>(1, 0.0));
    std::vector<RF> inletsz_default =  param.get<std::vector<RF> >(cname+".inletsz", std::vector<RF>(1, 0.0));
    std::vector<RF> inletsize_default =  param.get<std::vector<RF> >(cname+".inletsize", std::vector<RF>(1, 0.0));
    std::vector<RF> initial_values_default =  param.get<std::vector<RF> >(cname+".initial", std::vector<RF>(1, 0.0));
    std::vector<RF> input_values_default =  param.get<std::vector<RF> >(cname+".input", std::vector<RF>(1, 0.0));


    if (verbosity_level)
      std::cout << "create time spans for inlets " << std::endl;


    for(; it!=eit; ++it){
      const Dune::ParameterTree sub = intervals.Dune::ParameterTree::sub(*it);

      const std::vector<RF> flux {sub.get<RF>(std::string("leftflux"),0.),
          sub.get<RF>(std::string("rightflux"),0.),
          sub.get<RF>(std::string("topflux"),0.),
          sub.get<RF>(std::string("bottomflux"),0.),
          sub.get<RF>(std::string("frontflux"),0.),
          sub.get<RF>(std::string("backflux"),0.) };

      SIParameter p(sub.get<RF>(std::string("start")),
                    sub.get<RF>(std::string("end")),

                    sub.get<std::vector<RF> >(std::string(cname+"inletsx"),inletsx_default),
                    sub.get<std::vector<RF> >(std::string(cname+"inletsy"),inletsy_default),
                    sub.get<std::vector<RF> >(std::string(cname+"inletsz"),inletsz_default),
                    sub.get<std::vector<RF> >(std::string(cname+"inletsize"),inletsize_default),
                    sub.get<std::vector<RF> >(std::string(cname+"initial"),initial_values_default),
                    sub.get<std::vector<RF> >(std::string(cname+"input"),input_values_default),
                    flux);

      if(inlets_list.size() && p.start != inlets_list.back().end)
        {DUNE_THROW(Dune::Exception,"Found inconsistent time spans. One ends at "
                    << inlets_list.back().end << " while the next begins at "
                    << p.start);}

      inlets_list.push_back(p);

      if (verbosity_level)
        {
          std::cout << "\nTime stepping span (" <<*it<< "): "<< p.start<<"-"<<p.end
                    << "\n  inlets = ";;
          for(typename std::vector<RF>::size_type n = 0; n < p.inletsx.size(); ++n)
            {
              std::cout  << "\n  x = " << p.inletsx[n] << " y = " <<  p.inletsy[n];
              if (dim>2)
                std::cout << " z = " << p.inletsz[n];
              std::cout << " with inletsize " << p.inletsize[n] << " initial " << p.initial[n] << " input " << p.input[n];
            }
          std::cout << "\nfluxes ";
          for(typename std::vector<RF>::size_type n = 0; n < p.flux.size(); ++n)
            std::cout << flux[n] << " ";
          std::cout << std::endl;

        }
    }

    inlets_iterator = inlets_list.begin();
    time = inlets_list.front().start;

    if (cname == "InletsDefault")
      cname = "twophase";
  }


  //! boundary mapper
  /**
   * Returns the boundary mapper.  If there is no boundary mapper yet,
   * initialize it first.
   */
  virtual const Dune::shared_ptr<BoundaryMapper> &boundaryMapper() {
    if(!boundaryMapper_) {
      typedef typename GV::template Codim<0>::Iterator EIterator;
      typedef typename GV::IntersectionIterator IIterator;

      boundaryMapper_.reset(new BoundaryMapper(gv.grid()));
      const EIterator &eend = gv.template end<0>();
      for(EIterator eit = gv.template begin<0>(); eit != eend; ++eit) {
        const IIterator &iend = gv.iend(*eit);
        for(IIterator iit = gv.ibegin(*eit); iit != iend; ++iit)
          if(iit->boundary() && !iit->neighbor())
            {
              // for a universal mapper, this will create a new map entry
              boundaryMapper_->map(*eit, iit->indexInInside(), 1);
            }
      }
    }

    return boundaryMapper_;
  }

  template<class T>
  inline void setTime(T time_)
  {
    // std::cout << "setTime to " << time_ <<  std::endl;
    time = time_;

    if ( (time<inlets_iterator->start) || (time>inlets_iterator->end))
      {
        updated = false;

        while(time < inlets_iterator->start){
          if(inlets_iterator == inlets_list.begin())
            {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
          --inlets_iterator;
          update();
        }

        while(time >= inlets_iterator->end){
          if (verbosity_level)
            std::cout << "inlets timespan end " << inlets_iterator->end << std::endl;
          if(inlets_iterator == inlets_list.end())
            {
              if (verbosity_level)
                std::cout << "Last Time Span, time " << time << " timespan end " << inlets_iterator->end << std::endl;
              {DUNE_THROW(Dune::Exception,"Time" << time <<" is out of any given interval");}
            }
          else
            {
              ++inlets_iterator;
              std::cout << "update inlets span to (" << inlets_iterator->start << "-" <<  inlets_iterator->end << ") "<< std::endl;
              update();
            }
        }
      }

  }



  //! visualize boundary condition to vtk file
  void visualize()
  {
    /*
    typedef P0BoundaryGridFunction<GV, RF,1,BoundaryMapper,std::vector<RF> > BTFunction;
    BTFunction btFunction(gv, *this->boundaryMapper(), fluxvector);

    std::string bname;
    bname = cname + "_" + param.get<std::string>("VTKname").c_str() + "_r" + std::to_string(static_cast<long long int>(domain.refine)) + "_t" + std::to_string(static_cast<long long int>(time));

    std::string prefix = static_cast<std::string> (bname);
    if(prefix != "") {
      if (verbosity_level)
        std::cout  << "Writing boundary condition type to vtk-file " << prefix
                   << std::endl;
      Dune::VTK::NonConformingBoundaryWriter<GV> writer(gv);
      writer.addCellData(new Dune::PDELab::VTKBoundaryGridFunctionAdapter<BTFunction> (btFunction), "boundary value");

      std::string path = "boundary";
      int status;
      status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw;

      writer.pwrite(prefix,path,"",outputtype);
      if (verbosity_level)
        std::cout <<    "Writing boundary condition type to vtk-files... done"
                  << std::endl;
    }
    */
  }

  //! update the construction values for a given timestep
  virtual void update()
  {
    DUNE_THROW(Dune::NotImplemented, "You have to implement update function for inlets");
  }

  /*! \brief evaulate the reconstructed flux with in a certain element

    \param e intersection to evaluate in
    \param x local position in the intersection
    \param time
  */
  typename Traits::RangeFieldType
  evaluate (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time=0.) const
  {

    if (!updated)
      DUNE_THROW(Dune::Exception, "New inlets, update necessary");
    //std::cout << "id " << boundaryMapper_->map( *(is.inside()),  is.indexInInside(), dim-1)<< " value " << fluxvector[boundaryMapper_->map( *(is.inside()),  is.indexInInside(), dim-1)] << std::endl;
    return  fluxvector[boundaryMapper_->map( *(is.inside()),  is.indexInInside(), 1)];

  }

  /* go through all inlets, read x,y and z coordinates of inlets
   * and volume of inlets and compute function f at given coordinate x
   */
  template<class X>
  inline int
  face (const X& x) const
  {

    if (x[0]<eps)
      return left;
    if (x[0]>domain.width-eps)
      return right;
    if (x[dim-1]<eps)
      return bottom;
    if (x[dim-1]>domain.height-eps)
      return top;
    if (dim > 2 && x[1] < eps)
      return front;
    if (dim > 2 && x[1] > domain.depth - eps)
      return back;
    DUNE_THROW(Dune::Exception,"coordinates in face detection" << x << "are not at the boundary of the domain");
    return -1;
  }

  const Dune::ParameterTree & param;

protected:
SpanInletsList inlets_list;
  typename SpanInletsList::const_iterator inlets_iterator;


  const GV& gv;
  const HeleshawDomain<dim> domain;
  std::string cname;
  std::vector<RF> fluxvector;

  RF time;
  const RF eps;
  bool updated;
  int verbosity_level;
  bool verbosity_visualize;
  Dune::Timer watch;
  Dune::shared_ptr<BoundaryMapper> boundaryMapper_;

};



//! for two phases inlets class
template <class Traits, class RF>
class TwophaseInlets
  : public InletsBase<Traits,RF>
{

public:
  typedef  InletsBase<Traits,RF> Base;
  typedef typename Traits::GridViewType GV;
  //! \brief Export type for domain field

  typedef typename GV::Grid::ctype DF;
  typedef typename Traits::DomainType DomainType;

  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef Dune::PDELab::IntersectionGeometry<Intersection> IG;


  enum {dim=Traits::GridViewType::Grid::dimension};


  TwophaseInlets(const GV& gv_,  const Dune::ParameterTree & param_) DUNE_DEPRECATED :
    Base(gv_, param_),
    type(param_.get<std::string>("InletsDefault.type", "default")),
    sigma(param_.get<RF>("InletsDefault.sigma",-1.)),
    bump_k(param_.get<RF>("InletsDefault.bumpk",-1.)),
    bump_n(param_.get<RF>("InletsDefault.bumpn",-1.)),
    height_watertable( param_.get<RF>("Setup.waterheight")),// initial hight of water
    porosity(param_.sub(param_.sub("Setup").get<std::string>("material")).get<RF>("porosity")), // medium porosity
    waterheightcorrection(param_.get<bool>("InletsDefault.waterheightcorrection",false))
  {
    if (verbosity_level)
      std::cout <<"water height control " << waterheightcorrection << std::endl;

    update();
    if (verbosity_level)
      std::cout << "inlets for twophase constructed" << std::endl;
  }


  /* go through all inlets, read x,y and z coordinates of inlets
   * and volume of inlets and compute function f at given coordinate x
   */
  template<class X>
  inline typename Traits::RangeFieldType
  f (const X& x) const
  {

    RF value(0);
    X c(0);

    for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
      {
        c[0] = inlets_iterator->inletsx[i];
        c[dim-1] = inlets_iterator->inletsy[i];
        if (dim>2)
          c[1] = inlets_iterator->inletsz[i];
        RF r = inlets_iterator->inletsize[i]/2.;
        value+=f_value(x,c,r);
      }
    return value;
  }

  template <class X>
  inline typename Traits::RangeFieldType
  f_value (const X& x, X& c, const RF r) const
  {
    c-=x;
    RF d = c.two_norm();
    if (type == "default")
      return (d<=r)?(1.0*water_correction(x)):0.0;
    else if (type == "circle")
      return (d<=r)?(water_correction(x)*std::sqrt(r*r-d*d)):0.0;
    else if (type == "bump")
      return  (d<=r)?(water_correction(x)*std::pow(1./bump_k,bump_n)*std::exp(-1.*std::pow(1./bump_k,2.)/(1-std::pow(d/r,2.)))):0.0;
    if (type == "gauss")
      {
        if (sigma==-1)
          DUNE_THROW(Dune::Exception, "inlet type gauss need to define sigma!");
        return (d<=r)?(water_correction(x)*exp(-(d)*(d)/(2.*sigma*sigma))):0.0;
      }
    DUNE_THROW(Dune::Exception, "inlet type should be default, circle, bump or gauss!!");
  }

  template<class X>
  inline typename Traits::RangeFieldType
  water_correction (const X& x) const
  {
    size_t face_id = face(x);
    if (!waterheightcorrection)
      return 1.;
    if (face_id !=2 && face_id !=3)
      {
        RF v1 = (patm - std::pow(-1.,face_id) *height_watertable*Rho_l*9.81*porosity)/patm;
        if (x[dim-1] > height_watertable)
          return 1;
        else
          return (1.-v1)/height_watertable*x[dim-1]+v1;
      }
    else return 1.;
  }



  //! update the construction values for a given timestep
  void update()
  {
    if (updated) {
      if (verbosity_level)
        std::cout << "inlets update not necessary" << std::endl;
      return;
    }

    if (verbosity_level){
      std::cout << "updating inlets";
      if (gv.comm().size()>1)
        std::cout << " on processor " << gv.comm().rank();
      std::cout << std::endl;
    }

    watch.reset();


    std::vector<RF> scale;
    scale.resize(4);
    if (dim>2)
      scale.resize(6);

    std::fill(scale.begin(), scale.end(), 0.);
    scale_factor(scale);
    // if (gv.comm().size()>1)
    //   for (size_t i=0;i<scale.size();i++)
    //    scale[i] =  gv.comm().sum(scale[i]);

    // vector for control of fluxes
    std::vector<RF> fluxvectorcontrol(fluxvector.size(),0.);
    std::vector<RF> control(6,0.0);

    if (verbosity_level)
      {
        if (gv.comm().size()>1 && verbosity_level>1)
          std::cout << "scale on processor " << gv.comm().rank() << std::endl;
        std::cout << "scale left is " << scale[left] << std::endl;
        std::cout << "scale right is " << scale[right] << std::endl;
        std::cout << "scale bottom is " << scale[bottom] << std::endl;
        std::cout << "scale top is " << scale[top] << std::endl;
        if (dim>2)
          {
            std::cout << "scale front is " << scale[front] << std::endl;
            std::cout << "scale back is " << scale[back] << std::endl;
          }
      }

    std::fill(fluxvector.begin(), fluxvector.end(), 0.);

    // loop over cels
    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        unsigned int intersection_index = 0;
        IntersectionIterator iit = gv.ibegin(*it);
        IntersectionIterator eiit = gv.iend(*it);

        // loop over faces
        for(; iit!=eiit; ++iit, ++intersection_index)
          {
            // skeleton term
            if(!iit->boundary()) continue;

            // intersection
            typedef typename Dune::PDELab::IntersectionGeometry<Intersection> IG;
            IG ig(*iit,intersection_index);

            // face geometry
            const Dune::FieldVector<DF,IG::dimension-1>&
              face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            RF face_volume = ig.geometry().volume();

            Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
              global = ig.intersection().geometry().global(face_local);
            size_t id = boundaryMapper_->map(*ig.inside(), ig.indexInInside(), 1);
            RF faceintegral = average_at_intersection(ig)*scale[face(global)]*inlets_iterator->flux[face(global)];

            fluxvector[id] -= faceintegral;
            if (inlets_iterator->flux[face(global)]!=0.)
              fluxvectorcontrol[id]-=faceintegral*face_volume/inlets_iterator->flux[face(global)];
          }
      }


    updated=true;
    if (verbosity_level)
      std::cout << "inlets timespan update ... done : " << watch.elapsed() << " s" << std::endl;

    if (verbosity_visualize){
      if (verbosity_level)
        std::cout << "visualize new boundary conditions " << std::endl;
      this->visualize();
    }

    // loop over cels
    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        unsigned int intersection_index = 0;
        IntersectionIterator iit = gv.ibegin(*it);
        IntersectionIterator eiit = gv.iend(*it);

        // loop over faces
        for(; iit!=eiit; ++iit, ++intersection_index)
          {
            // skeleton term
            if(!iit->boundary()) continue;

            // intersection
            typedef typename Dune::PDELab::IntersectionGeometry<Intersection> IG;
            IG ig(*iit,intersection_index);


            // face geometry
            const Dune::FieldVector<DF,IG::dimension-1>&
              face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

            Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
              global = ig.intersection().geometry().global(face_local);

            size_t id = boundaryMapper_->map(*ig.inside(), ig.indexInInside(), 1);
            control[face(global)] +=fluxvectorcontrol[id];
          }
      }

    if (gv.comm().size()>1)
      for (size_t i=0;i<control.size();i++)
        control[i] =  gv.comm().sum(control[i]);

    RF control_volume(0.);
    for (auto it = fluxvectorcontrol.begin(); it != fluxvectorcontrol.end(); ++it)
      {
        control_volume+=*it;
      }

     control_volume = gv.comm().sum(control_volume);



    for (auto it = control.begin(); it != control.end(); ++it)
      {
        if (verbosity_level>1)
        std::cout << "rank " << gv.comm().rank() << " face " << it-control.begin() << " value " << *it << "\n";
      }

    if (verbosity_level)
    std::cout << "\ncontrol " << control_volume << std::endl;

}

private:

  //! computes an average value of function f at intersection ig
  template<class IG>
  RF average_at_intersection(const IG& ig) const
  {

    // face geometry
    const Dune::FieldVector<DF,IG::dimension-1>&
      face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
    RF face_volume = ig.geometry().volume();

    Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
      global = ig.intersection().geometry().global(face_local);

    RF integral_on_face = 0.;

    Dune::GeometryType gt = ig.geometry().type();
    const int intorder = 50;
    const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,intorder);
    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
      {

        Dune::FieldVector<DF,dim> position = ig.geometry().global(qit->position());
        RF fval = f(ig.geometry().global(qit->position()));
        RF weight = qit->weight();
        RF detjac = ig.geometry().integrationElement(qit->position());
        integral_on_face +=fval*weight*detjac;
      }
    return  (integral_on_face/face_volume);
  }

  //! scales function f that the integral of f over face is equal to volume of face
  void scale_factor(std::vector<RF>& scale) const
  {
    std::fill(scale.begin(), scale.end(), 0.);

    // vector, for each inlet stores integral over inlet for function f
    std::vector<RF> volume(scale.size(),0.);

    // loop over cels
    for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
      {
        if (it->partitionType()==Dune::OverlapEntity)
          continue;

        unsigned int intersection_index = 0;
        IntersectionIterator iit = gv.ibegin(*it);
        IntersectionIterator eiit = gv.iend(*it);

        // loop over faces
        for(; iit!=eiit; ++iit, ++intersection_index)
          {
            // skeleton term
            if(!iit->boundary()) continue;

            // intersection
            typedef typename Dune::PDELab::IntersectionGeometry<Intersection> IG;
            IG ig(*iit,intersection_index);

            // face geometry
            const Dune::FieldVector<DF,IG::dimension-1>&
              face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
            //   RF face_volume = ig.geometry().volume();

            Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
              global = ig.intersection().geometry().global(face_local);

            int s = face(global);

            Dune::GeometryType gt = ig.geometry().type();
            const int intorder = 50;
            const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,intorder);
            // loop over quadrature points
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
              {

                Dune::FieldVector<DF,dim> position = ig.geometry().global(qit->position());

                RF fval = f(position);
                RF weight = qit->weight();
                RF detjac = ig.geometry().integrationElement(qit->position());
                scale[s] +=fval*weight*detjac;
              }
          }
      }

    std::vector<RF> inlets_volume(scale.size(),0.0);

    //it is global for all processors!!!
    for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
      {
        Dune::FieldVector<RF,dim> c;
        c[0] = inlets_iterator->inletsx[i];
        c[dim-1] = inlets_iterator->inletsy[i];
        if (dim>2)
          c[1] = inlets_iterator->inletsz[i];
        RF r = inlets_iterator->inletsize[i];
        RF S = 2*r/2.;
        if (dim>2)
          S = M_PI * r/2. * r/2.;
        inlets_volume[face(c)]+=S;
      }


    std::cout << "rank " << gv.comm().rank() << " scale " << scale[bottom] << " volume " << inlets_volume[bottom] << std::endl;
    for (size_t i = 0; i<scale.size();i++)
      {
        // std::cout << "rank " << gv.comm().rank() << " scale " << scale[i] << " volume " << inlets_volume[i] << std::endl;

        if (gv.comm().size()>1)
          {
            scale[i]=gv.comm().sum(scale[i]);
          }
      }

    for (size_t i = 0; i<scale.size();i++)
      {
        if (std::abs(scale[i])<eps)
          scale[i]=0;
        else
          scale[i]=(inlets_volume[i]/scale[i]);

        //if (verbosity_level>1)
        //std::cout << "rank " << gv.comm().rank() << " scale " << i << " " << scale[i] << std::endl;
      }
  }


private:
  std::string type;
  const RF sigma;
  const RF bump_k;
  const RF bump_n;
  const RF height_watertable;
  const RF porosity;
  const bool waterheightcorrection;
  using Base::inlets_iterator;
  using Base::gv;
  using Base::domain;
  using Base::fluxvector;
  using Base::time;
  using Base::eps;
  using Base::updated;
  using Base::verbosity_level;
  using Base::verbosity_visualize;
  using Base::watch;

  using Base::boundaryMapper_;
};


//! original inlets class
template <class Traits, class RF>
class TransportInlets:
  public InletsBase<Traits,RF>
{
public:
  typedef  InletsBase<Traits,RF> Base;
  typedef typename Traits::GridViewType GV;
  //! \brief Export type for domain field

  typedef typename GV::Grid::ctype DF;
  typedef typename Traits::DomainType DomainType;

  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename Dune::PDELab::ElementGeometry<Element> EG;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;
  typedef Dune::PDELab::IntersectionGeometry<Intersection> IG;


  enum {dim=Traits::GridViewType::Grid::dimension};


  TransportInlets(const GV& gv_,  const Dune::ParameterTree & param_, const std::string cname_) DUNE_DEPRECATED :
    Base(gv_,param_,cname_),
    igstorage()

  {
    if (verbosity_level)
      std::cout << "create object for transport inlets" << std::endl;
    for (ElementIterator it = gv.template begin<0>();
         it!=gv.template end<0>(); ++it)
      {
        //typedef typename IntersectionIterator::Intersection Intersection;
        //IntersectionGeometry<Intersection>(*iit,intersection_index)
        unsigned int intersection_index = 0;
        // skeleton term
        IntersectionIterator endit = gv.iend(*it);
        for (IntersectionIterator iit = gv.ibegin(*it); iit!=endit; ++iit,++intersection_index)
          {
            if (!iit->boundary())
              continue;
            else
              {
                igstorage.push_back(iit);
                // IG ig(*(iit),intersection_index);
                // size_t id = boundaryMapper_->map(*ig.inside(), ig.indexInInside(), 1);
                //   std::cout << "ading side " << id  << std::endl;
              } // boundary
          } // skeleton
      } // elements
    update();
     if (verbosity_level)
      std::cout << "object for transport inlets was created" << std::endl;
  }




  //! update the construction values for a given timestep
  void update()
  {

    if (updated) {
      if (verbosity_level)
        std::cout << "inlets update not necessary" << std::endl;
      return;
    }

    unsigned int intersection_index = 0;
    for(auto pit=igstorage.begin();pit!=igstorage.end();++pit)
      {
        IG ig(*(*pit),intersection_index);
        // face geometry
        const Dune::FieldVector<DF,IG::dimension-1>&
          face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);

        //  RF face_volume = ig.geometry().volume();

        Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
          global = ig.intersection().geometry().global(face_local);
        size_t id = boundaryMapper_->map(*ig.inside(), ig.indexInInside(), 1);
        fluxvector[id]=f(global);
      }
    updated = true;

    if (verbosity_visualize){
      if (verbosity_level)
      std::cout << "visualize new boundary conditions " << std::endl;
      this->visualize();
    }
  }

  /* go through all inlets, read x,y and z coordinates of inlets
   * and volume of inlets and compute function f at given coordinate x
   */
  template<class X>
  inline typename Traits::RangeFieldType
  f (const X& x) const
  {

    RF value(0);
    X c(0);

    /*for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
      {
        c[0] = inlets_iterator->inletsx[i];
        c[dim-1] = inlets_iterator->inletsy[i];
        if (dim>2)
          c[1] = inlets_iterator->inletsz[i];
        RF r = inlets_iterator->inletsize[i]/2.;
        c-=x;
        RF d = c.two_norm();
        // on inflow init or input!
        value+=((d<=r)?(inlets_iterator->input[i]):(0.0));
      }
    */

    value = inlets_iterator->input[0];
    // if (value<eps)
    // return inlets_iterator->initial[0];
    return value;
  }

private:


  std::vector<IntersectionIterator> igstorage;
  using Base::inlets_iterator;
  using Base::gv;
  using Base::domain;
  using Base::fluxvector;
  using Base::time;
  using Base::eps;
  using Base::updated;
  using Base::verbosity_level;
  using Base::verbosity_visualize;

  using Base::boundaryMapper_;


};


#endif
