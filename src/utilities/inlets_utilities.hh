// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
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
#include <dune/pdelab/common/vtkexport.hh>

#include"../parameters/heleshawdomain.hh"
#include"gridinfo.hh"
#include"gridfunction_utilities.hh"
#include <dune/pdelab/common/geometrywrapper.hh>

namespace Dune {
  namespace Dycap{

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

      NewSpanInletsParameter(RF start_, RF end_, std::vector<RF> inletsx_,std::vector<RF> inletsy_,std::vector<RF> inletsz_,std::vector<RF> inletsize_, std::vector<RF> initial_, std::vector<RF> input_, const std::vector<RF> flux_)
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

      InletsBase(const GV& gv_,  const Dune::ParameterTree & param_, std::string cname_, std::string subintervals_, const int codim_ = 1,  bool test_boundary=true, bool silent = false):
        param(param_),
        gv(gv_),
        domain(param.sub("Domain")),
        cname(cname_),
        subintervals(subintervals_),
        codim(codim_),
        fluxvector(0.),
        eps(1.e-7),
        updated(false),
        verbosity_level(param.get<int>("InletsDefault.verbosity",0.)),
        verbosity_visualize(param.get<bool>("InletsDefault.visualize",false))

      {
        if (gv.comm().rank()>0)
          verbosity_level = 0;

        if (silent)
          verbosity_level = 0;

        if (verbosity_level)
          std::cout << "creating InletsBase class for " << cname << std::endl;

        boundaryMapper(codim);
        fluxvector.resize(boundaryMapper_->size());
        if (verbosity_level)
          std::cout <<"boundary mapper created, size of fluxvector is " << fluxvector.size() << std::endl;

        typedef typename Dune::ParameterTree::KeyVector KeyVector;
        const Dune::ParameterTree intervals = param.Dune::ParameterTree::sub(subintervals);
        const KeyVector keyvector = intervals.getSubKeys();
        KeyVector::const_iterator it = keyvector.begin();
        KeyVector::const_iterator eit = keyvector.end();

        std::vector<RF> inletsx_default =  param.get<std::vector<RF> >(cname+".inletsx", std::vector<RF>(1, 0.0));
        std::vector<RF> inletsy_default =  param.get<std::vector<RF> >(cname+".inletsy", std::vector<RF>(1, 0.0));
        std::vector<RF> inletsz_default =  param.get<std::vector<RF> >(cname+".inletsz", std::vector<RF>(1, 0.0));
        std::vector<RF> inletsize_default =  param.get<std::vector<RF> >(cname+".inletsize", std::vector<RF>(1, 0.0));
        std::vector<RF> initial_values_default =  param.get<std::vector<RF> >(cname+".initial", std::vector<RF>(1, 0.0));
        std::vector<RF> input_values_default =  param.get<std::vector<RF> >(cname+".input", std::vector<RF>(1, 0.0));


        if (param.hasKey(cname+".sideinlets"))
          {
            std::vector<RF> sideinlets =  param.get<std::vector<RF> >(cname+".sideinlets");
            inletsy_default = sideinlets;
            std::vector<RF> sideinletsize =  param.get<std::vector<RF> >(cname+".sideinletsize");
            inletsize_default = sideinletsize;


            inletsy_default.insert (inletsy_default.end(),sideinlets.begin(),sideinlets.end());
            std::vector<RF> side0(sideinlets.size(),0);
            inletsx_default=side0;
            std::fill(side0.begin(), side0.end(), domain.width);
            inletsx_default.insert (inletsx_default.end(),side0.begin(),side0.end());

            if (verbosity_level)
              {
                std::cout << "sideinlets\n";
                for (size_t i = 0; i<inletsx_default.size();++i)
                  std::cout << inletsx_default[i] << " " << inletsy_default[i] << std::endl;
              }

          }
        if (cname=="InletsDefault") cname="";
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
              std::cout << "\nfluxes in base class ";
              for(typename std::vector<RF>::size_type n = 0; n < p.flux.size(); ++n)
                std::cout << flux[n] << " ";
              std::cout << std::endl;

            }
        }

        inlets_iterator = inlets_list.begin();
        time = inlets_list.front().start;

        if (test_boundary)
        for (auto lit=inlets_list.begin();lit!=inlets_list.end();++lit)
          {
            for (auto vit = (*lit).inletsx.begin(); vit!=(*lit).inletsx.end();++vit)
              {
                std::size_t n = vit-(*lit).inletsx.begin();
                Dune::FieldVector<RF,dim> coordinate;
                coordinate[0]=(*lit).inletsx[n];
                coordinate[dim-1]=(*lit).inletsy[n];
                RF isize = (*lit).inletsize[n]/2.0;
                if (dim>2)
                  coordinate[1]=(*lit).inletsz[n];

                if (face_from_coordinates(coordinate)==0 || face_from_coordinates(coordinate)==1)
                  {
                    if (coordinate[dim-1]-isize<0. || coordinate[dim-1]+isize>(domain.height))
                      {DUNE_THROW(Dune::Exception,"x inlet in y-direction with index " << n << " is invalid (out of the domain)");}
                    if (dim>2)
                      if (coordinate[1]-isize<0. || coordinate[1]+isize>(domain.depth))
                        {DUNE_THROW(Dune::Exception,"x inlet in z-direction with index " << n << " is invalid (out of the domain)" << coordinate[1] << " " << isize);}
                  }

                if (face_from_coordinates(coordinate)==2 || face_from_coordinates(coordinate)==3)
                  {
                    if (coordinate[0]-isize<0. || coordinate[0]+isize>(domain.width))
                      {DUNE_THROW(Dune::Exception,"y inlet in x-direction with index " << n << " is invalid (out of the domain)");}
                    if (dim>2)
                      if (coordinate[1]-isize<0. || coordinate[1]+isize>(domain.depth))
                        {DUNE_THROW(Dune::Exception,"y inlet in z-direction with index " << n << " is invalid (out of the domain)");}
                  }

                if (dim>2)
                  if (face_from_coordinates(coordinate)==4 || face_from_coordinates(coordinate)==5)
                    {
                      if (coordinate[0]-isize<0. || coordinate[0]+isize>(domain.width))
                        {DUNE_THROW(Dune::Exception,"z inlet in x-direction with index " << n << " is invalid (out of the domain)");}
                      if (coordinate[dim-1]-isize<0. || coordinate[dim-1]+isize>(domain.depth))
                        {DUNE_THROW(Dune::Exception,"z inlet in y-direction with index " << n << " is invalid (out of the domain)");}
                    }
              }
          }


        if (cname == "InletsDefault")
          cname = "twophase";
      }

      //! boundary mapper
      /**
       * Returns the boundary mapper.  If there is no boundary mapper yet,
       * initialize it first.
       */
      virtual const Dune::shared_ptr<BoundaryMapper> &boundaryMapper(const int codim) {
        if(!boundaryMapper_) {
          typedef typename GV::template Codim<0>::Iterator EIterator;
          typedef typename GV::IntersectionIterator IIterator;
          //typedef typename GV::Traits::Grid G;
          //typedef typename G::Traits::LocalIdSet IdS;
          // typedef typename G::Traits::LocalIdSet::IdType IdType;

          boundaryMapper_.reset(new BoundaryMapper(gv.grid()));
          const EIterator &eend = gv.template end<0>();
          for(EIterator eit = gv.template begin<0>(); eit != eend; ++eit) {

            if (codim==0)
              // for a universal mapper, this will create a new map entry
              boundaryMapper_->index(eit);

            const IIterator &iend = gv.iend(*eit);
            for(IIterator iit = gv.ibegin(*eit); iit != iend; ++iit)
              if(iit->boundary() && !iit->neighbor())
                {
                  if (codim==1)
                  // for a universal mapper, this will create a new map entry
                  boundaryMapper_->subIndex(eit, iit->indexInInside(), codim);
                }
          }
        }

        return boundaryMapper_;
      }

      template<class T>
      inline void setTime(T time_)
      {
        //std::cout << "setTime to " << time_ <<  std::endl;
        time = time_;

        if ( (time<inlets_iterator->start) || (time>inlets_iterator->end))
          {
            updated = false;

            while(time < inlets_iterator->start){
              if(inlets_iterator == inlets_list.begin())
                {DUNE_THROW(Dune::Exception,"Time is out of any given interval");}
              --inlets_iterator;
              if(verbosity_level)
                std::cout << "update inlets span to (" << inlets_iterator->start << "-" <<  inlets_iterator->end << ") "<< std::endl;
              update();
            }

            while(time > inlets_iterator->end){
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
                  if(verbosity_level)
                    std::cout << "update inlets span to (" << inlets_iterator->start << "-" <<  inlets_iterator->end << ") "<< std::endl;
                  update();
                }
            }
          }

      }



      //! visualize boundary condition to vtk file
      void visualize(std::string name = "")
      {
        /*

        std::string bname;
        bname = cname + "_" + param.get<std::string>("VTKname").c_str() + "_r" + std::to_string(static_cast<long long int>(domain.refine)) + "_t" + std::to_string(static_cast<long long int>(time));

        std::string prefix = static_cast<std::string> (bname);
        if(prefix != "") {

          if (codim==1)
            {
              typedef P0BoundaryGridFunction<GV, RF,1,BoundaryMapper,std::vector<RF> > BTFunction;
              BTFunction btFunction(gv, *this->boundaryMapper(codim), fluxvector);

          if (verbosity_level)
            std::cout  << "Writing boundary condition type to vtk-file " << prefix
                       << std::endl;
          Dune::VTK::NonConformingBoundaryWriter<GV> writer(gv);
          writer.addCellData(new Dune::PDELab::VTKBoundaryGridFunctionAdapter<BTFunction> (btFunction), name);

          std::string path = "boundary";
          mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
          Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw;

          writer.pwrite(prefix,path,"",outputtype);
          if (verbosity_level)
            std::cout << "Writing boundary condition type to vtk-files... done"
                      << std::endl;
            }
          else
        {
          typedef P0GridFunction<GV, RF,dim,BoundaryMapper,std::vector<RF> > BTFunction;
          BTFunction btFunction(gv, *this->boundaryMapper(codim), fluxvector);

          if (verbosity_level)
            std::cout  << "Writing sources type to vtk-file " << prefix
                       << std::endl;
          Dune::VTKWriter<GV> writer(gv);
          writer.addCellData(new Dune::PDELab::VTKGridFunctionAdapter<BTFunction> (btFunction, name));

          std::string path = "sources";
          mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
          Dune::VTK::OutputType outputtype = Dune::VTK::appendedraw;

          writer.pwrite(prefix,path,"",outputtype);
          if (verbosity_level)
            std::cout << "Writing boundary condition type to vtk-files... done"
                      << std::endl;
            }
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
        return  fluxvector[boundaryMapper_->subIndex(is.inside(),  is.indexInInside(), 1)];

      }

        typename Traits::RangeFieldType
      evaluate (const typename Traits::ElementType& e, const typename Traits::DomainType& x, typename Traits::RangeFieldType time=0.) const
      {

        if (!updated)
          DUNE_THROW(Dune::Exception, "New inlets, update necessary");
        return  fluxvector[boundaryMapper_->map(e)];

      }

      /* go through all inlets, read x,y and z coordinates of inlets
       * and volume of inlets and compute function f at given coordinate x
       */
      template<class X>
      inline int
      face_from_coordinates (const X& x) const
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
        DUNE_THROW(Dune::Exception,"coordinates in face detection" << x << "are not at the boundary of the domain or you do not use cube");
        return -1;
      }

      const Dune::ParameterTree & param;

    protected:
      SpanInletsList inlets_list;
      typename SpanInletsList::const_iterator inlets_iterator;


      const GV& gv;
      const HeleshawDomain<dim> domain;
      std::string cname;
      std::string subintervals;
      const int codim;
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


      TwophaseInlets(const GV& gv_,  const Dune::ParameterTree & param_):
        Base(gv_, param_, "InletsDefault", "intervals"),
        type(param_.get<std::string>("InletsDefault.type", "default")),
        fluxtype(param_.get<std::string>("InletsDefault.fluxtype","none")),
        sigma(param_.get<RF>("InletsDefault.sigma",-1.)),
        bump_k(param_.get<RF>("InletsDefault.bumpk",-1.)),
        bump_n(param_.get<RF>("InletsDefault.bumpn",-1.)),
        height_watertable( param_.get<RF>("Setup.waterheight")),// initial hight of water
        porosity(param_.sub(param_.sub("Setup").get<std::string>("material")).get<RF>("porosity")), // medium porosity
        waterheightcorrection(param_.get<bool>("InletsDefault.waterheightcorrection",false)),
        inlets_volume(0.)
      {
        if (gv_.comm().rank()>0)
          verbosity_level=0;

        if (verbosity_level)
          std::cout <<"water height control " << waterheightcorrection << std::endl;

        if (gv_.comm().rank()==0)
          {
            if (fluxtype != "darcy" && fluxtype!="flux")
              {
                DUNE_THROW(Dune::Exception, "Twophase inlets fluxtype must be darcy [m/s] or flux [m^3/s]");
              }
            else
              {
                if (fluxtype=="darcy")
                  std::cout << "Twophase inlets: units are in Darcy velocity [m/s]" << std::endl;
                else
                  std::cout << "Twophase inlets: units are in flux [m^3/s]" << std::endl;
              }
          }

        update();
        if (verbosity_level)
          std::cout << "inlets for twophase constructed" << std::endl;

        control();
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
        size_t face_id = this->face_from_coordinates(x);
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

        update_inlets();
        watch.reset();

        std::vector<RF> flux = inlets_iterator->flux;



        for(typename std::vector<RF>::size_type n = 0; n < inlets_iterator->flux.size(); ++n)
          {
            std::cout << "n " << n << " flux " << flux[n] << " inletsflux " << inlets_iterator->flux[n] << std::endl;
            // flux[n]/=inlets_volume[n];
            //  if (dim==2)
            //    flux[n]/=domain.depth;
            //   RF flux = inlets_iterator->flux[n];
            //  std::cout << "\nDarcy flux velocity " << flux/inlets_volume[n] << " m/s = pore velocity " << flux/porosity/inlets_volume[n]  << " m/s = " << flux/porosity*86400/inlets_volume[n]  << " m/d" ;

          }
        std::vector<RF> scale =  scale_factor();
        // if (gv.comm().size()>1)
        //   for (size_t i=0;i<scale.size();i++)
        //    scale[i] =  gv.comm().sum(scale[i]);


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
                //  RF face_volume = ig.geometry().volume();

                Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
                  global = ig.intersection().geometry().global(face_local);
                size_t id = boundaryMapper_->subIndex(ig.inside(), ig.indexInInside(), 1);
                size_t fid = this->face_from_coordinates(global);

                RF faceintegral=average_at_intersection(ig)*scale[fid]*flux[fid]/(domain.depth*inlets_volume[fid]);
                if (dim==3)
                  faceintegral=average_at_intersection(ig)*scale[fid]*flux[fid]/(inlets_volume[fid]);
                if (inlets_volume[fid]<1.e-10)
                  faceintegral = 0;

                if (fluxtype=="darcy")
                  faceintegral = average_at_intersection(ig)*scale[fid]*flux[fid];


                //            std::cout << "side " << face_from_coordinates(global) << " value " << faceintegral << " " <<  average_at_intersection(ig) << " " << flux[this->face_from_coordinates(global)]<< std::endl;

                fluxvector[id] -= faceintegral;
              }
          }


        updated=true;
        if (verbosity_level)
          std::cout << "inlets timespan update ... done : " << watch.elapsed() << " s" << std::endl;

        if (verbosity_visualize){
          if (verbosity_level)
            std::cout << "visualize new boundary conditions " << std::endl;
          std::string name;

          this->visualize("darcy velocity in m/s");
        }

      }

    private:

      void control()
      {

        std::vector<RF> inflow_control;
        inflow_control.resize(4);
        if (dim>2)
          inflow_control.resize(6);

        std::fill(inflow_control.begin(), inflow_control.end(), 0.);

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
                size_t id = boundaryMapper_->subIndex(ig.inside(), ig.indexInInside(), 1);
                // std::cout << fluxvector[id] << " " << face_volume << " " << domain.depth << std::endl;
                RF lcontrol = fluxvector[id]*1.e6*3600.*face_volume;
                /*
                  in fluxvector is actually the darcy velocity (

                */
                if (dim<3)
                  lcontrol*=domain.depth;
                inflow_control[this->face_from_coordinates(global)]+=lcontrol;
              }
          }

        if (gv.comm().size()>1)
          for (size_t i=0;i<inflow_control.size();i++)
            inflow_control[i]= gv.comm().sum(inflow_control[i]);

        if (gv.comm().rank()==0)
          for (auto it = inflow_control.begin();it!=inflow_control.end();++it)
            std::cout << "Flow in side " << it - inflow_control.begin() << " is " << *it << " ml/h = " << *it/3600./1.e6 << " m^3/s\n";
      }


      //! computes an average value of function f at intersection ig
      template<class IG>
      RF average_at_intersection(const IG& ig) const
      {

        // face geometry
        //const Dune::FieldVector<DF,IG::dimension-1>&
        //face_local = Dune::ReferenceElements<RF,IG::dimension-1>::general(ig.geometry().type()).position(0,0);
        RF face_volume = ig.geometry().volume();

        // Dune::FieldVector<typename Traits::IntersectionType::ctype,Traits::IntersectionType::dimension>
        // global = ig.intersection().geometry().global(face_local);

        RF integral_on_face = 0.;

        Dune::GeometryType gt = ig.geometry().type();
        const int intorder = 50;
        const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gt,intorder);
        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
          {

            //Dune::FieldVector<DF,dim> position = ig.geometry().global(qit->position());
            RF fval = f(ig.geometry().global(qit->position()));
            RF weight = qit->weight();
            RF detjac = ig.geometry().integrationElement(qit->position());
            integral_on_face +=fval*weight*detjac;
          }
        return  (integral_on_face/face_volume);
      }

      //! scales function f that the integral of f over face is equal to volume of face
      std::vector<RF> scale_factor() const
      {
        std::vector<RF> scale;
        scale.resize(4);
        if (dim>2)
          scale.resize(6);
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

                int s = this->face_from_coordinates(global);

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


        //std::cout << "rank " << gv.comm().rank() << " scale " << scale[bottom] << " volume " << inlets_volume[bottom] << std::endl;

        for (size_t i = 0; i<scale.size();i++)
          {
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

          }
        return scale;
      }

      void update_inlets()
      {
        inlets_volume.resize(4);
        if (dim>2)
          inlets_volume.resize(6);

        std::fill(inlets_volume.begin(), inlets_volume.end(), 0.);

        //it is global for all processors!!!
        for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
          {
            Dune::FieldVector<RF,dim> c;
            c[0] = inlets_iterator->inletsx[i];
            c[dim-1] = inlets_iterator->inletsy[i];
            if (dim>2)
              c[1] = inlets_iterator->inletsz[i];
            RF r = inlets_iterator->inletsize[i];
            if (dim>2)
              r = 1./4. * M_PI * r * r;
            inlets_volume[this->face_from_coordinates(c)]+=r;

          }

      }


    private:
      std::string type;
      std::string fluxtype; //darcy[m/s] or flux [m^3/s]
      const RF sigma;
      const RF bump_k;
      const RF bump_n;
      const RF height_watertable;
      const RF porosity;
      const bool waterheightcorrection;
      std::vector<RF> inlets_volume;
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


      TransportInlets(const GV& gv_,  const Dune::ParameterTree & param_, const std::string cname_, const std::string intervals_="intervals"):
        Base(gv_,param_,cname_,intervals_),
        igstorage(), minh(computeH<GV,true>(gv))
      {
        if (gv_.comm().rank()>0)
          verbosity_level=0;

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
            size_t id = boundaryMapper_->subIndex(ig.inside(), ig.indexInInside(), 1);
            fluxvector[id]=f(global);
          }
        updated = true;

        if (verbosity_visualize){
          if (verbosity_level)
            std::cout << "visualize new boundary conditions " << std::endl;
          this->visualize("concentration");
        }
      }

      /* go through all inlets, read x,y and z coordinates of inlets
       * and volume of inlets and compute function f at given coordinate x
       */
      template<class X>
      inline typename Traits::RangeFieldType
      f (const X& x) const
      {
        if (dim==1)
          return inlets_iterator->input[0];

        X c(0);

        for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
          {
            c[0] = inlets_iterator->inletsx[i];
            c[dim-1] = inlets_iterator->inletsy[i];
            if (dim>2)
              c[1] = inlets_iterator->inletsz[i];
            RF r = inlets_iterator->inletsize[i]/2.+minh/2.0;
            c-=x;
            RF d = c.two_norm();

            // std::cout << "inlet input at " << x << " is " << inlets_iterator->input[i] << " with d " << d << " and r " << r << std::endl;
            // on inflow init or input!
            if (d<=r)
              {
                // std::cout << "inlet input at " << x << " is " << inlets_iterator->input[i] << " with r " << r << std::endl;
                return inlets_iterator->input[i];
              }
          }


        //value = inlets_iterator->input[0];
        // if (value<eps)
        return inlets_iterator->initial[0];
        //return value;
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
      RF minh;
      using Base::boundaryMapper_;


    };


   //! for two phases inlets class
    template <class Traits, class RF>
    class WaterSource
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


      WaterSource(const GV& gv_,  const Dune::ParameterTree & param_, const std::string name = "WaterSource"):
        Base(gv_, param_, name, "intervals",0,false,false),
        type(param_.template get<std::string>(name+".type", "default")),
        fluxtype(param_.template get<std::string>(name+".fluxtype","none")),
        sigma(param_.template get<RF>(name+".sigma",-1.)),
        bump_k(param_.template get<RF>(name+".bumpk",-1.)),
        bump_n(param_.template get<RF>(name+".bumpn",-1.)),
        height_watertable( param_.get<RF>("Setup.waterheight")),// initial hight of water
        porosity(param_.sub(param_.sub("Setup").get<std::string>("material")).get<RF>("porosity")), // medium porosity
        waterheightcorrection(param_.template get<bool>(name+".waterheightcorrection",false)),
        onlyinflow(param_.template get<bool>(name+".onlyinflow",true)),
        source_volume(0.)
      {
        if (gv_.comm().rank()>0)
          verbosity_level=0;

        if (verbosity_level)
          std::cout <<"water height control " << waterheightcorrection << std::endl;
        /*
        if (gv_.comm().rank()==0)
          {
            if (fluxtype != "darcy" && fluxtype!="flux")
              {
                DUNE_THROW(Dune::Exception, "Twophase sources fluxtype must be darcy [m/s] or flux [m^3/s]");
              }
            else
              {
                if (fluxtype=="darcy")
                  std::cout << "Twophase sources: units are in Darcy velocity [m/s]" << std::endl;
                else
                  std::cout << "Twophase sources: units are in flux [m^3/s]" << std::endl;
              }
          }
        */

        update();
        if (verbosity_level)
          std::cout << "sources for twophase constructed" << std::endl;

        //control();
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

            /* go through all inlets, read x,y and z coordinates of inlets
       * and volume of inlets and compute function f at given coordinate x
       */
      template<class X>
      inline typename Traits::RangeFieldType
      fluxvalue (const X& x) const
      {
         X c(0);

        for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
          {
            c[0] = inlets_iterator->inletsx[i];
            c[dim-1] = inlets_iterator->inletsy[i];
            if (dim>2)
              c[1] = inlets_iterator->inletsz[i];
            RF r = inlets_iterator->inletsize[i]/2.;
            c-=x;
            RF d = c.two_norm();
            if (d<=r)
              {
                return inlets_iterator->input[i];
              }
          }
        return 0.0;
      }


      template<class X>
      inline typename Traits::RangeFieldType
      water_correction (const X& x) const
      {
        if (!waterheightcorrection)
          return 1.;

        RF v1 = (patm - height_watertable*Rho_l*9.81*porosity)/patm;
        if (x[dim-1] > height_watertable)
          return 1;
        else
          return (1.-v1)/height_watertable*x[dim-1]+v1;
        return 1.;
      }



      //! update the construction values for a given timestep
      void update()
      {
        if (updated) {
          if (verbosity_level)
            std::cout << "sources update not necessary" << std::endl;
          return;
        }

        if (verbosity_level){
          std::cout << "updating sources";
          if (gv.comm().size()>1)
            std::cout << " on processor " << gv.comm().rank();
          std::cout << std::endl;
        }

        update_inlets();
        watch.reset();

        std::vector<RF> flux = inlets_iterator->input;

        for(typename std::vector<RF>::size_type n = 0; n < inlets_iterator->flux.size(); ++n)
          {
            std::cout << "n " << n << " flux " << flux[n]<< std::endl;
          }
        RF scale =  scale_factor();

        if (verbosity_level)
          {
            if (gv.comm().size()>1 && verbosity_level>1)
              std::cout << "scale on processor " << gv.comm().rank() << std::endl;
            std::cout << "scale is " << scale << " source volume : " << source_volume <<  std::endl;
          }

        std::fill(fluxvector.begin(), fluxvector.end(), 0.);

        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename Dune::PDELab::ElementGeometry<Element> EG;
        // loop over cels
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {


            EG eg(*it);

            Dune::FieldVector<DF,dim>
              global = eg.geometry().center();
            size_t id = boundaryMapper_->map(*it);

            RF elementintegral=average_at_intersection(eg)*scale*fluxvalue(global)/(domain.depth*source_volume);
            //std::cout << "average " << average_at_intersection(eg) << " scale " << scale << " fluxvalue " << fluxvalue(global) << std::endl;
            if (dim==3)
              elementintegral=average_at_intersection(eg)*scale*fluxvalue(global)/(source_volume);
            if (source_volume<1.e-14)
              elementintegral = 0;

            fluxvector[id] -= elementintegral;
          }

        updated=true;
        if (verbosity_level)
          std::cout << "sources timespan update ... done : " << watch.elapsed() << " s" << std::endl;

        if (verbosity_visualize){
          if (verbosity_level)
            std::cout << "visualize new boundary conditions " << std::endl;
          std::string name;

          this->visualize("darcy velocity in m/s");
        }

        control();
      }

    private:


      void control()
      {

        RF inflow;
        RF flow;

        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename Dune::PDELab::ElementGeometry<Element> EG;

        // loop over cels
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {
            EG eg(*it);
            size_t id = boundaryMapper_->map(*it);

            RF volume = eg.geometry().volume();
            RF lcontrol = fluxvector[id]*1.e6*3600.*volume;
            if (dim<3)
              lcontrol = fluxvector[id]*1.e6*3600.*volume*domain.depth;
            if (lcontrol>0)
              inflow+=lcontrol;
            flow+=lcontrol;

          }

        if (gv.comm().size()>1){
          inflow= gv.comm().sum(inflow);
          flow= gv.comm().sum(flow);
        }

        if (gv.comm().rank()==0)
          {
          std::cout << "Flow in is " << " is " << flow << " ml/h = " << flow/3600./1.e6 << " m^3/s\n";
          std::cout << "Inflow in is " << " is " << inflow << " ml/h = " << inflow/3600./1.e6 << " m^3/s\n";
          }
      }


      //! computes an average value of function f element ig
      template<class IG>
      RF average_at_intersection(const IG& ig) const
      {

        RF face_volume = ig.geometry().volume();
        RF integral_on_face = 0.;
        Dune::GeometryType gt = ig.geometry().type();
        const int intorder = 50;
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
          {
            RF fval = f(ig.geometry().global(qit->position()));
            RF weight = qit->weight();
            RF detjac = ig.geometry().integrationElement(qit->position());
            integral_on_face +=fval*weight*detjac;
          }
        return  (integral_on_face/face_volume);
      }

      //! scales function f that the integral of f over face is equal to volume of face
      RF scale_factor() const
      {
        RF scale(0.);
        // RF volume(0.0);

        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename Dune::PDELab::ElementGeometry<Element> EG;


        // loop over cels
        for (ElementIterator it = gv.template begin<0>(); it!=gv.template end<0>(); ++it)
          {
            if (it->partitionType()==Dune::OverlapEntity)
              continue;

            EG eg(*it);
            Dune::GeometryType gt = eg.geometry().type();
            const int intorder = 10;
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);
            // loop over quadrature points
            for (typename Dune::QuadratureRule<DF,dim>::const_iterator qit=rule.begin(); qit!=rule.end(); ++qit)
              {

                Dune::FieldVector<DF,dim> position = eg.geometry().global(qit->position());

                RF fval = f(position);
                if (!onlyinflow && fluxvalue(position)>0)
                  fval=0;
                RF weight = qit->weight();
                RF detjac = eg.geometry().integrationElement(qit->position());
                scale +=fval*weight*detjac;
              }
          }


        if (gv.comm().size()>1)
          {
            scale=gv.comm().sum(scale);
          }

        if (std::abs(scale)<eps)
          scale=0;
        else
          scale=(source_volume/scale);


        return scale;
      }

      void update_inlets()
      {
        source_volume=0.0;
        //it is global for all processors!!!
        for (size_t i = 0; i<inlets_iterator->inletsx.size();i++)
          {
            Dune::FieldVector<RF,dim> c;
            c[0] = inlets_iterator->inletsx[i];
            c[dim-1] = inlets_iterator->inletsy[i];
            if (dim>2)
              c[1] = inlets_iterator->inletsz[i];
            RF r = inlets_iterator->inletsize[i];

            if (dim==2)
              r = 1./4. * M_PI * r * r;

            if (onlyinflow)
              source_volume+=r;
            else if ((inlets_iterator->input[i])<0)
              source_volume+=r;
          }
      }


    private:
      std::string type;
      std::string fluxtype; //darcy[m/s] or flux [m^3/s]
      const RF sigma;
      const RF bump_k;
      const RF bump_n;
      const RF height_watertable;
      const RF porosity;
      const bool waterheightcorrection;
      const bool onlyinflow;
      RF source_volume;
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



  } // end namespace Dycap
} // end namespace Dune


#endif
