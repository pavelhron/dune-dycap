// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_DYCAP_TWOPHASEVELOCITY_HH
#define DUNE_DYCAP_TWOPHASEVELOCITY_HH

#include<tgmath.h>
#include<src/utilities/rt0qfem.hh>
#include<dune/localfunctions/raviartthomas/raviartthomascube.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>

//! enum to choose between gas and liquid phase
enum TwoPhaseVelocityName
  {
    TwoPhaseVelLiquid = 0,
    TwoPhaseVelCapillary = 1
  };

template<typename TP, TwoPhaseVelocityName P>
class TwoPhaseParameterMapper;

template<typename TP>
struct TwoPhaseParameterMapper<TP, TwoPhaseVelLiquid>
{
  typedef typename TP::Traits Traits;
  TwoPhaseParameterMapper(const TP& tp) : tp_(tp) {}

  typename Traits::RangeFieldType
  s (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType pc) const
  {
    return tp_.s_l(e,x,pc);
  }

  typename Traits::RangeFieldType
  kr (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType s_l) const
  {
    return tp_.kr_l(e,x,s_l);
  }

  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType p_l) const
  {
    return tp_.mu_l(e,x,p_l);
  }

  typename Traits::RangeFieldType
  nu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType p_l) const
  {
    return tp_.nu_l(e,x,p_l);
  }

  typename Traits::RangeFieldType
  rho (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType p_l) const
  {
    return tp_.rho_l(e,x,p_l);
  }

  int
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.bc_l(is,x,time);
  }

  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.g_l(is,x,time);
  }

  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.j_l(is,x,time);
  }

private:
  const TP& tp_;
};

template<typename TP>
struct TwoPhaseParameterMapper<TP, TwoPhaseVelCapillary>
{
  typedef typename TP::Traits Traits;
  TwoPhaseParameterMapper(const TP& tp) : tp_(tp) {}

  typename Traits::RangeFieldType
  s (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
     typename Traits::RangeFieldType pc) const
  {
    return 1.0 - tp_.s_l(e,x,pc);
  }

  typename Traits::RangeFieldType
  kr (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType s_g) const
  {
    return tp_.kr_g(e,x,s_g);
  }

  typename Traits::RangeFieldType
  mu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType p_g) const
  {
    return tp_.mu_g(e,x,p_g);
  }

  typename Traits::RangeFieldType
  nu (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
      typename Traits::RangeFieldType p_g) const
  {
    return tp_.nu_g(e,x,p_g);
  }

  typename Traits::RangeFieldType
  rho (const typename Traits::ElementType& e, const typename Traits::DomainType& x,
       typename Traits::RangeFieldType p_g) const
  {
    return tp_.rho_g(e,x,p_g);
  }

  int
  bc (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.bc_g(is,x,time);
  }

  typename Traits::RangeFieldType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.g_g(is,x,time);
  }

  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, typename Traits::RangeFieldType time) const
  {
    return tp_.j_g(is,x,time);
  }

private:
  const TP& tp_;
};

/** \brief Provide velocity/flux field for phase P

    Uses RT0 interpolation on a cell.

    @tparam T  provides TwoPhaseParameterInterface
    @tparam PL P0 function for liquid phase pressure
    @tparam PG P0 function for gas phase pressure
    @tparam P  TwoPhaseVelocityName enum selecting the phase P
*/
template<typename  TP, typename PL, typename PC, TwoPhaseVelocityName P>
class TwoPhasePhaseVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,
                                                                           PL::Traits::GridViewType::dimension,
                                                                           Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                          TwoPhasePhaseVelocity<TP,PL,PC,P> >
{
  // extract useful types
  typedef typename PL::Traits::GridViewType GV;
  typedef typename GV::Grid::ctype DF;
  typedef typename PL::Traits::RangeFieldType RF;
  typedef typename PL::Traits::RangeType RangeType;
  enum { dim = PL::Traits::GridViewType::dimension };
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::IntersectionIterator IntersectionIterator;
  typedef typename IntersectionIterator::Intersection Intersection;

  const TP& tp;
  const TwoPhaseParameterMapper<TP, P> tpm;
  const PL& pl;
  const PC& pc;
  Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0> rt0fe;
  typename TP::Traits::RangeFieldType time;
  bool compFlux;

  typedef typename Dune::RaviartThomasCubeLocalFiniteElement<DF,RF,dim,0>::Traits::LocalBasisType::Traits::RangeType RT0RangeType;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

  typedef Dune::PDELab::GridFunctionBase<Traits,TwoPhasePhaseVelocity<TP,PL,PC,P> > BaseT;

  TwoPhasePhaseVelocity (const TP& tp_, const PL& pl_, const PC& pc_, bool compFlux_ = true) :
    tp(tp_), tpm(tp), pl(pl_), pc(pc_), time(0), compFlux(compFlux_) {}

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (typename TP::Traits::RangeFieldType time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    // cell geometry
    const Dune::FieldVector<DF,dim>&
      inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
      general(e.type()).position(0,0);
    Dune::FieldVector<DF,dim>
      inside_cell_center_global = e.geometry().global(inside_cell_center_local);


    // absolute permeability
    RF k_abs_inside = tp.k_abs(e,inside_cell_center_local);

    // pressure evaluation
    typename PL::Traits::RangeType pl_inside, pc_inside, p_inside;
    pl.evaluate(e,inside_cell_center_local,pl_inside);
    pc.evaluate(e,inside_cell_center_local,pc_inside);
    if (P == TwoPhaseVelLiquid)
      p_inside = pl_inside;
    else
      p_inside = pl_inside+pc_inside;

    // density evaluation
    RF nu_inside = tpm.nu(e,inside_cell_center_local,p_inside);

    // for coefficient computation
    RF vn[2*dim];    // normal velocities
    RF coeff[2*dim]; // RT0 coefficient
    Dune::FieldMatrix<typename Traits::DomainFieldType,dim,dim>
      B = e.geometry().jacobianInverseTransposed(x); // the transformation. Assume it is linear
    RF determinant = B.determinant();

    // loop over cell neighbors
    IntersectionIterator endit = pl.getGridView().iend(e);
    for (IntersectionIterator iit = pl.getGridView().ibegin(e); iit!=endit; ++iit)
      {
        // set to zero for processor boundary
        vn[iit->indexInInside()] = 0.0;

        // face geometry
        const Dune::FieldVector<DF,dim-1>&
          face_local = Dune::ReferenceElements<DF,dim-1>::general(iit->geometry().type()).position(0,0);


        // interior face
        if (iit->neighbor())
          {
            const Dune::FieldVector<DF,dim>&
              outside_cell_center_local = Dune::ReferenceElements<DF,dim>::
              general(iit->outside().type()).position(0,0);
            Dune::FieldVector<DF,dim>
              outside_cell_center_global = iit->outside().geometry().global(outside_cell_center_local);

            // distance of cell centers
            Dune::FieldVector<DF,dim> d(outside_cell_center_global);
            d -= inside_cell_center_global;
            RF distance = d.two_norm();

            // absolute permeability
            RF k_abs_outside = tp.k_abs(iit->outside(),outside_cell_center_local);

            // gravity times normal
            RF gn = tp.gravity()*iit->unitOuterNormal(face_local);

            // pressure evaluation
            typename PL::Traits::RangeType pl_outside, pc_outside, p_outside;
            pl.evaluate((iit->outside()),outside_cell_center_local,pl_outside);
            pc.evaluate((iit->outside()),outside_cell_center_local,pc_outside);
            if (P == TwoPhaseVelLiquid)
              p_outside = pl_outside;
            else
              p_outside = pc_outside+pl_outside;

            // needed for both phases
            RF pc_upwind, s_upwind;

            // P phase calculation
            RF rho_inside = tpm.rho(e,inside_cell_center_local,p_inside);
            RF rho_outside = tpm.rho(iit->outside(),outside_cell_center_local,p_outside);
            RF w = (p_inside-p_outside)/distance + aavg(rho_inside,rho_outside)*gn; // determines direction
            if (w>=0) // upwind capillary pressure on face
              {
                pc_upwind = pc_inside;
                s_upwind = tpm.s(e,inside_cell_center_local,pc_upwind);
              }
            else
              {
                pc_upwind = pc_outside;
                s_upwind = tpm.s(iit->outside(),outside_cell_center_local,pc_upwind);
              }
            RF lambda_inside = tpm.kr(iit->inside(),inside_cell_center_local,s_upwind)/
              tpm.mu(iit->inside(),inside_cell_center_local,p_inside);
            RF lambda_outside = tpm.kr(iit->outside(),outside_cell_center_local,s_upwind)/
              tpm.mu(iit->outside(),outside_cell_center_local,p_outside);
            RF sigma = havg(lambda_inside*k_abs_inside,lambda_outside*k_abs_outside);

            // density
            RF nu_outside = tpm.nu(iit->outside(),outside_cell_center_local,p_outside);

            // set coefficient
            if (compFlux)
              vn[iit->indexInInside()] = aavg(nu_inside,nu_outside) * sigma * w;
            else
              vn[iit->indexInInside()] = sigma * w;
          }

        // boundary face
        else if (iit->boundary())
          {
            // distance of cell center to boundary
            Dune::FieldVector<DF,dim> d = iit->geometry().global(face_local);
            d -= inside_cell_center_global;
            RF distance = d.two_norm();

            // evaluate boundary condition type
            int bc = tpm.bc(*iit,face_local,time);

            // phase Dirichlet boundary
            if (bc==1)
              {
                RF rho_inside = tpm.rho(e,inside_cell_center_local,p_inside);
                RF g = tpm.g(*iit,face_local,time);
                RF gn = tp.gravity()*iit->unitOuterNormal(face_local);
                RF w = (p_inside-g)/distance + rho_inside*gn;
                RF s = tpm.s(e,inside_cell_center_local,pc_inside);
                RF lambda_inside = tpm.kr(e,inside_cell_center_local,s)/
                  tpm.mu(e,inside_cell_center_local,p_inside);
                RF sigma = lambda_inside*k_abs_inside;

                if (compFlux)
                  vn[iit->indexInInside()] = nu_inside * sigma * w;
                else
                vn[iit->indexInInside()] = sigma * w;
              }

            // phase Neumann boundary
            if (bc==0)
              {
                RF j = tpm.j(*iit,face_local,time);
                if (compFlux)
                  vn[iit->indexInInside()] = j;
                else
                  vn[iit->indexInInside()] = j;// / nu_inside;
              }
          }

        // compute coefficient
        Dune::FieldVector<DF,dim> vstar=iit->unitOuterNormal(face_local); // normal on tranformef element
        vstar *= vn[iit->indexInInside()];
        Dune::FieldVector<RF,dim> normalhat(0); // normal on reference element
        if (iit->indexInInside()%2==0)
          normalhat[iit->indexInInside()/2] = -1.0;
        else
          normalhat[iit->indexInInside()/2] =  1.0;
        Dune::FieldVector<DF,dim> vstarhat(0);
        B.umtv(vstar,vstarhat); // Piola backward transformation
        vstarhat *= determinant;
        coeff[iit->indexInInside()] = vstarhat*normalhat;
      }

    // compute velocity on reference element
    std::vector<RT0RangeType> rt0vectors(rt0fe.localBasis().size());
    rt0fe.localBasis().evaluateFunction(x,rt0vectors);
    typename Traits::RangeType yhat(0);
    for (unsigned int i=0; i<rt0fe.localBasis().size(); i++)
      yhat.axpy(coeff[i],rt0vectors[i]);

    // apply Piola transformation
    B.invert();
    y = 0;
    B.umtv(yhat,y);
    y /= determinant;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }

private:

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

/** \brief Provide velocity/flux field for liquid phase

    Uses RT0 interpolation on a cell.

    - T  : provides TwoPhaseParameterInterface
    - PL : P0 function for liquid phase pressure
    - PG : P0 function for gas phase pressure
*/
template<typename  TP, typename PL, typename PC>
class VelocityLiquid
  : public TwoPhasePhaseVelocity<TP,PL,PC,TwoPhaseVelLiquid>
{
public:
  VelocityLiquid (const TP& tp_, const PL& pl_, const PC& pc_, bool compFlux_ = true) :
    TwoPhasePhaseVelocity<TP,PL,PC,TwoPhaseVelLiquid> (tp_, pl_, pc_, compFlux_) {}
};

/** \brief Provide velocity/flux field for gas phase

    Uses RT0 interpolation on a cell.

    - T  : provides TwoPhaseParameterInterface
    - PL : P0 function for liquid phase pressure
    - PC : P0 function for gas phase pressure
*/
template<typename  TP, typename PL, typename PC>
class VelocityGas
  : public TwoPhasePhaseVelocity<TP,PL,PC,TwoPhaseVelCapillary>
{
public:
  VelocityGas (const TP& tp_, const PL& pl_, const PC& pc_, bool compFlux_ = true) :
    TwoPhasePhaseVelocity<TP,PL,PC,TwoPhaseVelCapillary> (tp_, pl_, pc_, compFlux_) {}
};



template<typename PL>
class ZeroVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,
                                                                           PL::Traits::GridViewType::dimension,
                                                                           Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                          ZeroVelocity<PL> >
{
  // extract useful types
  typedef typename PL::Traits::GridViewType GV;
  typedef typename GV::Grid::ctype DF;
  typedef typename PL::Traits::RangeFieldType RF;
  typedef typename PL::Traits::RangeType RangeType;
  enum { dim = PL::Traits::GridViewType::dimension };

  const PL& pl;
  RF time;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

  ZeroVelocity ( const PL& pl_) :
    pl(pl_),time(0){}

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (RF time_)
  {
    time = time_;
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y=0.0;;
  }

  inline const typename Traits::GridViewType& getGridView ()
  {
    return pl.getGridView();
  }


};


/**
 * \brief grid function to compute unit velocity in m/s
 *
 * \tparam TP    TwoPhaseParameters
 */
template<typename  TP, typename PL, typename SL>
class UnitVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename PL::Traits::GridViewType,
                                                                           typename PL::Traits::RangeFieldType,
                                                                           PL::Traits::GridViewType::dimension,
                                                                           Dune::FieldVector<typename PL::Traits::RangeFieldType,PL::Traits::GridViewType::dimension> >,
                                          UnitVelocity<TP,PL,SL> >
{
  // extract useful types
  typedef typename PL::Traits::GridViewType GV;
  typedef typename GV::Grid::ctype DF;
  typedef typename PL::Traits::RangeFieldType RF;
  typedef typename PL::Traits::RangeType RangeType;
  enum { dim = PL::Traits::GridViewType::dimension };

  const TP& tp;
  SL& sl;
  RF value;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

  UnitVelocity (const TP& tp_, SL& sl_, RF value_) :
    tp(tp_), sl(sl_), value(value_) {}

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (typename TP::Traits::RangeFieldType time_)
  {
    time = time_;
  }

  //! return constant value over the domain multiplied with saturation
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    const Dune::FieldVector<DF,dim>&
      inside_cell_center_local = Dune::ReferenceElements<DF,dim>::
      general(e.type()).position(0,0);
    Dune::FieldVector<DF,dim>
      inside_cell_center_global = e.geometry().global(inside_cell_center_local);


    y=0.;
    // set coefficient
    //if (inside_cell_center_global[dim-1]<0.2)
    typename PL::Traits::RangeType slvalue;
    sl.evaluate(e,inside_cell_center_global,slvalue);

    y[0] = value*slvalue;
  }

  //! get GridView
  inline const typename Traits::GridViewType& getGridView ()
  {
    return tp.getGridView();
  }

};


/**
 * \brief grid function to compute constant velocity in m/s
 *
 * \tparam S
 */
template<typename GV, typename RF>
class ConstantVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::dimension,Dune::FieldVector<RF,GV::dimension> >, ConstantVelocity<GV,RF> >
{

  enum { dim = GV::dimension };

  GV gv;
  RF xvalue,yvalue,zvalue;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

  ConstantVelocity (const GV& gv_, RF xvalue_, RF yvalue_, RF zvalue_ = 0.0) :
    gv(gv_), xvalue(xvalue_), yvalue(yvalue_), zvalue(zvalue_) {}

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (typename Traits::RangeFieldType time_)
  {
  }

  //! return constant value over the domain multiplied with saturation
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    y[0] = xvalue;
    y[dim-1] = yvalue;
    if (dim>2)
      y[1] = zvalue;
  }

  //! get GridView
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }

};


/**
 * \brief grid function to compute constant velocity in m/s
 *
 * \tparam S
 */
template<typename GV, typename RF>
class RotateVelocity
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,RF,GV::dimension,Dune::FieldVector<RF,GV::dimension> >, RotateVelocity<GV,RF> >
{

  enum { dim = GV::dimension };

  GV gv;
  RF xvalue,yvalue,zvalue;

public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,dim,Dune::FieldVector<RF,dim> > Traits;

  RotateVelocity (const GV& gv_, RF xvalue_, RF yvalue_, RF zvalue_ = 0.0) :
    gv(gv_), xvalue(xvalue_), yvalue(yvalue_), zvalue(zvalue_) {}

  // set time where operator is to be evaluated (i.e. end of the time intervall)
  void set_time (typename Traits::RangeFieldType time_)
  {
  }

  //! return constant value over the domain multiplied with saturation
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    Dune::FieldVector<RF,dim> global = e.geometry().global(x);
    global-=0.5;
    RF r = global.two_norm();
    RF theta = std::atan2(global[1],global[0]);//in grades
    //RF theta1 = std::atan(global[0]/global[1]);//in grades
    y[0] = -r*std::sin(theta);
    y[dim-1] = r*std::cos(theta);
    /*Dune::FieldVector<RF,dim> global = e.geometry().global(x);
    global-=0.5;
    RF r = global.two_norm();
    RF theta = std::atan(global[1]/global[0]);//in grades
    y[0] = -r*std::sin(theta);
    y[dim-1] = r*std::cos (theta);*/
  }

  //! get GridView
  inline const typename Traits::GridViewType& getGridView ()
  {
    return gv;
  }

};


/**
 * \brief grid function to visualize pore velocity
 if darcy velocity is in m/s, pore velocity is in m/d
 * \tparam VEL        darcy velocity
 */
template<typename  VEL, typename TP>
class PoreVelocity
  : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<typename VEL::Traits::GridViewType,
                                   typename VEL::Traits::RangeFieldType,
                                   VEL::Traits::dimRange,
                                   typename VEL::Traits::RangeType>,
  PoreVelocity<VEL,TP> >
{

  const VEL& vel;
  const TP& tp;

public:

  //! exort the GFS Traits
  typedef Dune::PDELab::GridFunctionTraits<typename VEL::Traits::GridViewType,
                                           typename VEL::Traits::RangeFieldType,
                                           VEL::Traits::dimRange,
                                           typename VEL::Traits::RangeType> Traits;
  //! \brief constructor for PoreVelocity
  /**
   * \tparam VEL  TwophaseVelocity
   * \tparam TP   TwophaseParameters, porosity is needed
   */
  PoreVelocity (const VEL& vel_, const TP& tp_) : vel(vel_), tp(tp_) {}

  //! return pore velocity in [m/d]
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {
    typename Traits::RangeType velo;
    vel.evaluate(e,x, velo);
    velo*=(86400./tp.phi(e,x));
    y = velo;
  }

  //! get GriView
  inline const typename Traits::GridViewType& getGridView ()
  {
    return vel.getGridView();
  }
};


#endif
