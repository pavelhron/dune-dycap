// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
#ifndef DUNE_PM_ODE_HH
#define DUNE_PM_ODE_HH

#include <iostream>
#include<iomanip>
#include<vector>
#include<map>
#include<memory>
#include<cmath>
#include<dune/common/exceptions.hh>




/** @brief Explicit Euler method as an example for an ODE solver

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/

template<class M>
class OdeBase
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

  typedef typename M::V V;

  typedef typename M::Matrix Matrix;

  //! constructor stores reference to the model
  OdeBase(const M& model_)
    : model(model_), u(model.size()), f(model.size())
  {
    model.initialize(t,u);
    dt = 0.1;
    dt_min = 1.e-12;
    TOL = time_type(0.0001);
  }

  virtual void step() = 0;

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dt_;
  }


  //! set current state
  void set_state (time_type t_, const V & u_)
  {
    t = t_;
    u = u_;
  }

  //! get current state
  const V & get_state () const
  {
    return u;
  }

  //! get current time
  time_type get_time () const
  {
    return t;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_dt () const
  {
    return dt;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_TOL () const
  {
    return TOL;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_dtmin () const
  {
    return dt_min;
  }

  //! set tolerance for adaptive computation
  void set_TOL (time_type TOL_)
  {
    TOL = TOL_;
  }

  //! set tolerance for adaptive computation
  void set_dt_min (time_type dt_min_)
  {
    dt_min = dt_min_;
  }


  //! return consistency order of the method
  virtual size_type get_order () const = 0;

  //! return consistency order of the method
  virtual bool adaptive () const = 0;


protected:
  const M& model;
  time_type t, dt, dt_min, TOL;
  V u;
  V f;
};

/** @brief A dummy ode base class

    \tparam M the model type
*/

template<class M>
class DummyOdeBase : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;

  /** \brief export size_type */
  typedef typename M::size_type size_type;

  //! constructor stores reference to the model
  DummyOdeBase(const M& model)
    : Base(model)
  {
  }

  virtual void step()
  {DUNE_THROW(Dune::NotImplemented,"Illegal member call of DummyOdeBase class!");}

  //! return consistency order of the method
  virtual size_type get_order () const
  {DUNE_THROW(Dune::NotImplemented,"Illegal member call of DummyOdeBase class!");}

  //! return consistency order of the method
  virtual bool adaptive () const
  {DUNE_THROW(Dune::NotImplemented,"Illegal member call of DummyOdeBase class!");}


};




template<class M>
class ExplicitEuler : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  ExplicitEuler(const M& model)
    : Base(model)
  {
  }

  //! do one step
  void step ()
  {
    model.f(t,u,f);   // evaluate model
    u.axpy(dt,f);   // advance state
    t += dt;          // advance time
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 1;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::t;
  using Base::dt;
  using Base::u;
  using Base::f;
};


/** @brief Modified Euler method (order 2 with 2 stages)

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class ModifiedEuler : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;

  //! constructor stores reference to the model
  ModifiedEuler (const M& model )
    : Base(model), w(model.size()), k1(model.size()), k2(model.size())
  {
    c2 = 0.5;
    a21 = 0.5;
    b2 = 1.0;
    model.initialize(t,u);
    dt = 0.1;
  }

  //! do one step
  void step ()
  {
    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // final
    u.axpy(dt*b2,k2);
    t += dt;
  }

  //! return consistency order of the method
  size_type get_order () const
  {
    return 2;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::t;
  using Base::dt;
  time_type c2,a21,b2;
  using Base::u;
  V w;
  V k1,k2;
};



/** @brief Heun method (order 2 with 2 stages)

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class Heun2 : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  Heun2 (const M& model)
    : Base(model), w(model.size()), k1(model.size()), k2(model.size())
  {
    c2 = 1.0;
    a21 = 1.0;
    b1 = 0.5;
    b2 = 0.5;
    model.initialize(t,u);
    dt = 0.1;
  }

  //! do one step
  void step ()
  {
    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // final
    u.axpy(dt*b1,k1);
    u.axpy(dt*b2,k2);
    t += dt;
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 2;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type c2,a21,b1,b2;
  using Base::u;
  V w;
  V k1,k2;
};


/** @brief Heun method (order 3 with 3 stages)

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class Heun3 : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  Heun3 (const M& model)
    : Base(model), w(model.size()), k1(model.size()),
      k2(model.size()), k3(model.size())
  {
    c2 = time_type(1.0)/time_type(3.0);
    c3 = time_type(2.0)/time_type(3.0);
    a21 = time_type(1.0)/time_type(3.0);
    a32 = time_type(2.0)/time_type(3.0);
    b1 = 0.25;
    b2 = 0.0;
    b3 = 0.75;
    model.initialize(t,u);
    dt = 0.1;
  }


  //! do one step
  void step ()
  {
    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // stage 3
    w = u;
    w.axpy(dt*a32,k2);
    model.f(t+c3*dt,w,k3);

    // final
    u.axpy(dt*b1,k1);
    u.axpy(dt*b3,k3);
    t += dt;
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 3;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type c2,c3,a21,a31,a32,b1,b2,b3;
  using Base::u;
  V w;
  V k1,k2,k3;
};

/** @brief Kutta method (order 3 with 3 stages)

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class RungeKutta3 : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;

  //! constructor stores reference to the model
  RungeKutta3 (const M& model)
    : Base(model), w(model.size()), k1(model.size()),
      k2(model.size()), k3(model.size())
  {
    c2 = 0.5;
    c3 = 1.0;
    a21 = 0.5;
    a31 = -1.0;
    a32 = 2.0;
    b1 = time_type(1.0)/time_type(6.0);
    b2 = time_type(4.0)/time_type(6.0);
    b3 = time_type(1.0)/time_type(6.0);
    model.initialize(t,u);
    dt = 0.1;
  }

  //! do one step
  void step ()
  {
    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // stage 3
    w = u;
    w.axpy(dt*a31,k1);
    w.axpy(dt*a32,k2);
    model.f(t+c3*dt,w,k3);

    // final
    u.axpy(dt*b1,k1);
    u.axpy(dt*b2,k2);
    u.axpy(dt*b3,k3);
    t += dt;
  }

  //! return consistency order of the method
  size_type get_order () const
  {
    return 3;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type c2,c3,a21,a31,a32,b1,b2,b3;
  using Base::u;
  V w;
  V k1,k2,k3;
};

/** @brief classical Runge-Kutta method (order 4 with 4 stages)

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class RungeKutta4 : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  RungeKutta4 (const M& model)
    : Base(model), w(model.size()), k1(model.size()),
      k2(model.size()), k3(model.size()), k4(model.size())
  {
    c2 = 0.5;
    c3 = 0.5;
    c4 = 1.0;
    a21 = 0.5;
    a32 = 0.5;
    a43 = 1.0;
    b1 = time_type(1.0)/time_type(6.0);
    b2 = time_type(2.0)/time_type(6.0);
    b3 = time_type(2.0)/time_type(6.0);
    b4 = time_type(1.0)/time_type(6.0);
    model.initialize(t,u);
    dt = 0.1;
  }


  //! do one step
  void step ()
  {
    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // stage 3
    w = u;
    w.axpy(dt*a32,k2);
    model.f(t+c3*dt,w,k3);

    // stage 4
    w = u;
    w.axpy(dt*a43,k3);
    model.f(t+c4*dt,w,k4);

    // final
    u.axpy(dt*b1,k1);
    u.axpy(dt*b2,k2);
    u.axpy(dt*b3,k3);
    u.axpy(dt*b4,k4);
    t += dt;
  }

  //! return consistency order of the method
  size_type get_order () const
  {
    return 4;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type c2,c3,c4,a21,a32,a43,b1,b2,b3,b4;
  using Base::u;
  V w;
  V k1,k2,k3,k4;
};

/** @brief Adaptive Runge-Kutta-Fehlberg method

    \tparam M the model type
*/
template<class M>
class RKF45 : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  RKF45 (const M& model)
    : Base(model), w(model.size()), ww(model.size()), k1(model.size()),
      k2(model.size()), k3(model.size()), k4(model.size()), k5(model.size()), k6(model.size()),
      steps(0), rejected(0)
  {
    rho = time_type(0.8);
    alpha = time_type(0.25);
    beta = time_type(4.0);

    c2 = time_type(1.0)/time_type(4.0);
    c3 = time_type(3.0)/time_type(8.0);
    c4 = time_type(12.0)/time_type(13.0);
    c5 = time_type(1.0);
    c6 = time_type(1.0)/time_type(2.0);

    a21 = time_type(1.0)/time_type(4.0);

    a31 = time_type(3.0)/time_type(32.0);
    a32 = time_type(9.0)/time_type(32.0);

    a41 = time_type(1932.0)/time_type(2197.0);
    a42 = time_type(-7200.0)/time_type(2197.0);
    a43 = time_type(7296.0)/time_type(2197.0);

    a51 = time_type(439.0)/time_type(216.0);
    a52 = time_type(-8.0);
    a53 = time_type(3680.0)/time_type(513.0);
    a54 = time_type(-845.0)/time_type(4104.0);

    a61 = time_type(-8.0)/time_type(27.0);
    a62 = time_type(2.0);
    a63 = time_type(-3544.0)/time_type(2565.0);
    a64 = time_type(1859.0)/time_type(4104.0);
    a65 = time_type(-11.0)/time_type(40.0);

    b1 = time_type(25.0)/time_type(216.0);
    b2 = time_type(0.0);
    b3 = time_type(1408.0)/time_type(2565.0);
    b4 = time_type(2197.0)/time_type(4104.0);
    b5 = time_type(-1.0)/time_type(5.0);

    bb1 = time_type(16.0)/time_type(135.0);
    bb2 = time_type(0.0);
    bb3 = time_type(6656.0)/time_type(12825.0);
    bb4 = time_type(28561.0)/time_type(56430.0);
    bb5 = time_type(-9.0)/time_type(50.0);
    bb6 = time_type(2.0)/time_type(55.0);

    model.initialize(t,u);
    dt = 0.1;
  }


  //! do one step
  void step ()
  {
    steps++;

    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // stage 3
    w = u;
    w.axpy(dt*a31,k1);
    w.axpy(dt*a32,k2);
    model.f(t+c3*dt,w,k3);

    // stage 4
    w = u;
    w.axpy(dt*a41,k1);
    w.axpy(dt*a42,k2);
    w.axpy(dt*a43,k3);
    model.f(t+c4*dt,w,k4);

    // stage 5
    w = u;
    w.axpy(dt*a51,k1);
    w.axpy(dt*a52,k2);
    w.axpy(dt*a53,k3);
    w.axpy(dt*a54,k4);
    model.f(t+c5*dt,w,k5);

    // stage 6
    w = u;
    w.axpy(dt*a61,k1);
    w.axpy(dt*a62,k2);
    w.axpy(dt*a63,k3);
    w.axpy(dt*a64,k4);
    w.axpy(dt*a65,k5);
    model.f(t+c6*dt,w,k6);

    // compute order 4 approximation
    w = u;
    w.axpy(dt*b1,k1);
    w.axpy(dt*b2,k2);
    w.axpy(dt*b3,k3);
    w.axpy(dt*b4,k4);
    w.axpy(dt*b5,k5);

    // compute order 5 approximation
    ww = u;
    ww.axpy(dt*bb1,k1);
    ww.axpy(dt*bb2,k2);
    ww.axpy(dt*bb3,k3);
    ww.axpy(dt*bb4,k4);
    ww.axpy(dt*bb5,k5);
    ww.axpy(dt*bb6,k6);

    // estimate local error
    w -= ww;
    time_type error(w.two_norm()+1.E-30);
    time_type tolerr = TOL/error;
    time_type dt_opt(dt*rho*std::pow(tolerr,0.2));
    dt_opt = std::min(beta*dt,std::max(alpha*dt,dt_opt));
    // std::cout << "est. error=" << error << " dt_opt=" << dt_opt << " " << dt*rho*std::pow(tolerr,0.2) << std::endl;

    if (error<=TOL)
      {
        t += dt;
        u = ww;
        dt = dt_opt;
      }
    else
      {
        rejected++;
        dt = dt_opt;
        if (dt>dt_min) step();
        else
          {
          std::cerr << "RKF45 error " << error << " " << dt_opt << std::endl;
          for (auto it:ww)
            std::cerr << it << " ";
          DUNE_THROW(Dune::Exception, "ODEGridOperator::RKF45 error");
          }
      }
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 5;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return true;
  }

  //! print some information
  void get_info () const
  {
    std::cout << "RE: steps=" << steps << " rejected=" << rejected << std::endl;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type rho,alpha,beta;
  using Base::TOL;
  using Base::dt_min;
  time_type c2,c3,c4,c5,c6;
  time_type a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65;
  time_type b1,b2,b3,b4,b5; // 4th order
  time_type bb1,bb2,bb3,bb4,bb5,bb6; // 5th order
  using Base::u;
  V w,ww;
  V k1,k2,k3,k4,k5,k6;
  mutable size_type steps, rejected;
};



/** @brief Adaptive Runge-Kutta-Fehlberg method

    \tparam M the model type
*/
template<typename M>
class RKF23  : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;

  //! constructor stores reference to the model
  RKF23 (const M& model)
    : Base(model), w(model.size()), ww(model.size()),
      k1(model.size()), k2(model.size()), k3(model.size()),
      steps(0), rejected(0)
  {
    rho = time_type(0.8);
    alpha = time_type(0.5);
    beta = time_type(2.0);

    c2 = time_type(1.0);
    c3 = time_type(1.0)/time_type(2.0);

    a21 = time_type(1.0);

    a31 = time_type(1.0)/time_type(4.0);
    a32 = time_type(1.0)/time_type(4.0);

    b1 = time_type(1.0)/time_type(2.0);
    b2 = time_type(1.0)/time_type(2.0);

    bb1 = time_type(1.0)/time_type(6.0);
    bb2 = time_type(1.0)/time_type(6.0);
    bb3 = time_type(4.0)/time_type(6.0);

    model.initialize(t,u);
    dt = 0.1;
  }

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dt_;
  }

  //! do one step
  void step ()
  {
    steps++;

    // stage 1
    model.f(t,u,k1);

    // stage 2
    w = u;
    w.axpy(dt*a21,k1);
    model.f(t+c2*dt,w,k2);

    // stage 3
    w = u;
    w.axpy(dt*a31,k1);
    w.axpy(dt*a32,k2);
    model.f(t+c3*dt,w,k3);

    // compute order 2 approximation
    w = u;
    w.axpy(dt*b1,k1);
    w.axpy(dt*b2,k2);


    // compute order 3 approximation
    ww = u;
    ww.axpy(dt*bb1,k1);
    ww.axpy(dt*bb2,k2);
    ww.axpy(dt*bb3,k3);

    // estimate local error
    w -= ww;
    time_type error(w.two_norm()+1.E-40);
    time_type dt_opt(dt*rho*pow(TOL/error,1./3.));
    dt_opt = std::min(beta*dt,std::max(alpha*dt,dt_opt));
    //std::cout << "est. error=" << error << " dt_opt=" << dt_opt << std::endl;

    if (error<=TOL)
      {
        t += dt;
        u = ww;
        dt = dt_opt;
      }
    else
      {
        rejected++;
        dt = dt_opt;
        if (dt>dt_min) step();
        else
          std::cout << error << " " << dt_opt << std::endl;
      }

  }

  //! get current state
  const V& get_state () const
  {
    return u;
  }

  //! get current time
  time_type get_time () const
  {
    return t;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_dt () const
  {
    return dt;
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 3;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return true;
  }

  //! print some information
  void get_info () const
  {
    std::cout << "RE: steps=" << steps << " rejected=" << rejected << std::endl;
  }

private:
  using Base::model;
  using Base::dt;
  using Base::t;
  time_type rho,alpha,beta;
  using Base::TOL;
  using Base::dt_min;
  time_type c2,c3;
  time_type a21,a31,a32;
  time_type b1,b2; // 2nd order
  time_type bb1,bb2,bb3; // 3rd order
  using Base::u;
  V w,ww;
  V k1,k2,k3;
  mutable size_type steps, rejected;
};


/** @brief Adaptive one-step method using Richardson extrapolation

    \tparam M a model
    \tparam S any of the (non-adaptive) one step methods (solving model M)
*/
template<class M, class S>
class RE : public OdeBase<M>
{
public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;


  //! constructor stores reference to the model
  RE (const M& model, S& solver_)
    : Base(model), solver(std::make_shared<S>(solver_)),
      wlow(model.size()), whigh(model.size()), ww(model.size()),
      steps(0), rejected(0)
  {
    model.initialize(t,u); // initialize state
    dt = 0.1;              // set initial time step
    two_power_m = 1.0;
    for (size_type i=0; i<solver->get_order(); i++)
      two_power_m *= 2.0;
    TOL = time_type(0.0001);
    rho = time_type(0.8);
    alpha = time_type(0.25);
    beta = time_type(4.0);
    dt_min = 1E-12;
  }


  //! constructor stores reference to the model
  RE (const M& model, std::shared_ptr<S> solver_)
    : Base(model), solver(solver_),
      wlow(model.size()), whigh(model.size()), ww(model.size()),
      steps(0), rejected(0)
  {
    model.initialize(t,u); // initialize state
    dt = 0.1;              // set initial time step
    two_power_m = 1.0;
    for (size_type i=0; i<solver->get_order(); i++)
      two_power_m *= 2.0;
    TOL = time_type(0.0001);
    rho = time_type(0.8);
    alpha = time_type(0.25);
    beta = time_type(4.0);
    dt_min = 1E-12;
  }

  //! do one step
  void step ()
  {
    // count steps done
    steps++;

    // do 1 step with 2*dt
    time_type H(2.0*dt);
    solver->set_state(t,u);
    solver->set_dt(H);
    solver->step();
    wlow = solver->get_state();

    // do 2 steps with dt
    solver->set_state(t,u);
    solver->set_dt(dt);
    solver->step();
    solver->step();
    whigh = solver->get_state();

    // estimate local error
    ww = wlow;
    ww -= whigh;
    time_type error(ww.two_norm()/(std::pow(H,1.0+solver->get_order())*(1.0-1.0/two_power_m))+1.E-40);
    time_type dt_opt(std::pow(rho*TOL/error,1.0/((time_type)solver->get_order())));
    dt_opt = std::min(beta*dt,std::max(alpha*dt,dt_opt));
    //std::cout << "est. error=" << error << " dt_opt=" << dt_opt << std::endl;

    if (dt<=dt_opt)
      {
        t += H;
        u = whigh;
        u *= two_power_m;
        u -= wlow;
        u /= two_power_m-1.0;
        dt = dt_opt;
      }
    else
      {
        rejected++;
        dt = dt_opt;
        if (dt>dt_min) step();
      }
  }



  //! return consistency order of the method
  size_type get_order () const
  {
    return solver->get_order()+1;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return true;
  }

  //! print some information
  void get_info () const
  {
    std::cout << "RE: steps=" << steps << " rejected=" << rejected << std::endl;
  }

private:
  using Base::model;
  std::shared_ptr<S> solver;
  using Base::dt;
  using Base::t;
  using Base::dt_min;
  time_type two_power_m;
  using Base::u;
  V wlow,whigh,ww;
  time_type rho,alpha,beta;
  using Base::TOL;
  mutable size_type steps, rejected;
};

/** @brief Implicit Euler using Newton's method to solve nonlinear system

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
    \tparam S nonlinear solver
*/
template<class M, class S>
class ImplicitEuler : public OdeBase<M>
{

public:
  typedef OdeBase<M> Base;
  typedef typename Base::size_type size_type;
  typedef typename Base::time_type time_type;
  typedef typename Base::number_type number_type;
  typedef typename Base::V V;
  typedef typename Base::Matrix Matrix;


private:
  //! class providing nonlinear problem to be solved
  // h_n f(t_n, y_n) - y_n + y_{n-1} = 0
  class NonlinearProblem
  {
  public:

    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    typedef typename M::V V;

    typedef typename M::Matrix Matrix;

    //! constructor stores parameter lambda
    NonlinearProblem (const M& model_, const V & yold_,
                      time_type tnew_, time_type dt_)
      :model(model_), yold(yold_), tnew(tnew_), dt(dt_)
    {}

    //! return number of componentes for the model
    std::size_t size () const
    {
      return model.size();
    }

    //! model evaluation
    void F (const V& x, V& result) const
    {
      model.f(tnew,x,result);
      result *= dt;
      result -= x;
      result += yold;
    }

    //! jacobian evaluation needed for implicit solvers
    void F_x (const V& x, Matrix & result) const
    {
      model.f_x(tnew,x,result);
      result *= dt;
      for (size_type i=0; i<model.size(); i++) result[i][i] -= number_type(1.0);
    }

    void set_tnew_dt (typename M::time_type tnew_, typename M::time_type dt_)
    {
      tnew = tnew_;
      dt = dt_;
    }

  private:
    const M& model;
    const V& yold;
    typename M::time_type tnew;
    typename M::time_type dt;
  };

public:

  //! constructor stores reference to the model
  ImplicitEuler (const M& model, const S& newton_)
    : Base(model), newton(newton_), unew(model.size()), verbosity(0)
  {
    model.initialize(t,u);
    dt = dtmax = 0.1;
  }

  //! set time step for subsequent steps
  void set_dt (time_type dt_)
  {
    dt = dtmax = dt_;
  }

  //! set verbosity level
  void set_verbosity (size_type verbosity_)
  {
    verbosity = verbosity_;
  }

  //! do one step
  void step ()
  {
    if (verbosity>=2)
      std::cout << "IE: step" << " t=" << t << " dt=" << dt << std::endl;
    NonlinearProblem nlp(model,u,t+dt,dt);
    bool reduced = false;
    error = false;
    while (1)
      {
        unew = u;
        newton.solve(nlp,unew);
        if (newton.has_converged())
          {
            u = unew;
            t += dt;
            if (!reduced && dt<dtmax-dt_min/10.)
              {
                dt = std::min(2.0*dt,dtmax);
                if (verbosity>0)
                  std::cout << "IE: increasing time step to " << dt << std::endl;
              }
            return;
          }
        else
          {
            if (dt<dt_min)
              {
                DUNE_THROW(Dune::Exception, "time step too small in implicit Euler");
                error = true;
                break;
              }
            dt *= 0.5;
            reduced = true;
            nlp.set_tnew_dt(t+dt,dt);
            if (verbosity>0) std::cout << "IE: reducing time step to " << dt << std::endl;
          }
      }
  }


  //! return consistency order of the method
  size_type get_order () const
  {
    return 1;
  }

  //! return consistency order of the method
  bool adaptive () const
  {
    return false;
  }

  //! print some information
  void get_info () const
  {
  }

private:
  using Base::model;
  const S& newton;
  using Base::dt;
  using Base::t;
  using Base::dt_min;
  time_type dtmax;
  number_type reduction;
  size_type linesearchsteps;
  using Base::u;
  V unew;
  size_type verbosity;
  mutable bool error;
};



/** @brief Solve nonlinear problem using a damped Newton method

    The Newton solver is parametrized by a model. The model also
    exports all relevant types for types.

*/
class ODENewton
{

  typedef std::size_t size_type;

public:
  //! constructor stores reference to the model
  ODENewton ()
    : maxit(25), linesearchsteps(10), verbosity(3),
      reduction(1e-14), abslimit(1e-30), converged(false)
  {}

  //! maximum number of iterations before giving up
  void set_maxit (size_type n)
  {
    maxit = n;
  }

  //! maximum number of steps in linesearch before giving up
  void set_linesearchsteps (size_type n)
  {
    linesearchsteps = n;
  }

  //! control output given 0=nothing, 1=summary, 2=every step, 3=include line search
  void set_verbosity (size_type n)
  {
    verbosity = n;
  }

  //! basolute limit for defect
  void set_abslimit (double l)
  {
    abslimit = l;
  }

  //! reduction factor
  void set_reduction (double l)
  {
    reduction = l;
  }

  //! do one step
  template <class M>
  void solve (const M & model, typename M::V & x) const
  {
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    //typedef typename M::time_type time_type;

    /** \brief export number_type */
    //typedef typename M::number_type number_type;

    typedef typename M::V V;

    typedef typename M::Matrix Matrix;


    typedef typename M::number_type N;
    V r(model.size());              // residual
    Matrix A(0.); // Jacobian matrix
    V y(model.size());              // temporary solution in line search
    V z(model.size());              // solution of linear system
    V s(model.size());              // scaling factors

    model.F(x,r);                                     // compute nonlinear residual
    N R0(r.two_norm());                          // norm of initial residual
    N R(R0);                                // current residual norm
    if (verbosity>=1)
      {
        std::cout << "Newton "
                  << "   norm=" << std::scientific << std::showpoint
                  << std::setprecision(4) << R0
                  << std::endl;
      }

    converged = false;
    for (size_type i=1; i<=maxit; i++)                // do Newton iterations
      {
        // check absolute size of residual
        if (R<=abslimit)
          {
            converged = true;
            return;
          }

        // solve Jacobian system for update
        model.F_x(x,A);// compute Jacobian matrix
        std::cout << A << std::endl;
        std::cout << "r " << r << std::endl;
        std::cout << "z " << z << std::endl;
        A.solve(z,r);
        std::cout << "z " << z << std::endl;

        // line search
        N lambda(1.0);                      // start with lambda=1
        for (size_type k=0; k<linesearchsteps; k++)
          {
            y = x;
            y.axpy(-lambda,z);                       // y = x+lambda*z
            model.F(y,r);                             // r = F(y)
            N newR(r.two_norm());                // compute norm
            if (verbosity>=3)
              {
                std::cout << "    line search "  << std::setw(2) << k
                          << " lambda=" << std::scientific << std::showpoint
                          << std::setprecision(4) << lambda
                          << " norm=" << std::scientific << std::showpoint
                          << std::setprecision(4) << newR
                          << " red=" << std::scientific << std::showpoint
                          << std::setprecision(4) << newR/R
                          << std::endl;
              }
            if (newR<(1.0-0.25*lambda)*R)            // check convergence
              {
                if (verbosity>=2)
                  {
                    std::cout << "  step"  << std::setw(3) << i
                              << " norm=" << std::scientific << std::showpoint
                              << std::setprecision(4) << newR
                              << " red=" << std::scientific << std::showpoint
                              << std::setprecision(4) << newR/R
                              << std::endl;
                  }
                x = y;
                R = newR;
                break;                                // continue with Newton loop
              }
            else lambda *= 0.5;                       // reduce damping factor
            if (k==linesearchsteps-1)
              {
                if (verbosity>=3)
                  std::cout << "    line search not converged within " << linesearchsteps << " steps" << std::endl;
                return;
              }
          }

        // check convergence
        if (R<=reduction*R0)
          {
            if (verbosity>=1)
              {
                std::cout << "Newton converged in "  << i << " steps"
                          << " reduction=" << std::scientific << std::showpoint
                          << std::setprecision(4) << R/R0
                          << std::endl;
              }
            converged = true;
            return;
          }
        if (i==maxit)
          {
            if (verbosity>=2)
              std::cout << "Newton not converged within " << maxit << " iterations" << std::endl;
            converged = false;
            return;
          }
      }
  }

  bool has_converged () const
  {
    return converged;
  }

private:
  size_type maxit;
  size_type linesearchsteps;
  size_type verbosity;
  double reduction;
  double abslimit;
  mutable bool converged;
};

template<typename ODEProblem>
class OdeTimeSteppingMethods
{
public:
  typedef OdeBase<ODEProblem> ODESolver;

  OdeTimeSteppingMethods(const ODEProblem & rmodel)
  {
    methods.clear();
    // explicit methods
    methods.insert(std::make_pair("ExplicitEuler",std::shared_ptr<ODESolver>(new ExplicitEuler<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("ModifiedEuler",std::shared_ptr<ODESolver>(new ModifiedEuler<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("Heun2",std::shared_ptr<ODESolver>(new Heun2<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("Heun3",std::shared_ptr<ODESolver>(new Heun3<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("RungeKutta3",std::shared_ptr<ODESolver>(new RungeKutta3<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("RungeKutta4",std::shared_ptr<ODESolver>(new RungeKutta4<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("RKF45",std::shared_ptr<ODESolver>(new RKF45<ODEProblem>(rmodel))));
    methods.insert(std::make_pair("RKF23",std::shared_ptr<ODESolver>(new RKF23<ODEProblem>(rmodel))));

    // implicit methods
    ODENewton snewton;
    methods.insert(std::make_pair("ImplicitEuler",std::shared_ptr<ODESolver>(new ImplicitEuler<ODEProblem,ODENewton>(rmodel,snewton))));
  }

  template<typename OSM>
  void setTimestepMethod(OSM& osm, std::string method_name)
  {

    if (methods.count(method_name) )
      {
        std::shared_ptr<ODESolver> method = methods.find(method_name)->second;
        osm.setMethod(*method);
      }
    else  DUNE_THROW(Dune::Exception,"method " << method_name << "can NOT be found");
  }

private:
  std::map<std::string, std::shared_ptr<ODESolver> > methods;
};

template<typename ODEProblem>
std::shared_ptr<OdeBase<ODEProblem> > setOdeSolver(const ODEProblem & rmodel, std::string type)
{
  typedef OdeBase<ODEProblem> ODESolver;
  if (type == "ExplicitEuler") // both, but more for liquid

    //    std::cout <<
    return std::shared_ptr<ODESolver>(new ExplicitEuler<ODEProblem>(rmodel));
  else if (type == "ModifiedEuler") // liquid
    {
      return std::shared_ptr<ODESolver>(new ModifiedEuler<ODEProblem>(rmodel));
    }
  else if (type == "Heun2") //gas
    {
      return std::shared_ptr<ODESolver>(new Heun2<ODEProblem>(rmodel));
    }
  else if (type == "Heun3") //gas
    {
      return std::shared_ptr<ODESolver>(new Heun3<ODEProblem>(rmodel));
    }
  else if (type == "RungeKutta3") // both
    {
      return std::shared_ptr<ODESolver>(new RungeKutta3<ODEProblem>(rmodel));
    }
  else if (type == "RungeKutta4") // both
    {
      return std::shared_ptr<ODESolver>(new RungeKutta4<ODEProblem>(rmodel));
    }
  else if (type == "RKF45") // both
    {
      return std::shared_ptr<ODESolver>(new RKF45<ODEProblem>(rmodel));
    }
  else if (type == "RKF23") // both
    {
      return std::shared_ptr<ODESolver>(new RKF23<ODEProblem>(rmodel));
    }
  else if (type == "RE") // both
    {
      typedef RungeKutta4<ODEProblem> RK;
      std::shared_ptr<RK> rk(new RK(rmodel));

      RK rkk(rmodel);
      return std::shared_ptr<ODESolver>(new RE<ODEProblem,RK >(rmodel,rkk));
    }
  else if (type == "ImplicitEuler") // both
    {
      ODENewton snewton;
      return  std::shared_ptr<ODESolver>(new ImplicitEuler<ODEProblem,ODENewton>(rmodel,snewton));
    }
  else
    DUNE_THROW(Dune::Exception, "ODE Solver type " << type << " is not known.");
  return std::shared_ptr<ODESolver>();
}

#endif
