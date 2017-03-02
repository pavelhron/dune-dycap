// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_RESIDUUM_UTILITIES_HH
#define DUNE_DYCAP_RESIDUUM_UTILITIES_HH

#include<iostream>
#include<vector>
#include<cmath>
#include<map>
#include<list>
#include<dune/common/exceptions.hh>


/**
 * \brief class to control the residuum

 * \tparam GO        grid operator
 * \tparam TrlV      test vector
 */
template<class GO, class TrlV>
class Residuum
{
  //! \brief export type for grid operator
  typedef GO GridOperator;

  //! \brief export type for test vector
  typedef TrlV TestVector;

  //! \brief export type of test vector range field
  typedef typename TestVector::ElementType RF;


public:

  //! set verbosity
  void setVerbosityLevel(unsigned int verbosity_level_)
  {
    //verobisity for each processor except 1 is zero
    if (gridoperator.trialGridFunctionSpace().gridView().comm().rank()>0)
      verbosity_level = 0;
    else
      verbosity_level = verbosity_level_;
  }

  //! default constructor
  Residuum(GridOperator& go, unsigned int verbosity_level_ = 0)
    : gridoperator(go)
    , verbosity_level(verbosity_level_)
    , r(gridoperator.testGridFunctionSpace(),0.)

  {
    if (gridoperator.trialGridFunctionSpace().gridView().comm().rank()>0)
      verbosity_level = 0;
  }

  //! get the residuum with given test vector v
  TestVector defect(TestVector& v)
  {
    r = 0.0;
    gridoperator.residual(v, r);
    RF twonorm  =  r.two_norm();

    if (!std::isfinite(twonorm))
      DUNE_THROW(Dune::Exception,
                 "Residuum::defect(): Non-linear defect is NaN or Inf");
    if (verbosity_level)
      std::cout << "residual 2 norm is " << twonorm << std::endl;
    return r;

  }

private:
  const GridOperator& gridoperator;
  unsigned int verbosity_level;
  TestVector r;
};

#endif
