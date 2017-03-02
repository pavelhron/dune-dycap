//***********************************************************************
//***********************************************************************
// diffusion problem with time dependent coefficients
//***********************************************************************
//***********************************************************************

template<typename GV, typename RF>
class GenericTransportProblem
{
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;
  RF LX,LY,LZ;

public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericTransportProblem (Dune::ParameterTree& ptree) : time(0.0)
  {
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? 0 : 0;
    RF angle = ptree.get<RF>("problem.angle");
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      v[i] = 0.0;
    v[0] = cos(M_PI*angle/180.0);
    v[1] = sin(M_PI*angle/180.0);
    y0 = ptree.get<RF>("problem.y0");
    X = ptree.get<RF>("problem.X");
    Y = ptree.get<RF>("problem.Y");
    H = ptree.get<RF>("problem.H");
    LX = ptree.get<RF>("grid.structured.LX");
    LY = ptree.get<RF>("grid.structured.LY");
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return v;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    typename Traits::DomainType x = is.geometry().global(xlocal);
    if (is.outerNormal(xlocal)*v<0)
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    else
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow;
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& xlocal) const
  {
    typename Traits::DomainType x = e.geometry().global(xlocal);

    // initial condition
    if (time<1e-8)
      {
	typename Traits::DomainType center = e.geometry().center();
	if (center[0]>X && center[0]<X+H && center[1]>Y && center[1]<Y+H)
	  return 1.0;
	return 0.0;
      }

    // check for interior (happens in parallel)
    if (x[0]>1e-6 && x[0]<LX-1e-6 && x[1]>1e-6 && x[1]<LY-1e-6)
      return 0.0;

    // boundary condition for time >0
    if (x[1] < v[1]*x[0]+y0)
      return 0.0;
    else
      return 1.0;
    return 0.0;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  typename Traits::PermTensorType I;
  typename Traits::RangeType v;
  typename Traits::RangeFieldType y0;
  typename Traits::RangeFieldType X;
  typename Traits::RangeFieldType Y;
  typename Traits::RangeFieldType H;
  RF time;
};


