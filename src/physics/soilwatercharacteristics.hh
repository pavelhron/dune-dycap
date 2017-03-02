#ifndef SOILWATERCHARACTERISTICS_HH
#define SOILWATERCHARACTERISTICS_HH

#include <string>
#include <cmath>

#include <dune/common/parametertree.hh>
#include "integration_utilities.hh"

using namespace std;

/** \file

    Implementations of the Brooks-Corey, Van Genuchten and
    Simplified Van Genuchten parameterizations for the soil water characteristic
    and Mualem-Brooks-Corey, Mualem-Van Genuchten and Mualem-Simplified Van Genuchten
    parameterizations for the conductivity.
    The condictivity of the Van Genuchten parametrization (with m) is precalculated with
    a numerical Gauss-Kronrod integration and stored in a lookup table. Values between
    two calulated conductivities are interpolated linearly.

    The implementation was done according to the definitions and
    notations in the "Soil Physics - Lecture Notes" by Kurt Roth.
*/

//############################
// Capacity parameter classes
//############################

class CapacityBase
{
public:
  virtual double capacity(const double head) const {return 0;}
  virtual std::string parameterizationName() {return std::string("");}
};


class BrooksCoreyParameterization : public CapacityBase
{
private:
  double h0, lambda;

public:

  BrooksCoreyParameterization(const double h0_, const double lambda_)
    : h0(h0_), lambda(lambda_)
  {}

  BrooksCoreyParameterization(const Dune::ParameterTree & p)
    : h0(1./p.get<double>("alpha")), lambda(p.get<double>("lambda"))
  { }

  double capacity(const double head) const
  {
    if(head < h0)
      return std::pow(head/h0,-lambda);
    else
      return 1.;
  }

  std::string parameterizationName()
  {
    return std::string("BrooksCorey");
  }

  double getH0() const {return h0;}

  double getLambda() const {return lambda;}

  void setH0(double h0_) {h0 = h0_;}

  void setLambda(double lambda_) {lambda = lambda_;}
};


class SimplifiedVanGenuchtenParameterization : public CapacityBase
{
private:
  double alpha, n;
public:

  SimplifiedVanGenuchtenParameterization(const double alpha_, const double n_)
    : alpha(alpha_), n(n_)
  { }

  SimplifiedVanGenuchtenParameterization(const Dune::ParameterTree & p)
    : alpha(p.get<double>("alpha")), n(p.get<double>("n"))
  { }

  double capacity(const double head) const
  {
    if(head >= 0.)
      return 1.;

    double v = std::pow(1. + std::pow(alpha * head,n),-1.+1./n);
#ifndef NDEBUG
    if(!(v==v)){
      std::cerr << "Error in Simplified Van Genuchten: head: " << head
                << " alpha: " << alpha << " n: " << n << std::endl;
      exit(1);
    }
#endif

    return v;
  }

  std::string parameterizationName()
  {
    return std::string("SimplifiedVanGenuchten");
  }

  double getAlpha() const {return alpha;}

  double getN() const {return n;}

  void setAlpha(double alpha_) {alpha = alpha_;}

  void setN(double n_) {n = n_;}

};


class VanGenuchtenParameterization : public CapacityBase
{
private:
  double alpha, m, n;
public:

  VanGenuchtenParameterization(const double alpha_, const double m_, const double n_)
    : alpha(alpha_), m(m_), n(n_)
  { }

  VanGenuchtenParameterization(const Dune::ParameterTree & p)
    : alpha(p.get<double>("alpha")), m(p.get<double>("m")), n(p.get<double>("n"))
  { }

  double capacity(const double head) const
  {
    if(head >= 0.)
      return 1.;

    double v = std::pow(1. + std::pow(alpha * head,n),-m);
#ifndef NDEBUG
    if(!(v==v)){
      std::cerr << "Error in Van Genuchten: head: " << head
                << " alpha: " << alpha << " m: " << m << " n: " << n << std::endl;
      exit(1);
    }
#endif

    return v;
  }

  std::string parameterizationName()
  {
    return std::string("VanGenuchten");
  }

  double getAlpha() const {return alpha;}

  double getM() const {return m;}

  double getN() const {return n;}

  void setAlpha(double alpha_) {alpha = alpha_;}

  void setM(double m_) {m = m_;}

  void setN(double n_) {n = n_;}

};


//################################
// Conductivity parameter classes
//################################

class ConductivityBase
{
public:
  virtual double conductivity(const double saturation) const {return 0;}
  virtual std::string parameterizationName() {return std::string("");}
};


class MualemBrooksCoreyParameterization : public ConductivityBase
{
private:
  double a, lambda, K0, bias_factor;
public:

  MualemBrooksCoreyParameterization(const Dune::ParameterTree & p)
    : a(p.get<double>("a")), lambda(p.get<double>("lambda")), K0(p.get<double>("K")),
      bias_factor(p.get<double>("bias_factor"))
  {}


  MualemBrooksCoreyParameterization(const double a_, const double lambda_, const double K0_)
    : a(a_), lambda(lambda_), K0(K0_)
  {}

  double conductivity(const double saturation) const
  {
    if(saturation >= 1.)
      return K0 + bias_factor * K0;

    double v = K0 * std::pow(saturation, a+2.+2./lambda);

#ifndef NDEBUG
    if(!(v==v)){
      std::cerr << "Error in Mualem Brooks-Corey: saturation: " << saturation
                << " lambda: " << lambda
                << " a: " << a << " K: " << K0 << std::endl;
      exit(1);
    }
#endif
    v += bias_factor * K0;
    return v;
  }

  std::string parameterizationName()
  {
    return std::string("MualemBrooksCorey");
  }

  double getA() const {return a;}

  double getLambda() const {return lambda;}

  double getK0() const {return K0;}

  double getBiasFactor() const {return bias_factor;}

  void setA(double a_) {a = a_;}

  void setLambda(double lambda_) {lambda = lambda_;}

  void setK0(double K0_) {K0 = K0_;}

  void setBiasFactor(double bias_factor_) {bias_factor = bias_factor_;}
};


class MualemSimplifiedVanGenuchtenParameterization : public ConductivityBase
{
private:
  double alpha, a, n, K0, bias_factor;

public:

  MualemSimplifiedVanGenuchtenParameterization
  (const double alpha_, const double a_, const double n_, const double K_, const double bias_factor_)
    : alpha(alpha_), a(a_), n(n_), K0(K_), bias_factor(bias_factor_)
  {
  }

  MualemSimplifiedVanGenuchtenParameterization
  (const Dune::ParameterTree & p)
    : alpha(p.get<double>("alpha")), a(p.get<double>("a")),
      n(p.get<double>("n")), K0(p.get<double>("K")),
      bias_factor(p.get<double>("bias_factor"))
  { }

  double conductivity(const double saturation) const
  {
    if(saturation >= 1.)
      return K0 + bias_factor * K0;

    double v = K0 * std::pow(saturation,a) * std::pow(1. - std::pow(1. - std::pow(saturation,n/(n-1.)),1.-1./n),2.);
#ifndef NDEBUG
    if(!(v==v)){
      std::cerr << "Error in Mualem Simplified Van Genuchten: saturation: " << saturation
                << " alpha: " << alpha << " n: " << n
                << " n: " << n << " a: " << a << " K: " << K0 << std::endl;
      exit(1);
    }
#endif
    v += bias_factor * K0;
    return v;
  }

  std::string parameterizationName()
  {
    return std::string("MualemSimplifiedVanGenuchten");
  }

  double getAlpha() const {return alpha;}

  double getA() const {return a;}

  double getN() const {return n;}

  double getK0() const {return K0;}

  double getBiasFactor() const {return bias_factor;}

  void setAlpha(double alpha_) {alpha = alpha_;}

  void setA(double a_) {a = a_;}

  void setN(double n_) {n = n_;}

  void setK0(double K0_) {K0 = K0_;}

  void setBiasFactor(double bias_factor_) {bias_factor = bias_factor_;}
};


class MualemVanGenuchtenParameterization : public ConductivityBase
{
private:
  double alpha, a, m, n, K0, bias_factor;

  double norm;  // integration over saturation [0,1] should be done only once, this is not a mualem parameter
  double maxIntegrationError;
  int satSteps; // how much steps from saturation 0 to 1


  std::map<unsigned int,std::pair<double,double> > tableSatK;

public:

  MualemVanGenuchtenParameterization
  (const double alpha_, const double a_, const double m_, const double n_, const double K_, const double bias_factor_, const double maxIntegrationError_ = 1e-6, const double satSteps_ = 10000)
    : alpha(alpha_), a(a_), m(m_), n(n_), K0(K_), bias_factor(bias_factor_), maxIntegrationError(maxIntegrationError_), satSteps(satSteps_), oneOverHm(this)
  {
    buildLookUpTable();
  }

  MualemVanGenuchtenParameterization
  (const Dune::ParameterTree & p)
    : alpha(p.get<double>("alpha")), a(p.get<double>("a")),
      m(p.get<double>("m")), n(p.get<double>("n")), K0(p.get<double>("K")),
      bias_factor(p.get<double>("bias_factor")), maxIntegrationError(p.get<double>("max_integration_error")),
      satSteps(p.get<double>("saturation_steps")),
      oneOverHm(this)
  {
    buildLookUpTable();
  }

  double conductivity(const double saturation) const
  {
    const double & s  = saturation;

    if(saturation >= 1.)
      return K0 + bias_factor * K0;

    if(saturation <= 0.)
      return 0.0000000000000000000000000;

    unsigned int index = s*satSteps + 5;
    if (index > (satSteps-2))
      index = satSteps-2;

    while (tableSatK.at(index).first > s)
      --index;

    // linear interpolation between table values
    double m = (tableSatK.at(index + 1).second - tableSatK.at(index).second)
        / (tableSatK.at(index + 1).first - tableSatK.at(index).first);

    double c = tableSatK.at(index).second - tableSatK.at(index).first * m;

    return m * s + c;
  }

  std::string parameterizationName()
  {
    return std::string("MualemVanGenuchten");
  }

  double getAlpha() const {return alpha;}

  double getA() const {return a;}

  double getM() const {return m;}

  double getN() const {return n;}

  double getK0() const {return K0;}

  double getBiasFactor() const {return bias_factor;}

  void setAlpha(double alpha_) {alpha = alpha_;}

  void setA(double a_) {a = a_;}

  void setM(double m_) {m = m_;}

  void setN(double n_) {n = n_;}

  void setK0(double K0_) {K0 = K0_;}

  void setBiasFactor(double bias_factor_) {bias_factor = bias_factor_;}

private:

  void buildLookUpTable()
  {
    std::vector<double> saturation;
    std::vector<double> conductivity;

    double h = 1.0000000000000000/(double)satSteps;
    for (int i = 0; i < (satSteps-1); ++i)
      saturation.push_back(i*h);
    saturation.push_back(1.0000000000000000);

    norm = integrateG7K15(oneOverHm,0.0000000000000000,1.0000000000000000,maxIntegrationError);

    for (size_t i = 0; i < saturation.size(); ++i)
      conductivity.push_back(calcMualemConductivity(saturation[i]));

#ifndef NDEBUG
    // saturation and conductivity vector must have the same size
    if (saturation.size() != conductivity.size())
    {
      std::cout << "ERROR in MualemVanGenuchtenParameterization.buildLookUpTable():" << std::endl;
      std::cout << "Saturation vector and Conductivity vector don't have the same size." << std::endl;
      exit(1);
    }
#endif

    for (unsigned int i = 0; i < saturation.size(); ++i)
      tableSatK[i] = std::make_pair(saturation[i],conductivity[i]);
  }

  double calcMualemConductivity(double & sat)
  {
    if (sat <= 0.)
      return 0.;

    if (sat >= 1.0000000000000000)
      return K0 + bias_factor * K0;

    double cond = integrateG7K15(oneOverHm,0,sat,maxIntegrationError)/norm;
    cond = K0 * std::pow(sat,a) * cond * cond;

#ifndef NDEBUG
    if(!(cond==cond)){
      std::cerr << "Error in Mualem Van Genuchten: saturation: " << sat
                << " alpha: " << alpha << " m: " << m
                << " n: " << n << " a: " << a << " K: " << K0
                << " maxIntegrationError: " << maxIntegrationError << std::endl;
      exit(1);
    }
#endif

    cond += bias_factor * K0;

    return cond;
  }

  // integrand: 1/matric_head(saturation)
  class OneOverHm
  {
  public:
    OneOverHm(MualemVanGenuchtenParameterization* MvGP) : self(MvGP) {}

    double evaluate(double sat) const
    {
      return self->alpha * std::pow(std::pow(sat,-1.0000000000000000/self->m)-1.0000000000000000, -1.0000000000000000/self->n);
    }
    MualemVanGenuchtenParameterization* self;
  };

  friend class MualemVanGenuchtenParameterization::OneOverHm;

  OneOverHm oneOverHm;
};

/**
    This is a container class which provides a way to manage multiple
    parameter sets (corresponding e.g. to multiple soil materials)
    according to the Van Genuchten and Brooks-Corey soil model.
    Conductivity is calculated with Mualem's approach.

*/

class SoilWaterModel
{
public:
  typedef std::map<int,int> IndexMap;

private:
  std::vector<CapacityBase*> capacity_functions;
  std::vector<ConductivityBase*> conductivity_functions;
  std::vector<double> theta_r, theta_s;
  IndexMap material_map;
  int nsize;

public:

  SoilWaterModel(const Dune::ParameterTree & p) : nsize(0)
  {
    typedef Dune::ParameterTree::KeyVector KeyVector;
    Dune::ParameterTree material_keys = p.sub("materials");
    const KeyVector keyvector = material_keys.getSubKeys();
    KeyVector::const_iterator it = keyvector.begin();
    KeyVector::const_iterator eit = keyvector.end();
    for(; it!=eit; ++it){
      const Dune::ParameterTree sub = material_keys.Dune::ParameterTree::sub(*it);
      add(sub);
    }
  }

  void add(const Dune::ParameterTree & material)
  {
    const std::string parameterizationType = material.get<std::string>("parametrization_type");
    if(parameterizationType.find("brooks_corey") != string::npos)
    {
      capacity_functions.push_back(new BrooksCoreyParameterization(material));
      conductivity_functions.push_back(new MualemBrooksCoreyParameterization(material));
    }
    else if(parameterizationType.find("simplified_van_genuchten") != string::npos)
    {
      capacity_functions.push_back(new SimplifiedVanGenuchtenParameterization(material));
      conductivity_functions.push_back(new MualemSimplifiedVanGenuchtenParameterization(material));
    }
    else if(parameterizationType.find("van_genuchten") != string::npos)
    {
      capacity_functions.push_back(new VanGenuchtenParameterization(material));
      conductivity_functions.push_back(new MualemVanGenuchtenParameterization(material));
    }
    else
    {
      std::cout << "ERROR in RichardsSoilModel.add()" << std::endl
                << "Parameterization: \"" << parameterizationType << "\" unknown or not implemented."
                << std::endl;
      exit(1);
    }
    theta_r.push_back(material.get<double>("theta_r"));
    theta_s.push_back(material.get<double>("theta_s"));
    material_map[material.get<int>("physical_index")] = nsize;

    nsize++;
  }

  const IndexMap & getIndexMap()
  { return material_map; }

  double saturation(const int material_index, const double matric_head) const
  {
    return capacity_functions[material_index]->capacity(matric_head);
  }

  double conductivity(const int material_index, const double saturation) const
  {
    double p = conductivity_functions[material_index]->conductivity(saturation);
    return p;
  }

  double waterContent(const int material_index, const double saturation) const
  {
    return saturation * (theta_s[material_index] - theta_r[material_index]) + theta_r[material_index];
  }

  int size() const
  {
    return nsize;
  }
};

#endif
