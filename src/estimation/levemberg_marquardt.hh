// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef LEVEMBERGMARQUARDT_HH_
#define LEVEMBERGMARQUARDT_HH_

#include<iostream>
#include "iterativefitclass.hh"

class LevMarqClass : public IterativeFitClass
{
private:
  FLOAT lambda_;

  virtual FLOAT fitStep(FLOAT **sensitivityMatrix, vector<FLOAT> &weightedResiduum);

public:
  LevMarqClass(FitDataClass &data, FitModelClass &fitmodel);
  virtual ~LevMarqClass();
  inline void SetInitialLambda(FLOAT lambda)
  {
    lambda_=lambda;
  }

};

#undef  _DEBUG_SENSITIVITY_

LevMarqClass::LevMarqClass(FitDataClass &data, FitModelClass &fitModel) :
  IterativeFitClass(data,fitModel,10.),
  lambda_(0.001)
{
}


LevMarqClass::~LevMarqClass()
{
}


FLOAT LevMarqClass::fitStep(FLOAT **sensitivityMatrix,vector<FLOAT> &weightedResiduum)
{
  vector<FLOAT> yModelTry(fitData_.GetNumValues());
  vector<FLOAT> residuumTry(fitData_.GetNumValues());
  vector<FLOAT> oldParam(fitModel_.GetParameters());
  FLOAT chiSqr;
  std::cout << "fitStep" << std::endl;

  do
    {
      cout << "lambda: " << std::scientific << lambda_ << endl;

      vector<FLOAT> correction(fitModel_.GetNumParameter());
      CalculateCorrection(sensitivityMatrix, correction, weightedResiduum, lambda_);
      CheckParameter(correction);

      for (size_t i=0;i<fitModel_.GetNumParameter();++i)
        fitModel_.SetParameter(i,oldParam[i]+correction[i]);
      cout << "new parameters: " << endl;
      fitModel_.OutputParameter(cout);
      cout << endl;
      try
        {
          chiSqr=AssembleResiduum(yModelTry,residuumTry);
        }
      catch(timestep_too_small)
        {
          cout << "time step too small, increasing lambda" << endl;
          chiSqr=chiSqr_;
          goto exit;
        }
      cout << "new residual: " << left << scientific << setw(17) << setprecision(8) << chiSqr;
      if (chiSqr >= chiSqr_)
        {
          cout << "   no improvement, new parameters rejected" << endl;
        exit:
          fitModel_.SetParameter(oldParam);
          lambda_ *= SCALE_FACTOR_;//*1.5;
        }
    } while ((chiSqr >= chiSqr_)&&(lambda_<1e20));
  FLOAT improvement=chiSqr_-chiSqr;
  if (chiSqr < chiSqr_)
    {
      lambda_ /= SCALE_FACTOR_;
      chiSqr_=chiSqr;
      yModel_=yModelTry;
      residuum_=residuumTry;
      cout << "   improvement: " << std::scientific << improvement << "   new parameters accepted" << endl << endl;
      fitData_.OutputSolution(yModel_);

      AssembleSensitivityMatrix(sensitivityMatrix,weightedResiduum,yModel_,residuum_);
    }
  cout << endl;
  return(improvement);
}

#endif
