// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef GAUSSNEWTON_HH_
#define GAUSSNEWTON_HH_

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<map>
#include"utilities.hh"
#include "fitclass.hh"
#include "iterativefitclass.hh"

class GaussNewtonClass : public IterativeFitClass
{
private:
  virtual FLOAT fitStep(FLOAT **sensitivityMatrix, vector<FLOAT> &weightedResiduum);

public:
  GaussNewtonClass(FitDataClass &data, FitModelClass &fitmodel);
  virtual ~GaussNewtonClass();
};

GaussNewtonClass::GaussNewtonClass(FitDataClass &data, FitModelClass &fitModel) : IterativeFitClass(data,fitModel,0.5)
{
}


GaussNewtonClass::~GaussNewtonClass()
{
}


FLOAT GaussNewtonClass::fitStep(FLOAT **sensitivityMatrix,vector<FLOAT> &weightedResiduum)
{
  const INT maxLineSearch=8;
  vector<vector<FLOAT> > yModelTry(maxLineSearch);
  vector<vector<FLOAT> > residuumTry(maxLineSearch);
  vector<FLOAT> rho(maxLineSearch);
  vector<FLOAT> oldParam(fitModel_.GetParameters());
  vector<FLOAT> correction(fitModel_.GetNumParameter());
  FLOAT chiSqr;

  FLOAT stepSize=1.0;
  CalculateCorrection(sensitivityMatrix, correction, weightedResiduum,0.);

  INT step=0;
  bool accept=false;

  do
    {
      vector<FLOAT> newCorrection(correction);
      CheckParameter(newCorrection,stepSize);
      for (size_t i=0;i<fitModel_.GetNumParameter();i++)
        fitModel_.SetParameter(i,oldParam[i]+stepSize*newCorrection[i]);
      try
        {
          yModelTry[step].resize(fitData_.GetNumValues());;
          residuumTry[step].resize(fitData_.GetNumValues());;
          chiSqr=AssembleResiduum(yModelTry[step],residuumTry[step]);
        }
      catch(timestep_too_small)
        {
          rho[step]=1e300;
          goto redo;
        }
      rho[step] = chiSqr/chiSqr_;

      {
        ostringstream buffer;
        buffer.precision(4);
        buffer << " ++ ls=" << setw(2) << step+1 << ", s=" << scientific << setw(12) << chiSqr;
        buffer << ", rho=" << setw(8) << rho[step];
        buffer << ", lambda=" << setw(8) << stepSize << endl;
        cout << buffer.str().c_str();
      }

      if (rho[step]<=1.-0.25*fabs(stepSize))
        {
          accept=true;
          break;
        }
      /* else reduce lambda */
    redo:   stepSize = SCALE_FACTOR_*stepSize;
      fitModel_.SetParameter(oldParam);
      step++;
    } while (step<maxLineSearch);

  if (!accept)
    {
      INT bestLineSearch = 0;
      for (step=0; step<maxLineSearch; step++)
        if (!(rho[step]==rho[step]))
          {
            bestLineSearch = step;
            break;
          }
      FLOAT rhomin = rho[bestLineSearch];
      for (step=bestLineSearch; step<maxLineSearch; step++)
        {
          if (rhomin>rho[step])
            {
              rhomin = rho[step];
              bestLineSearch = step;
            }
        }
      step=bestLineSearch;
      if (rhomin<1.)
        {
          ostringstream buffer;
          buffer << " ++ accepting linesearch " << bestLineSearch+1 << endl;
          cout << buffer.str().c_str();
        }
    }

  chiSqr=chiSqr_*rho[step];
  FLOAT improvement=chiSqr_-chiSqr;
  if (improvement>0.)
    {
      fitModel_.SetParameter(oldParam);
      stepSize = std::pow(SCALE_FACTOR_,step);
      CheckParameter(correction,stepSize);
      yModel_=yModelTry[step];
      residuum_=residuumTry[step];

      for (size_t i=0;i<fitModel_.GetNumParameter();i++)
        fitModel_.SetParameter(i,oldParam[i]+stepSize*correction[i]);

      cout << "new parameters: " << endl;
      fitModel_.OutputParameter(cout);
      cout << endl;
      chiSqr_=chiSqr;
      cout << "new residual: " << left << scientific << setw(17) << setprecision(8) << chiSqr_ << "   improvement: " << improvement << endl << endl;
      fitData_.OutputSolution(yModel_);
      AssembleSensitivityMatrix(sensitivityMatrix,weightedResiduum,yModel_,residuum_);
    }
  else
    fitModel_.SetParameter(oldParam);

  return(improvement);
}

#endif
