// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef ITERATIVECLASS_HH_
#define ITERATIVECLASS_HH_

#include<cmath>
#include<iomanip>

#include <gsl/gsl_math.h>  // Header for gsl-Lib
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define  _DEBUG_SENSITIVITY_
#include "utilities.hh"
#include "fitclass.hh"



class IterativeFitClass : public FitClass
{
protected:
  enum solverType{QRP,QRLS,SVD,GJ};
  FLOAT relLimit_;
  FLOAT absLimit_;
  INT maxIterations_;
  INT numIterations_;
  const FLOAT EPS_;
  const FLOAT INCREMENT_;
  const FLOAT SCALE_FACTOR_;
  const solverType ALGORITHM_;

  void PseudoInvers(FLOAT **covariance,INT numFitParam);
  void SolveGaussJordan(FLOAT **A, INT n, FLOAT *b);
  virtual void Initialize(FLOAT **sensitivityMatrix,vector<FLOAT> &weightedResiduum);
  virtual void AssembleSensitivityMatrix(FLOAT **matrix, vector<FLOAT> &weightedResiduum, vector<FLOAT> &yModel, vector<FLOAT> &residuum);
  virtual void CalculateCorrection(FLOAT **sensitivityMatrix, vector<FLOAT> &correction, vector<FLOAT> &weightedResiduum, FLOAT lambda);
  virtual FLOAT fitStep(FLOAT **sensitivityMatrix,vector<FLOAT> &weightedResiduum) = 0;
  FLOAT AssembleResiduum(vector<FLOAT> &yModel, vector<FLOAT> &residuum);
  void AssembleCovariance(FLOAT **sensitivityMatrix, INT numParam, FLOAT **covariance);
  virtual void CheckParameter(vector<FLOAT> &correction, FLOAT stepSize=1.);

public:
  IterativeFitClass(FitDataClass &data, FitModelClass &fitModel, FLOAT scaleFactor);
  virtual ~IterativeFitClass();
  FLOAT fit(const FLOAT precission);
  inline void SetMaxIterations(INT maxIt)
  {
    maxIterations_=maxIt;
  }
  inline void SetAbsLimit(FLOAT limit)
  {
    absLimit_=limit;
  }
  inline INT getNumberOfIterations()
  {
    return numIterations_;
  }
};

IterativeFitClass::IterativeFitClass(FitDataClass &data, FitModelClass &fitModel, FLOAT scaleFactor) :
  FitClass(data,fitModel), relLimit_(0.), absLimit_(0.001), maxIterations_(25),numIterations_(0), EPS_(1e-16), INCREMENT_(1e-4),
  SCALE_FACTOR_(scaleFactor), ALGORITHM_(SVD)
{
  for(size_t i=0;i<fitModel_.GetNumParameter();++i)
    {
      fitModel_.SetUpperBoundary(i,fitModel_.GetUpperBoundary(i)*(1.-INCREMENT_));
    }
}


IterativeFitClass::~IterativeFitClass()
{
}


void IterativeFitClass::Initialize(FLOAT **sensitivityMatrix,vector<FLOAT> &weightedResiduum)
{
  cout << "initial parameters: " << endl;
  fitModel_.OutputParameter(cout);
  cout << "end of initial parameters " << endl;
  chiSqr_=AssembleResiduum(yModel_,residuum_);
  cout << "initial residual: " << left << scientific << setw(17) << setprecision(8) << chiSqr_ << endl << endl;

  fitData_.OutputSolution(yModel_);

  relLimit_=chiSqr_;
  if (fitModel_.GetNumParameter()<1)
    {
      cerr << "IterativeFitClass error: no parameters to fit!" << endl;
      throw error_abort();
    }
  std::cout << "call sensitivity matrix \n";
  AssembleSensitivityMatrix(sensitivityMatrix,weightedResiduum,yModel_,residuum_);
  std::cout << "assemble sensitivity matrix assembled in initialize\n";
}


void IterativeFitClass::AssembleCovariance(FLOAT **sensitivityMatrix, INT numParam, FLOAT **covariance)
{
  // compute pseudo-inverse of matrix
  gsl_matrix* matrixGsl = gsl_matrix_alloc(numParam,numParam);
  gsl_matrix* vSvdGsl = gsl_matrix_alloc(numParam,numParam);
  gsl_vector* sSvdGsl = gsl_vector_alloc(numParam);
  gsl_vector* workSvdGsl = gsl_vector_alloc(numParam);
  gsl_matrix* dSvdGsl = gsl_matrix_calloc(numParam,numParam); // that's a diagonal matrix !!!
  gsl_matrix* inverseSvdGsl = gsl_matrix_alloc(numParam,numParam);

  // fill matrix with original
  for (INT j=0;j<numParam;j++)
    for (INT i=0;i<numParam;i++)
      gsl_matrix_set (matrixGsl,j,i, sensitivityMatrix[j][i]);

  // solve singular value decomposition
  gsl_linalg_SV_decomp( matrixGsl, vSvdGsl, sSvdGsl,workSvdGsl);

  gsl_matrix_transpose (matrixGsl);

  for(INT i=0;i<numParam;i++)
    gsl_matrix_set(dSvdGsl,i,i,1.0/(gsl_vector_get(sSvdGsl,i)));  // fill diagonal matrix

  // calculate D*U
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, dSvdGsl, matrixGsl, 0.0, inverseSvdGsl);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, vSvdGsl, inverseSvdGsl, 0.0, matrixGsl);

  // fill matrix with inverse
  for(INT i=0;i<numParam;i++)
    for(INT j=0;j<numParam;j++)
      covariance[i][j] = gsl_matrix_get(matrixGsl,i,j);

  gsl_matrix_free(matrixGsl);
  gsl_matrix_free(vSvdGsl);
  gsl_vector_free(sSvdGsl);
  gsl_vector_free(workSvdGsl);
  gsl_matrix_free(dSvdGsl);
  gsl_matrix_free(inverseSvdGsl);
}


void IterativeFitClass::CalculateCorrection(FLOAT **sensitivityMatrix, vector<FLOAT> &correction,
                                            vector<FLOAT> &weightedResiduum, FLOAT lambda)
{
  const INT numParam=correction.size();

  gsl_matrix* sensitivityMatrixGsl = gsl_matrix_alloc(numParam,numParam);
  gsl_vector* weightedResiduumGsl = gsl_vector_alloc(numParam);
  gsl_vector* correctionGsl = gsl_vector_alloc(numParam);
#ifdef _DEBUG_SENSITIVITY_
  std::cout << "Sensitivity Matrix: "<< std::endl;
  for (size_t i=0;i<fitModel_.GetNumParameter();i++)
    {
      for (size_t j=0;j<fitModel_.GetNumParameter();j++)
        std::cout << left << scientific << setw(17) << setprecision(8)  << sensitivityMatrix[i][j] << "  ";
      std::cout << std::endl;
    }
  std::cout << std::endl;
  std::cout << "Weighted Residuum: "<< std::endl;
  for (size_t i=0;i<fitModel_.GetNumParameter();i++)
    std::cout << left << scientific << setw(17) << setprecision(8) << weightedResiduum[i] << std::endl;
  std::cout << std::endl;
#endif

  for (INT j=0;j<numParam;j++)
    {
      for (INT k=0;k<numParam;k++)
        if (k==j)
          gsl_matrix_set (sensitivityMatrixGsl,j,k, sensitivityMatrix[j][k]*(1.+lambda));
        else
          gsl_matrix_set (sensitivityMatrixGsl,j,k, sensitivityMatrix[j][k]);
      gsl_vector_set (weightedResiduumGsl,j,weightedResiduum[j]);
    }

  switch(ALGORITHM_)
    {
    case QRP:
      { // QR mit Spaltenpivoting:
        int signumGsl;
        gsl_vector* tauGsl = gsl_vector_alloc(numParam);
        gsl_permutation* pGsl = gsl_permutation_alloc(numParam);
        gsl_vector* normGsl = gsl_vector_alloc(numParam);

        gsl_linalg_QRPT_decomp(sensitivityMatrixGsl,tauGsl,pGsl,&signumGsl,normGsl);
        gsl_linalg_QRPT_solve(sensitivityMatrixGsl,tauGsl,pGsl,weightedResiduumGsl,correctionGsl);

        gsl_permutation_free(pGsl);
        gsl_vector_free(normGsl);
        gsl_vector_free(tauGsl);
        break;
      }
    case QRLS :
      { // Least Squares QR
        gsl_vector* tauGsl = gsl_vector_alloc(numParam);
        gsl_vector* residualGsl = gsl_vector_alloc(numParam);

        gsl_linalg_QR_decomp(sensitivityMatrixGsl,tauGsl);
        gsl_linalg_QR_lssolve(sensitivityMatrixGsl,tauGsl,weightedResiduumGsl,correctionGsl,residualGsl);

        gsl_vector_free(tauGsl);
        gsl_vector_free(residualGsl);
        break;
      }
    case SVD  :
      { // Singulärwertzerlegung
        gsl_matrix* vSvdGsl = gsl_matrix_alloc(numParam,numParam);
        gsl_vector* sSvdGsl = gsl_vector_alloc(numParam);
        gsl_vector* workSvdGsl = gsl_vector_alloc(numParam);

        gsl_linalg_SV_decomp( sensitivityMatrixGsl, vSvdGsl, sSvdGsl,workSvdGsl);

        gsl_linalg_SV_solve(sensitivityMatrixGsl,vSvdGsl,sSvdGsl,weightedResiduumGsl,correctionGsl);

        gsl_matrix_free(vSvdGsl);
        gsl_vector_free(sSvdGsl);
        gsl_vector_free(workSvdGsl);
        break;
      }
    default:
      break;
    }

  for (INT j=0;j<numParam;j++)
    correction[j] =  gsl_vector_get(correctionGsl,j);
#ifdef _DEBUG_SENSITIVITY_
  std::cout << "Correction: "<< std::endl;
  for (size_t i=0;i<fitModel_.GetNumParameter();i++)
    std::cout << correction[i] << std::endl;
  std::cout << std::endl;
#endif

  gsl_matrix_free(sensitivityMatrixGsl);
  gsl_vector_free(weightedResiduumGsl);
  gsl_vector_free(correctionGsl);

  return;
}


void IterativeFitClass::AssembleSensitivityMatrix(FLOAT **sensitivityMatrix, vector<FLOAT> &weightedResiduum,
                                                  vector<FLOAT> &yModel, vector<FLOAT> &residuum)
{
  vector<FLOAT> yModelP(yModel.size());
#ifdef CENTRAL_DERIVATIVE
  vector<FLOAT> yModelM(yModel.size());
#endif
  INT size=fitModel_.GetNumParameter();
  FLOAT **dYdParam=new FLOAT *[size];
  if (dYdParam==0)
    cerr << "IterativeFitClass::AssembleSensitivityMatrix: out of memory" << endl;
  dYdParam[0]=new FLOAT[size*fitData_.GetNumValues()];
  if (dYdParam[0]==0)
    cerr << "IterativeFitClass::AssembleSensitivityMatrix: out of memory" << endl;
  for(int i=1;i<size;i++)
    dYdParam[i]=dYdParam[i-1]+fitData_.GetNumValues();


  for (INT i=0;i<size;i++)
    {
      FLOAT oldParam=fitModel_.GetParameter(i);
      FLOAT delta = INCREMENT_*fabs(oldParam)+EPS_;
      fitModel_.SetParameter(i,oldParam+delta);
      fitModel_.Model(fitData_,yModelP);
#ifdef CENTRAL_DERIVATIVE
      fitModel_.SetParameter(i,oldParam-delta);
      fitModel_.Model(fitData_,yModelM);
      for (INT j=0;j<fitData_.GetNumValues();j++)
        dYdParam[i][j]=(yModelP[j]-yModelM[j])/(2*delta);
#else
      for (INT j=0;j<fitData_.GetNumValues();j++)
        dYdParam[i][j]=(yModelP[j]-yModel[j])/delta;
#endif
      fitModel_.SetParameter(i,oldParam);
    }

  for (INT i=0;i<size;i++)
    {
      for (INT j=0;j<=i;j++)
        sensitivityMatrix[i][j]=0.0;
      weightedResiduum[i]=0.0;
    }

  for (INT j=0;j<fitData_.GetNumValues();j++)
    {
      FLOAT stdDevSqr=fitData_.GetStdDev(j);
      stdDevSqr*=stdDevSqr;
      for (INT i=0;i<size;i++)
        {
          FLOAT weight=dYdParam[i][j]/stdDevSqr;
          for (INT k=0;k<=i;k++)
            sensitivityMatrix[i][k] += weight*dYdParam[k][j];
          weightedResiduum[i] += residuum_[j]*weight;
        }
    }

  for (INT i=1;i<size;i++)
    for (INT j=0;j<i;j++)
      sensitivityMatrix[j][i]=sensitivityMatrix[i][j];

#ifdef _DEBUG_SENSITIVITY_
  std::cout << "Sensitivity Matrix: "<< std::endl;
  for (INT i=0;i<size;i++)
    {
      for (INT j=0;j<size;j++)
        std::cout << left << scientific << setw(17) << setprecision(8) << sensitivityMatrix[i][j] << "  ";
      std::cout << std::endl;
    }
  std::cout << std::endl;
#endif

  delete [] dYdParam[0];
  delete [] dYdParam;
}


FLOAT IterativeFitClass::AssembleResiduum(vector<FLOAT> &yModel, vector<FLOAT> &residuum)
{
  fitModel_.Model(fitData_,yModel);

  FLOAT result = 0.;
  FLOAT r2result = 0.0;
  FLOAT r3result = 0.0;
  FLOAT meanvalue = 0.0;
  map<string,FLOAT> partresiduum;

  for (INT i=0;i<fitData_.GetNumValues();++i)
    {
      FLOAT value = fitData_.GetValue(i);
      meanvalue+=value;
    }
  meanvalue/=static_cast<FLOAT>(fitData_.GetNumValues());


  for (INT i=0;i<fitData_.GetNumValues();++i)
    {

      FLOAT value = fitData_.GetValue(i);
      //to compute R2
      FLOAT r2value=0.0;
      FLOAT r3value=0.0;
      //do not consider negative values!!
      if (value>0.){
        residuum[i]=fitData_.GetValue(i)-yModel[i];
        r2value = fitData_.GetValue(i)-meanvalue;
        r3value = fitData_.GetValue(i);
      }
      //residuum[i]=(fitData_.GetValue(i)-yModel[i])/fitData_.GetValue(i);
      // std::cout << "fitdata[" << i << "] = " << std::setprecision(12) << fitData_.GetValue(i) << ", yModel[" << i << "] = "
      //        << yModel[i] << " --> " << residuum[i] << std::endl;
      if (fabs(fitData_.GetStdDev(i))>1e-50){
        // std::cout <<"type " <<fitData_.GetType(i)<< " fitdata[" << i << "] = " << std::setprecision(12) << fitData_.GetValue(i) << ", yModel[" << i << "] = "   << yModel[i] << " --> " << residuum[i] << std::endl;

        FLOAT res = residuum[i]*residuum[i]/(fitData_.GetStdDev(i)*fitData_.GetStdDev(i));
        FLOAT res2 = r2value*r2value/(fitData_.GetStdDev(i)*fitData_.GetStdDev(i));
        FLOAT res3 = r3value*r3value/(fitData_.GetStdDev(i)*fitData_.GetStdDev(i));
        result += res;
        r2result += res2;
        r3result += res3;
        //result +=residuum[i]*residuum[i];

        if (partresiduum.count(fitData_.GetType(i))>0)
          partresiduum.find(fitData_.GetType(i))->second +=res;
        else
          partresiduum.insert ( std::pair<string,FLOAT>(fitData_.GetType(i),res));

      }
      //result += std::abs(residuum[i]);


    }

  for (auto it=partresiduum.begin(); it!=partresiduum.end(); ++it)
    std::cout<< "partial residuum in "  << it->first << " is " << it->second << '\n';

  FLOAT r2 = 1 - result/r2result;
  FLOAT r3 = 1 - result/r3result;
  std::cout<< "coefficient of determination (R^2) is " << r2 << ", in % " << r2*100 << '\n';
  std::cout<< "modified (mean is zero) coefficient of determination (R^2) is " << r3 << ", in % " << r3*100 << '\n';
  return(result);
}


void IterativeFitClass::CheckParameter(vector<FLOAT> &correction, FLOAT stepSize)
{
  FLOAT oldStepSize=stepSize;
  vector<FLOAT> newParam(fitModel_.GetParameters());

  for (size_t i=0;i<newParam.size();++i)
    {
      if (fitModel_.LogScale(i))
        {
          if (fabs(stepSize*correction[i])>2.3)
            {
              cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " reducing correction" << endl;
              correction[i]=sign(correction[i])*0.69/stepSize;
              //                cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " keeping old value" << endl;
              //                correction[i]=0.;
              continue;
            }
        }
      else
        {
          if ((newParam[i]>DBL_EPS)&&((newParam[i]+stepSize*correction[i])/newParam[i])>10.)
            {
              cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " reducing correction" << endl;
              correction[i]=9.9*newParam[i]/stepSize;
              //                cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " keeping old value" << endl;
              //                correction[i]=0.0;
              continue;
            }
          else if ((newParam[i]>DBL_EPS)&&((newParam[i]+stepSize*correction[i])/newParam[i])<0.1)
            {
              cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " reducing correction" << endl;
              correction[i]=-0.101*newParam[i]/stepSize;
              //                cout << "Change too large in parameter " << fitModel_.GetParameterName(i) << " keeping old value" << endl;
              //                correction[i]=0.0;
              continue;
            }
        }
    }
  for (size_t i=0;i<newParam.size();++i)
    {
      if ((newParam[i]+stepSize*correction[i])<fitModel_.GetLowerBoundary(i))
        {
          cout << "Parameter " << fitModel_.GetParameterName(i) << " hit lower boundary: " << newParam[i]+stepSize*correction[i] << " <  " << fitModel_.GetLowerBoundary(i) << endl;
          newParam[i]=fitModel_.GetLowerBoundary(i);
        }
      else if ((newParam[i]+stepSize*correction[i])>fitModel_.GetUpperBoundary(i))
        {
          cout << "Parameter " << fitModel_.GetParameterName(i) << " hit upper boundary: " << newParam[i]+stepSize*correction[i]<<" > " << fitModel_.GetUpperBoundary(i) << endl;
          newParam[i]=fitModel_.GetUpperBoundary(i);
        }
      else
        newParam[i]=newParam[i]+stepSize*correction[i];
    }

  fitModel_.CheckParameter(newParam);

  for (size_t i=0;i<fitModel_.GetNumParameter();i++)
    correction[i]=(newParam[i]-fitModel_.GetParameter(i))/oldStepSize;
}


FLOAT IterativeFitClass::fit(const FLOAT precision)
{

  FLOAT **sensitivityMatrix(0);
  if (fitModel_.GetNumParameter()>0)
    {
      sensitivityMatrix= new FLOAT *[fitModel_.GetNumParameter()];
      sensitivityMatrix[0]= new FLOAT[fitModel_.GetNumParameter()*fitModel_.GetNumParameter()];
      for (size_t i=1;i<fitModel_.GetNumParameter();i++)
        sensitivityMatrix[i]= sensitivityMatrix[i-1]+fitModel_.GetNumParameter();
    }
  vector<FLOAT> weightedResiduum(fitModel_.GetNumParameter());

  std::cout << "fit method called, call initialize " << std::endl;
  Initialize(sensitivityMatrix,weightedResiduum);
  std::cout << "opt algorithm initialized" << std::endl;

  relLimit_*=precision;

  FLOAT improvement;
  numIterations_=1;
  do
    {
      cout << "step " << numIterations_++ << endl;
      improvement=fitStep(sensitivityMatrix,weightedResiduum);
    } while ((chiSqr_>relLimit_)&&(improvement>=relLimit_) && (improvement>=absLimit_) && (numIterations_<maxIterations_));
  if (chiSqr_<=relLimit_)
    cout << "Optimisation terminated, relative limit for residuum " << relLimit_ << " reached." << std::endl;
  if (improvement<relLimit_)
    cout << "Optimisation terminated, improvement " << improvement << " too small (relative limit)!" << std::endl;
  if (improvement<absLimit_)
    cout << "Optimisation terminated, improvement " << improvement << " too small (absolute limit)!" << std::endl;

  AssembleCovariance(sensitivityMatrix,fitModel_.GetNumParameter(),covariance_);

  delete [] sensitivityMatrix[0];
  delete [] sensitivityMatrix;

  if (chiSqr_>relLimit_)
    return(-chiSqr_);
  else
    return(chiSqr_);

}


#endif
