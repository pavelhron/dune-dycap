// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*- vi:
// set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARAMETERESTIMATION_FITCLASS_HH
#define DUNE_DYCAP_PARAMETERESTIMATION_FITCLASS_HH

#include<iostream>
#include<iomanip>
#include<vector>
#include<map>

using namespace std;

/**
   \page ParameterEstimation
   @ingroup ParameterEstimation

   \section ParameterEstimationDocIntroduction Introduction

   In the parameter estimation concept ...
*/


//! \addtogroup ParameterEstimation
//! \ingroup Dycap
//! \{

/** \brief Fit data class provides the types and virtual functions which are necessary for parameter estimation

 */
class FitDataClass
{
public:

  //! \struct for measure type
  struct measurementType
  {
    bool active; /*!< active variable */
    string type; /*!< type of variable */
    FLOAT value; /*!< variable value*/
    FLOAT stdDev; /*!< standart deviation */
    INT index; /*!< index */
    INT localindex;
  };

  //! The type of measurements types, map of vectors
  typedef map<FLOAT,vector<measurementType> > DataMapType;

  //! output solution
  virtual void OutputSolution(const vector<FLOAT> yModel)=0;

  //! get number of values we are interested in
  virtual INT GetNumValues()=0;

  //! get all the timesteps
  virtual vector<FLOAT> & GetTimes()=0;

  //! get value with index i
  virtual FLOAT GetValue(INT i)=0;

  //! get measurements for time
  virtual const vector<measurementType> GetMeasurements(FLOAT time)=0;

  //! get array with measurement
  virtual const DataMapType GetMeasurementArray(FLOAT oldTime,FLOAT newTime)=0;

  //! get standart deviation
  virtual FLOAT GetStdDev(INT i)=0;

  //! get standart deviation
  virtual string GetType(INT i)=0;

  //! get number of experiments
  virtual INT  getNumberOfExperiments()=0;

  //! default destructor
  virtual ~FitDataClass()
  {};
};

/** \brief Fit model class provides the types and virtual functions which are necessary for parameter estimation

 */
class FitModelClass
{
public:

  /*! \brief This function computes model with current \param data and returns values in \param yModel
   *
   * At first it is neccessary to set new parameters,
   * then to set initial solution for ODE problem and initial data,
   * compute the model and copy result to \param yModel
   */
  virtual void Model(FitDataClass &data, vector<FLOAT> &yModel) = 0;

  /*!
   * \brief Checks if the given parameter complies with the constraints; if not, it is modified accordingly
   *
   * \param newParameter
   */
  virtual void CheckParameter(vector<FLOAT> &newParameter){};

  /*!
   * \brief Writes the estimated parameters to the file. Is called after the computation.
   */
  virtual void OutputParameter(ostream &ostream)=0;

  //! get the number of parameters to be estimated
  virtual unsigned INT GetNumParameter()=0;

  //! get the current parameter with index i
  virtual FLOAT GetParameter(INT i)=0;

  //! get all the current parameters
  virtual vector<FLOAT> GetParameters()=0;

  //! set parameter with index parameter with value value
  virtual void SetParameter(INT parameter,FLOAT value)=0;

  //! set all the parameters with vector values
  virtual void SetParameter(vector<FLOAT> &values)=0;

  //! get name of parameter with index i
  virtual string GetParameterName(INT i)=0;

  //! get all experiments
  virtual map<string,FLOAT> getAllExperiments()=0;

  //! use logscale
  bool LogScale(INT i)
  {
    return false;
  }

  //! get lower boundary for parameter with index i
  virtual FLOAT GetLowerBoundary(INT i)=0;

  //! set lower boundary for parameter with index parameter to value
  virtual void SetLowerBoundary(INT parameter,FLOAT value)=0;

  //! get upper boundary for parameter with index i
  virtual FLOAT GetUpperBoundary(INT i)=0;

  //! set upper boundary for parameter with index parameter to value
  virtual void SetUpperBoundary(INT parameter,FLOAT value)=0;

  //! destructor
  virtual ~FitModelClass()
  {};
};


class FitClass
{
protected:
  FitDataClass &fitData_; /*!< data to be fitted */
  FitModelClass &fitModel_; /*!< model to be fitted */
  FLOAT chiSqr_; /*!< variance, first moment */
  FLOAT **covariance_; /*!<covariance, \f[ \sigma(x,y) = \operatorname{E}{\big[(x - \operatorname{E}[x])(y - \operatorname{E}[y])\big]} \f]*/
  vector<FLOAT> yModel_;
  vector<FLOAT> residuum_; /*!< residuum vector */


  //! assembles residuum from the model result
  virtual FLOAT AssembleResiduum(vector<FLOAT> &yModel, vector<FLOAT> &residuum) = 0;

  //! assembles covariance matrix from the sensitivity matrix
  virtual void AssembleCovariance(FLOAT **sensitivityMatrix, INT numParam, FLOAT **covariance) = 0;

public:

  //! constructor
  FitClass(FitDataClass &data, FitModelClass &fitModel): fitData_(data), fitModel_(fitModel), yModel_(data.GetNumValues()), residuum_(data.GetNumValues())
  {
    // if we have some parameters to estimate, allocate memory
    if (fitModel_.GetNumParameter()>0)
      {
        covariance_= new FLOAT *[fitModel_.GetNumParameter()];
        covariance_[0]= new FLOAT[fitModel_.GetNumParameter()*fitModel_.GetNumParameter()];
        for (size_t i=1;i<fitModel_.GetNumParameter();i++)
          covariance_[i]= covariance_[i-1]+fitModel_.GetNumParameter();
      }
    else
      covariance_=0;
  }

  virtual ~FitClass()
  {
    if (covariance_!=0)
      {
        delete [] covariance_[0];
        delete [] covariance_;
      }
  }

  //! fit function, important!
  virtual FLOAT fit(const FLOAT precission) = 0;

  //! some informations
  friend ostream& operator<< (ostream &ostr, FitClass &fitObj);
  friend void visualise(FitClass & fitObj, std::string fileName);
};

ostream& operator<< (ostream &ostr, FitClass &fitObj)
{
  ostr << endl << "Parameter: " << endl;

  // add parameters to ostr, name and value
  fitObj.fitModel_.OutputParameter(ostr);
  ostr << endl;

  // get the max length of parameter name
  size_t length=10;
  for(size_t i=0;i<fitObj.fitModel_.GetNumParameter();i++)
    length=max(fitObj.fitModel_.GetParameterName(i).size()+2,length);

  ostr << endl << "Standard Deviation: " << endl;
  // write sqrt of diagonal covariance (std deviation)
  for(size_t i=0;i<fitObj.fitModel_.GetNumParameter();i++)
    ostr << left << setw(length+3) << fitObj.fitModel_.GetParameterName(i) << scientific << setw(12) << setprecision(2) << sqrt(fitObj.covariance_[i][i]) << endl;
  ostr << fixed << setprecision(3);
  ostr << endl << endl;

  ostr << "Covariance Matrix: " << endl;
  // parameter names
  ostr << setw(length) << ' ';
  for(size_t i=0;i<fitObj.fitModel_.GetNumParameter()-1;i++)
    {
      if ((fitObj.fitModel_.GetParameterName(i).size()+1)<=12)
        ostr << right << setw(12) << fitObj.fitModel_.GetParameterName(i);
      else
        {
          string name;
          for (size_t j=0;j<11;++j)
            name+=fitObj.fitModel_.GetParameterName(i)[j];
          ostr << right << setw(12) << name;
        }
    }
  ostr << endl;
  // covariance matrix
  for(size_t i=1;i<fitObj.fitModel_.GetNumParameter();i++)
    {
      ostr << left << setw(length) << fitObj.fitModel_.GetParameterName(i);
      for(size_t j=0;j<i;j++)
        ostr << right << scientific << setw(12) << setprecision(2) << fitObj.covariance_[i][j];
      ostr  << endl;
    }
  ostr << endl << endl;

  ostr << "Correlation Coefficients: " << endl;

  ostr << setw(length) << ' ';
  for(size_t i=0;i<fitObj.fitModel_.GetNumParameter()-1;i++)
    {
      if ((fitObj.fitModel_.GetParameterName(i).size()+1)<=12)
        ostr << right << setw(12) << fitObj.fitModel_.GetParameterName(i);
      else
        {
          string name;
          for (size_t j=0;j<11;++j)
            name+=fitObj.fitModel_.GetParameterName(i)[j];
          ostr << right << setw(12) << name;
        }
    }
  ostr << endl;
  for(size_t i=1;i<fitObj.fitModel_.GetNumParameter();i++)
    {
      ostr << left << setw(length) << fitObj.fitModel_.GetParameterName(i);
      for(size_t j=0;j<i;j++)
        ostr << right << fixed << setw(12) << setprecision(3) << fitObj.covariance_[i][j]/sqrt(fitObj.covariance_[i][i]*fitObj.covariance_[j][j]);
      ostr  << endl;
    }
  ostr << endl;

  return ostr; // output-Stream weiterreichen (damit Verkettung moeglich)
}

//! \} group ParameterEstimation

#endif
