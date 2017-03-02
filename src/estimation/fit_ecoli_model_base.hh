// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*- vi:
// set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARAMETERESTIMATION_ECOLIMODELBASE_HH
#define DUNE_DYCAP_PARAMETERESTIMATION_ECOLIMODELBASE_HH

#include<dune/common/exceptions.hh>
#include<fstream>
#include "fitclass.hh"
#include "utilities.hh"

//! \addtogroup ParameterEstimation
//! \ingroup Dycap
//! \{

class FitEcoliModelBase : public FitModelClass
{
protected:

  //! \struct containing types of parameters
  struct parameterType
  {
    FLOAT value;
    bool fit;
    FLOAT lowerBoundary;
    FLOAT upperBoundary;
    bool logScale;
    string name;
    int index;
  };

  //! list of all parameters
  typedef map<string,parameterType> parameterListType;

  parameterListType curParameter_;
  map<string,string> curStringParameter_;
  vector<parameterType> fitParameter_;

public:

  //! constructor
  FitEcoliModelBase(string parameterFileName, string dataFile)
  {
    std::cout << "Construct class FitEcoliModelBase" << std::endl;
    try {
      ReadParameterFile(parameterFileName);
    }
    catch (Dune::Exception e) {
      std::cout << e << std::endl;
    }
  }

  ~FitEcoliModelBase()
  {

  }

  //! writes estimated parameters to file
  void WriteParameterFile(string filename)
  {
    ofstream outfile(filename.c_str());
    if (!outfile)
      DUNE_THROW(Dune::Exception, "FitEcoliAnaerobicModel::WriteParameterFile Error: Could not open parameterfile file "+filename );
    else
      cout << "Write parameters to file " << filename << endl;


    size_t length = 0;
    for (map<string, parameterType>::iterator param = curParameter_.begin(); param
           != curParameter_.end(); ++param)
      length = max(param->second.name.size(), length);


    for (map<string, parameterType>::iterator param = curParameter_.begin(); param
           != curParameter_.end(); ++param)
      outfile << left << setw(length + 3) << param->first << setw(17) << setprecision(8) << param->second.value << setw(5) << param->second.fit << setw(17) << setprecision(8) << param->second.lowerBoundary << setw(17) << setprecision(8) << param->second.upperBoundary<< std::endl;


    //  for (map<string, parameterType>::iterator param = curParameter_.begin(); param
    //         != curParameter_.end(); ++param)
    //    {
    //      outfile << left << setw(length + 3) << param->second.name;
    //      outfile << left << setw(5) << param->second.fit << setw(17) << setprecision(8)
    //              << param->second.lowerBoundary << param->second.upperBoundary;
    //      outfile << endl;
    //    }
    //  outfile << fixed << setprecision(3);
    outfile << endl;
  }



private:

  unsigned INT GetNumParameter()
  {
    return(fitParameter_.size());
  }

  FLOAT GetParameter(INT i)
  {
    return(fitParameter_[i].value);
  }

  vector<FLOAT> GetParameters()
  {
    vector<FLOAT> result(fitParameter_.size());
    for(size_t i=0;i<fitParameter_.size();++i)
      result[i]=fitParameter_[i].value;
    return(result);
  }

  void SetParameter(INT i,FLOAT value)
  {
    fitParameter_[i].value=value;
    curParameter_[fitParameter_[i].name].value=value;
  }

  void SetParameter(vector<FLOAT> &values)
  {
    for(size_t i=0;i<fitParameter_.size();++i)
      {
        fitParameter_[i].value=values[i];
        curParameter_[fitParameter_[i].name].value=values[i];
      }
  }

  string GetParameterName(INT i)
  {
    return(fitParameter_[i].name);
  }

  bool LogScale(INT i)
  {
    return(fitParameter_[i].logScale);
  }

  FLOAT GetLowerBoundary(INT i)
  {
    return(fitParameter_[i].lowerBoundary);
  }

  void SetLowerBoundary(INT i,FLOAT value)
  {
    fitParameter_[i].lowerBoundary=value;
    curParameter_[fitParameter_[i].name].lowerBoundary=value;
  }

  FLOAT GetUpperBoundary(INT i)
  {
    return(fitParameter_[i].upperBoundary);
  }

  void SetUpperBoundary(INT i,FLOAT value)
  {
    fitParameter_[i].upperBoundary=value;
    curParameter_[fitParameter_[i].name].upperBoundary=value;
  }



  void CheckParameter(vector<FLOAT> &newParameter)
  {
  }

  void OutputParameter(ostream &ostream)
  {
    size_t length = 0;
    for (vector<parameterType>::iterator param = fitParameter_.begin(); param != fitParameter_.end(); ++param)
      length = max(param->name.size(), length);
    for (vector<parameterType>::iterator param = fitParameter_.begin(); param != fitParameter_.end(); ++param)
      ostream << left << setw(length + 3) << param->name << scientific << setw(15) << setprecision(8)
              << param->value << fixed << setprecision(3) << "(" << param->value << ")" << endl;
    ostream << endl;
  }



  void ReadParameterFile(string filename)
  {
    ifstream infile(filename.c_str());
    if (!infile)
      DUNE_THROW(Dune::Exception,  "FitEcoliAnaerobicModel::ReadParameterFile Error: Could not open parameterfile file " +filename );
    else
      cout << "ReadParameterFile reads parameters from file " << filename << endl;

    curParameter_.clear();
    curStringParameter_.clear();
    SkipComment(infile);
    do
      {
        char buffer[256];
        infile.getline(buffer, 256);
        istringstream instream(buffer);
        string parameterName;
        instream >> parameterName;

        char testChar;
        do
          {
            instream.get(testChar);
          } while (testChar == ' ');
        instream.unget();

        // Attention: Only the first character of a parameter value is tested. In order to be recognized as
        // a string parameter, the parameter value has to start with a non-numeric character!
        if (!isNumber(testChar))
          {
            std::string value;
            instream >> value;
            curStringParameter_.insert(make_pair(parameterName, value));
            std::cout << "Inserting string parameter " << parameterName << " with initial value " << value << std::endl;
          }
        else
          {
            vector <FLOAT> values;
            FLOAT value;
            while (instream >> value)
              values.push_back(value);
            if (values.size() != 1 && values.size() != 4)
              DUNE_THROW(Dune::Exception, "Illegal number of parameters in line: " + static_cast<std::string>(buffer) );

            parameterType parameter;
            parameter.name = parameterName;
            parameter.value = values[0];
            if (values.size() == 4)
              {
                if (values[1] > 0.)
                  parameter.fit = true;
                else
                  parameter.fit = false;
                parameter.lowerBoundary = values[2];
                parameter.upperBoundary = values[3];
                parameter.logScale = false;
                if (values[0]>values[3] || values[0] < values[2])
                  DUNE_THROW(Dune::Exception,  "Illegal start gues of parameter " + parameter.name );

              } else
              {
                parameter.fit = false;
                parameter.lowerBoundary = parameter.upperBoundary = 0.;
              }
            std::cout << "parameter " << parameterName << " " << scientific << setw(15) << setprecision(8) << parameter.value << " fit " << parameter.fit << std::endl;
            curParameter_.insert(make_pair(parameterName, parameter));
          }
        SkipComment(infile);
      } while (!infile.eof());

    // Copy the parameters to be fit from curParameter_ to fitParameter_
    for (map<string, parameterType>::iterator param = curParameter_.begin(); param
           != curParameter_.end(); ++param)
      if (param->second.fit)
        {
          param->second.index = fitParameter_.size();
          fitParameter_.push_back(param->second);
        }
  }

  map<string,FLOAT> getAllExperiments()
  {
    std::map<string,FLOAT> mymap;
    for (auto param = curParameter_.begin(); param != curParameter_.end(); ++param)
      mymap.insert ( std::pair<string,FLOAT>(param->second.name,param->second.value) );
    return mymap;
  }


};

//! \} group ParameterEstimation


#endif
