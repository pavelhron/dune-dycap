// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set ts=4 sw=2 et sts=2:
#ifndef DUNE_DYCAP_PARAMETERESTIMATION_UTILITIES_HH
#define DUNE_DYCAP_PARAMETERESTIMATION_UTILITIES_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <limits>
#include <vector>
#include <sstream>

using namespace std;
#include "compiler.h"



// function reads a line from a given file and skipes all lines starting with an '#' character
void SkipComment(ifstream &infile)
{
  if (infile.eof())
    return;
  do
    {
      char c=0;
      if (infile.get(c))
        {
          if (c=='#')
            {
              std::string dummy;
              getline(infile,dummy);
            }
          else if ((c!='\a')&&(c!='\b')&&(c!='\f')&&(c!='\n')&&(c!='\r')&&(c!='\t')&&(c!='\v')&&(c!=' ')&&(c!='\'')&&(c!='\"'))
            {
              infile.unget(); // decrement get pointer
              break;
            }
        }
      //  else
      //   throw std::ios_base::failure("Could not read byte in SkipComment!");
    }
  while (!infile.eof());
}

// test if is number or not
bool isNumber(char testChar)
{
  if ((testChar>='0' && testChar<='9') || testChar=='+' || testChar=='-')
    return true;
  else
    return false;
}

// determine sign of number
template<class T>
char sign(T &a)
{
  if (a>numeric_limits<T>::epsilon())
    return 1;
  else if (numeric_limits<T>::is_signed && (a < -numeric_limits<T>::epsilon()))
    return -1;
  else
    return 0;
}

void tolower(std::string &word)
{
  for (size_t i=0;i<word.size();++i)
    word[i]=tolower(word[i]);
}

/*
  template<class T>
  FLOAT GetValue(T container,string key,bool abort=true,FLOAT defValue=0.)
  {
  typename T::iterator element=container.find(key);
  if (element!=container.end())
  return(element->second);
  else
  {
  if (abort)
  {
  cerr << "GetValue: key \"" << key << "\" not found";
  cerr << endl << endl << "Available keys: " << endl;
  for(element=container.begin();element!=container.end();++element)
  cerr << element->first << endl;
  throw error_abort();
  }
  else
  {
  cerr << "GetValue: key \"" << key << "\" not found";
  cerr << " using default value " << defValue << endl;
  }
  }
  return(defValue);
  }

  string GetStringValue(std::map<string,string> container,string key,bool abort,string defValue="false")
  {
  std::map<string,string>::iterator element=container.find(key);
  if (element!=container.end())
  {
  for (size_t i=0;i<element->second.length();++i)
  element->second[i]=tolower(element->second[i]);
  return(element->second);
  }
  if (abort)
  {
  cerr << "GetValue: key \"" << key << "\" not found";
  cerr << endl << endl << "Available keys: " << endl;
  for(element=container.begin();element!=container.end();++element)
  cerr << element->first << endl;
  throw error_abort();
  }
  else
  {
  cerr << "GetValue: key \"" << key << "\" not found";
  cerr << " using default value " << defValue << endl;
  }
  return(defValue);
  }

  template<class T, class T2>
  bool GetValue(T container,std::string key, T2 &value)
  {
  typename T::iterator element=container.find(key);
  if (element!=container.end())
  {
  value=element->second;
  return(true);
  }
  else
  return(false);
  }

  template<class DataType>
  pair<DataType,bool> CheckSplineDerivative(DataType derivative, DataType x1, DataType y1, DataType x2, DataType y2)
  {
  pair<DataType,bool> result(make_pair(derivative,false));
  DataType h = x2-x1;
  DataType delta = (y2-y1) / h;

  if ((abs(delta) < DBL_EPS)&&(abs(derivative) > DBL_EPS))
  {
  result.first = 0.0;
  result.second = true;
  }
  else
  {
  DataType alpha = derivative / delta;

  if ((alpha < 0.) || (alpha > 3.))
  {
  cerr << "CheckSplineDerivative: value " << derivative << " out of the DeBOOR-Schwartz-Box!!" << endl;
  result.first = min(3.,alpha) * delta;
  cerr << "adjusting value to " << result.first << endl;
  result.second = true;
  }
  }
  return result;
  }


  void tolower(std::string &word)
  {
  for (size_t i=0;i<word.size();++i)
  word[i]=tolower(word[i]);
  }


  template<class T>
  double Blend (T x, T x0, T x1)
  {
  if (x>=x1) return 1.0;
  if (x<=x0) return 0.0;
  x = (x-x0)/(x1-x0);
  if (x<=0.5)
  return(2.*x*x);
  else
  return(1.-2.*(1.-x)*(1.-x));
  }

  typedef map<string,FLOAT> parameterListType;
  typedef map<string,string> parameterStringListType;

  void ReadParameterFile(string filename, parameterListType &curParameter, parameterStringListType &curStringParam);

  void ReadParameterFile(string filename, parameterListType &curParameter, parameterStringListType &curStringParam)
  {
  ifstream infile(filename.c_str());
  if (!infile)
  {
  cerr << "ReadParameterFile Error: Could not open parameterfile file " << filename << "!" << endl;
  throw error_abort();
  }
  SkipComment(infile);
  do
  {
  char buffer[256];
  infile.getline(buffer,256);
  istringstream instream(buffer);
  string parameterName;
  instream >> parameterName;
  tolower(parameterName);
  char testChar;
  do
  {
  instream.get(testChar);
  } while (testChar==' '||testChar=='\t');
  instream.unget();
  if (!isNumber(testChar))
  {
  std::string value;
  instream >> value;
  tolower(value);
  curStringParam.insert(make_pair(parameterName,value));
  }
  else
  {
  vector<FLOAT> values;
  FLOAT value;
  while (instream>>value)
  values.push_back(value);
  if (values.size()!=1 && values.size()!=4)
  {
  cerr << "Illegal number of parameters in line:" << endl << "\"" << buffer << "\"" << endl;
  throw error_abort();
  }
  curParameter.insert(make_pair(parameterName,values[0]));
  }
  SkipComment(infile);
  }while (!infile.eof());
  int numMaterials=GetValue(curParameter,"num_materials");
  for (int i=0;i<numMaterials;++i)
  {
  ostringstream number;
  number << "_"<< (i+1);
  string logScale;
  if (GetValue(curStringParam,"spline_conductivityl_log"+number.str(),logScale) && (logScale=="true"))
  {
  int numPoints=GetValue(curParameter,"num_spline_points_conductivityl"+number.str());
  for (int j=0;j<numPoints;++j)
  {
  ostringstream pointNumber;
  number << "_"<< (j+1);
  map<string,double>::iterator value=curParameter.find("spline_conductivityl_y"+number.str()+pointNumber.str());
  value->second=log(value->second);
  #ifdef DOUBLE_LOGARITHMIC
  value=curParameter.find("spline_conductivityl_x"+number.str()+pointNumber.str());
  value->second=log(value->second);
  #endif
  }
  numPoints=GetValue(curParameter,"num_spline_points_thetal"+number.str());
  for (int j=0;j<numPoints;++j)
  {
  ostringstream pointNumber;
  number << "_"<< (j+1);
  map<string,double>::iterator value=curParameter.find("spline_thetal_y"+number.str()+pointNumber.str());
  value->second=log(value->second);
  }
  }
  }
  }
*/
#endif
