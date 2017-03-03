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


#ifdef __cplusplus
#include <iostream>
#include <climits>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#else
#include <limits.h>
#include <float.h>
#endif
#include <time.h>

#define INT int
#define SHORT short

#ifdef _LONG_
#define DOUBLE long double
#define FLOAT long double
#else
#define DOUBLE double
#define FLOAT double
#endif
#define DBL_EPS DBL_EPSILON

class timestep_too_small
{
};

class error_abort
{
};


using namespace std;


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
#endif
