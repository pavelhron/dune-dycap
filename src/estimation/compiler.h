#ifndef _COMPILER_H_
#define _COMPILER_H_
#undef  _LONG_

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
inline long double Pow(long double a, long double b)
{
  return(pow(a,b));
}
#else
#define DOUBLE double
#define FLOAT double
inline double Pow(double a, double b)
{
  return(pow(a,b));
}
#endif

#define DBL_EPS DBL_EPSILON

#define SMALL_FAC            10
#define SMALL_D         (DBL_EPS*SMALL_FAC)

#ifdef __cplusplus
using namespace std;
#endif

inline DOUBLE HARM_MEAN(DOUBLE x,DOUBLE y)
{
  if ((x+y)==0.)
    return(0.);
  else
    return(2.*x*y/(x+y));
}


inline DOUBLE HARM_MEAN(DOUBLE x,DOUBLE weightX,DOUBLE y,DOUBLE weightY)
{
  if ((x+y)==0.)
    return(0.);
  else
    return((weightX+weightY)*x*y/(x*weightY+y*weightX));
}

inline DOUBLE GEOM_MEAN(DOUBLE x,DOUBLE y)
{
  if ((x<=0.) or (y<=0.))
    return(0.);
  else
    return(sqrt(x*y));
}

inline DOUBLE ARITHM_MEAN (DOUBLE x,DOUBLE weightX,DOUBLE y,DOUBLE weightY)
{
  if ((x+y)==0.)
    return(0.);
  else
    return(weightX*x+weightY*y);
}


class timestep_too_small
{
};

class error_abort
{
};

#endif
