/*
Copyright (c) 2003, Parallel Computing Group, Interdisciplinary Center for Scientific Comput, University of Heidelberg
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

- Neither the name of the Parallel Computing Group nor the names of its contributors
  may be used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/


//
//  File: Interpol.hh
//  Description: This file contains the definition of interpolation classe
//    used by HydrParamClass.
//  Version: 0.9
//  Author(s): Olaf Ippisch
//  mail: olaf.ippisch@iwr.uni-heidelberg.de
//

#ifndef DUNE_PM_INTERPOLATION_HH
#define DUNE_PM_INTERPOLATION_HH

#include<vector>
#include<map>
#include<iostream>
#include<fstream>
#include<float.h>

template<class DataType>
class LinearInterpolationClass
{
  private:
    std::vector<std::vector<DataType> > interpolationTable_;
    DataType interpolationInterval_;

  public:
    void ReadInterpolationTable(std::string filename)
    {
        std::ifstream infile(filename.c_str());
        int numPoints;
        infile >> numPoints;
        infile >> interpolationInterval_;
        interpolationTable_.resize(numPoints);
        for (int i=0;i<numPoints;++i)
        {
            interpolationTable_[i].resize(3);
            for (int j=0;j<3;++j)
                infile >> interpolationTable_[i][j];
        }
    }
    void ReadInterpolationValues(std::string filename)
    {
        std::ifstream infile(filename.c_str());
        int numPoints;
        infile >> numPoints;
        DataType minX, maxX;
        infile >> minX >> maxX;
        std::vector<DataType> y(numPoints);
        for (int i=0;i<numPoints;++i)
            infile >> y[i];
        Init(minX,maxX,y);
    }
    void WriteInterpolationTable(std::string filename)
    {
        std::ofstream outfile(filename.c_str());
        outfile << interpolationTable_.size() << std::endl;
        outfile << interpolationInterval_ << std::endl;
        for (size_t i=0;i<interpolationTable_.size();++i)
        {
            for (int j=0;j<3;++j)
                outfile << interpolationTable_[i][j];
            outfile << std::endl;
        }
    }
    DataType GetMinX()
    {
        if (interpolationTable_.size()>0)
            return(interpolationTable_.front()[0]);
        else
            return(0.);
    }
    DataType GetMaxX()
    {
        if (interpolationTable_.size()>0)
            return(interpolationTable_.back()[0]);
        else
            return(0.);
    }
    void Init(DataType minX, DataType maxX, std::vector<DataType>& values)
    {
        if (values.size()<2)
        {
            std::cerr << "LinearInterpolationClass::Init: empty value std::vector" << std::endl;
        }
        interpolationTable_.resize(values.size());
        if (minX>maxX)
        {
            std::swap(minX,maxX);
            interpolationInterval_=(maxX-minX)/(values.size()-1);
            for (size_t i=0;i<values.size();++i)
            {
                interpolationTable_[i].resize(3);
                interpolationTable_[i][0]=minX+i*interpolationInterval_;
                interpolationTable_[i][1]=values[values.size()-i-1];
            }
        }
        else
        {
            interpolationInterval_=(maxX-minX)/(values.size()-1);
            for (size_t i=0;i<values.size();++i)
            {
                interpolationTable_[i].resize(3);
                interpolationTable_[i][0]=minX+i*interpolationInterval_;
                interpolationTable_[i][1]=values[i];
            }
        }
        for (size_t i=0;i<values.size()-1;++i)
            interpolationTable_[i][2]=(interpolationTable_[i+1][1]-interpolationTable_[i][1])/interpolationInterval_;
    }
    DataType GetValue(DataType x)
    {
        assert(interpolationTable_.size() > 0);
        if (x<interpolationTable_.front()[0])
            return(interpolationTable_.front()[1]);
        if (x>interpolationTable_.back()[0])
            return(interpolationTable_.back()[1]);
        int index = static_cast<int>((x-interpolationTable_.front()[0])/interpolationInterval_);
        return(interpolationTable_[index][1]+
                (x-interpolationTable_[index][0])*interpolationTable_[index][2]);
    }
    DataType GetDerivative(DataType x)
    {
        assert(interpolationTable_.size() > 0);
        if (x<interpolationTable_.front()[0])
            return(0.);
        if (x>interpolationTable_.back()[0])
            return(0.);
        int index =           static_cast<int>((x-interpolationTable_.front()[0])/interpolationInterval_);
        return(interpolationTable_[index][2]);
    }
};




template<class DataType>
class CubicInterpolationClass
{
  public:
    typedef enum {BUTLAND,BRODLIE,TWOPOINT,THREEPOINT} derivativeType;
    typedef enum {HOMOTHETIC,ORTHOGONAL} projectionType;

  private:
    std::map<DataType,std::vector<DataType> > interpolationTable_;
    derivativeType derivativeSwitch_;
    projectionType projectionSwitch_;

    void CalculateDelta(const std::vector<DataType>& x, const std::vector<DataType>& y,
                        std::vector<DataType>& h, std::vector<DataType>& delta);
    void CalculateDerivative(const std::vector<DataType>& x, const std::vector<DataType>& y,
                             std::vector<DataType>& h, std::vector<DataType>& delta,
                             std::vector<DataType>& derivative);
    void CheckAndModify(std::vector<DataType>& derivative, std::vector<DataType>& delta);
    void ComputeHermiteCoefficients(const std::vector<DataType>& x, const std::vector<DataType>& y,
                                    std::vector<DataType>& h, std::vector<DataType>& delta,
                                    std::vector<DataType>& derivative);

  public:
    CubicInterpolationClass() : derivativeSwitch_(BUTLAND),
                                projectionSwitch_(HOMOTHETIC)
    {};
    void ReadInterpolationTable(std::string filename)
    {
        std::ifstream infile(filename.c_str());
        int numPoints;
        infile >> numPoints;
        std::vector<std::pair<DataType,std::vector<DataType> > > coeff(numPoints);
        for (size_t i=0;i<coeff.size();++i)
        {
            infile >> coeff[i].first;
            coeff[i].second.resize(4);
            for (int i=0;i<4;++i)
                infile>> coeff[i].second[i];
        }
        if (interpolationTable_.size()>0)
            interpolationTable_.clear();
        interpolationTable_.insert(coeff.begin(),coeff.end());
    }
    void ReadInterpolationValues(std::string filename)
    {
        std::ifstream infile(filename.c_str());
        int numPoints;
        infile >> numPoints;
        std::vector<DataType> x(numPoints);
        std::vector<DataType> y(numPoints);
        for (int i=0;i<numPoints;++i)
            infile >> x[i] >> y[i];
        Init(x,y);
    }
    void WriteInterpolationTable(std::string filename)
    {
        std::ofstream outfile(filename.c_str());
        outfile << interpolationTable_.size() << std::endl;
        for (typename std::map<DataType,std::vector<DataType> >::iterator i=interpolationTable_.begin();i!=interpolationTable_.end();++i)
        {
            outfile << i->first;
            for (size_t j=0;j<i->second.size();++j)
                outfile << " " << i->second[j];
            outfile << std::endl;
        }
    }
    void SetDerivativeType(derivativeType derivativeSwitch)
    {
        derivativeSwitch_ = derivativeSwitch;
    }
    void SetProjectionType(projectionType projectionSwitch)
    {
        projectionSwitch_ = projectionSwitch;
    }
    DataType GetMinX()
    {
        if (interpolationTable_.size()>0)
            return(interpolationTable_.begin()->first);
        else
            return(0.);
    }
    DataType GetMaxX()
    {
        if (interpolationTable_.size()>0)
            return(interpolationTable_.rbegin()->first);
        else
            return(0.);
    }
    DataType GetValue(DataType x)
    {
        typename std::map<DataType,std::vector<DataType> >::iterator elem = interpolationTable_.lower_bound(x);
        if (elem==interpolationTable_.begin())
            return(elem->second[3]);
        if (elem==interpolationTable_.end())
            return(interpolationTable_.rbegin()->second[3]);
        elem--;
        x-=elem->first;
        return ((((elem->second[0] * x) + elem->second[1]) * x + elem->second[2]) * x +
                elem->second[3]);
    }
    DataType GetDerivative(DataType x)
    {
        typename std::map<DataType,std::vector<DataType> >::iterator elem = interpolationTable_.lower_bound(x);
        if (elem==interpolationTable_.begin())
            return(0.);
        if (elem==interpolationTable_.end())
            return(0.);
        elem--;
        x-=elem->first;
        return (((3 * elem->second[0] * x) + 2 * elem->second[1]) * x + elem->second[2]);
    }
    void Init(const std::vector<DataType>& x, const std::vector<DataType>& y);
    void Init(const std::vector<DataType>& x, const std::vector<DataType>& y, const DataType firstDerivative);
    void Init(const std::vector<DataType>& x, const std::vector<DataType>& y, const DataType firstDerivative,
              DataType lastDerivative);
    void Init(const std::vector<DataType>& x, const std::vector<DataType>& y, std::vector<DataType> derivative);
    void InitLast(const std::vector<DataType>& x, const std::vector<DataType>& y, const DataType lastDerivative);
};

template<class DataType>
void CubicInterpolationClass<DataType>::CalculateDelta(const std::vector<DataType>& x,
                const std::vector<DataType>& y, std::vector<DataType>& h, std::vector<DataType>& delta)
{
    h.resize(x.size());
    delta.resize(x.size());

    for (size_t i = 0; i < delta.size()-1;++i)
    {
        h[i] = x[i+1]-x[i];
        delta[i] = (y[i+1]-y[i]) / h[i];
    }
}

template<class DataType>
void CubicInterpolationClass<DataType>::CalculateDerivative(const std::vector<DataType>& x,
                const std::vector<DataType>& y, std::vector<DataType>& h, std::vector<DataType>& delta,
                std::vector<DataType>& derivative)
{
    derivative.resize(x.size());

    if (derivativeSwitch_!=THREEPOINT)
    {
        derivative[0] = delta[0];
        derivative[x.size()-1] = delta[x.size()-2];
    }
    else
    {
        derivative[0] = (((2.*h[0] + h[1]) * delta[0]) - (h[0] * delta[1])) /
            (h[0] + h[1]);
        derivative[x.size()-1] = (((2.*h[x.size()-2] + h[x.size()-1]) * delta[x.size()-2]) -
                (h[x.size()-2] * delta[x.size()-1])) / (h[x.size()-2] + h[x.size()-1]);
    }
    for (size_t i=0;i<x.size()-2; ++i)
    {
        switch(derivativeSwitch_)
        {
            case BRODLIE:
                if ((delta[i]*delta[i+1]) > 0)
                {
                    DataType weight = (h[i] + 2. * h[i+1]) / (3. * (h[i] + h[i+1]));
                    derivative[i+1] = (delta[i] * delta[i+1]) /
                        (weight * delta[i+1] + (1. - weight) * delta[i]);
                }
                else
                    derivative[i+1] = 0.0;
                break;
            case BUTLAND:
                if ((delta[i]*delta[i+1]) > 0)
                {
                    derivative[i+1] = (2. * delta[i] * delta[i+1]) /
                        (delta[i] + delta[i+1]);
                }
                else
                    derivative[i+1] = 0.0;
                break;
            case TWOPOINT:
                derivative[i+1] = (y[i+2] - y[i]) / (x[i+2] - x[i]);
                break;
            case THREEPOINT:
                derivative[i+1] = (h[i] * delta[i] + h[i+1] * delta[i+1]) / (h[i] + h[i+1]);
                break;
            default:
                break;
        }
    }
}


template<class DataType>
void CubicInterpolationClass<DataType>::CheckAndModify(std::vector<DataType>& derivative,
                                                       std::vector<DataType>& delta)
{
    for (size_t i=0; i<delta.size()-1;++i)
    {

        if (fabs(delta[i]) < DBL_EPSILON)
        {
            derivative[i] = 0.0;
            derivative[i+1] = 0.0;
        }
        else
        {
            DataType alpha = derivative[i] / delta[i];
            DataType beta = derivative[i+1] / delta[i];

            // Values outside of DeBOOR-SCHWARTZ-BOX?
            if ((alpha < 0.) || (alpha > 3.) || (beta < 0.) || (beta > 3.))
            {
                std::cerr << "CubicSplineClass::CheckAndModify: value " << i << " out of the DeBOOR-Schwartz-Box!!\n";

                //---------------------- pair mapping-------------------------------
                if (projectionSwitch_ == HOMOTHETIC) // homothetic projection
                {
                    if ((alpha >= 0.) && (beta >= 0.))
                    {
                        // Schnittpunkt mit der y=3-Gerade
                        if (alpha<beta) // intersection with beta=3 line
                        {
                            derivative[i] = 3. * alpha/beta * delta[i];
                            derivative[i+1] = 3. * delta[i];

                        }
                        else // intersection with alpha=3 line
                        {
                            derivative[i] = 3. * delta[i];
                            derivative[i+1] = 3. * beta/alpha * delta[i];
                        }
                    }
                    else
                    {
                        if (alpha < 0.)
                        {
                            std::cerr << "CubicSplineClass::CheckAndModify: alpha negative! index " << i
                                << "  alpha: " << alpha << " = " << derivative[i]
                                << " / " << delta[i] << std::endl;
                        }
                        else
                        {
                            std::cerr << "CubicSplineClass::CheckAndModify: beta negative! index " << i
                                << "  beta: " << beta << " = " << derivative[i+1]
                                << " / " << delta[i] << std::endl;
                        }
                        throw -1; // error_abort();
                    }
                }
                else // orthogonal projection
                {
                    derivative[i] = std::min(3.,alpha) * delta[i];
                    derivative[i+1] = std::min(3.,beta) * delta[i];
                }
            }
        }
    }
}


template<class DataType>
void CubicInterpolationClass<DataType>::ComputeHermiteCoefficients(const std::vector<DataType>& x,
                const std::vector<DataType>& y, std::vector<DataType>& h, std::vector<DataType>& delta,
                std::vector<DataType>& derivative)
{
    // computation of Hermite coefficients
    std::vector<std::pair<DataType,std::vector<DataType> > > coeff(x.size());
    for (size_t i=0;i<coeff.size()-1;++i)
    {
        coeff[i].second.resize(4);
        coeff[i].first = x[i];
        coeff[i].second[2] = derivative[i];
        coeff[i].second[1] = ((-derivative[i+1] - 2.*derivative[i]) + 3.*delta[i]) / h[i];
        coeff[i].second[0] = (delta[i] - derivative[i] - coeff[i].second[1]*h[i]) / (h[i]*h[i]);
        coeff[i].second[3] = y[i];
//      coeff[i].second[3] = y[i] - (((coeff[i].second[0]*x[i]) + coeff[i].second[1])*x[i] + coeff[i].second[2]) * x[i];
    }
    coeff[coeff.size()-1].second.resize(4);
    coeff[coeff.size()-1].first = x[coeff.size()-1];
    coeff[coeff.size()-1].second[3] = y[coeff.size()-1];
    coeff[coeff.size()-1].second[2] = 0.;
    coeff[coeff.size()-1].second[1] = 0.;
    coeff[coeff.size()-1].second[0] = 0.;
    if (interpolationTable_.size()>0)
        interpolationTable_.clear();
    interpolationTable_.insert(coeff.begin(),coeff.end());
}


template<class DataType>
void CubicInterpolationClass<DataType>::Init(const std::vector<DataType>& x, const std::vector<DataType>& y)
{
    std::vector<DataType> h, delta, derivative;

    // compute h and delta
    CalculateDelta(x,y,h,delta);
    // computation of derivatives
    CalculateDerivative(x,y,h,delta,derivative);
    // check and if necessary modification of derivatives
    if ((derivativeSwitch_!=BUTLAND) && (derivativeSwitch_!=BRODLIE))
        CheckAndModify(derivative,delta);

    ComputeHermiteCoefficients(x,y,h,delta,derivative);
}


template<class DataType>
void CubicInterpolationClass<DataType>::Init(const std::vector<DataType>& x, const std::vector<DataType>& y,
                                   const DataType firstDerivative)
{
    std::vector<DataType> h, delta, derivative;

    // compute h and delta
    CalculateDelta(x,y,h,delta);
    // computation of derivatives
    CalculateDerivative(x,y,h,delta,derivative);
    derivative[0]=firstDerivative;
    // check and if necessary modification of derivatives
    if ((derivativeSwitch_!=BUTLAND) && (derivativeSwitch_!=BRODLIE))
        CheckAndModify(derivative,delta);

    ComputeHermiteCoefficients(x,y,h,delta,derivative);
}


template<class DataType>
void CubicInterpolationClass<DataType>::InitLast(const std::vector<DataType>& x, const std::vector<DataType>& y,
                                   const DataType lastDerivative)
{
    std::vector<DataType> h, delta, derivative;

    // compute h and delta
    CalculateDelta(x,y,h,delta);
    // computation of derivatives
    CalculateDerivative(x,y,h,delta,derivative);
    derivative.back()=lastDerivative;
    // check and if necessary modification of derivatives
    if ((derivativeSwitch_!=BUTLAND) && (derivativeSwitch_!=BRODLIE))
        CheckAndModify(derivative,delta);

    ComputeHermiteCoefficients(x,y,h,delta,derivative);
}


template<class DataType>
void CubicInterpolationClass<DataType>::Init(const std::vector<DataType>& x, const std::vector<DataType>& y,
                                   const DataType firstDerivative, const DataType lastDerivative)
{
    std::vector<DataType> h, delta, derivative;

    // compute h and delta
    CalculateDelta(x,y,h,delta);
    // computation of derivatives
    CalculateDerivative(x,y,h,delta,derivative);
    derivative[0]=firstDerivative;
    derivative.back()=lastDerivative;
    // check and if necessary modification of derivatives
    if ((derivativeSwitch_!=BUTLAND) && (derivativeSwitch_!=BRODLIE))
        CheckAndModify(derivative,delta);

    ComputeHermiteCoefficients(x,y,h,delta,derivative);
}


template<class DataType>
void CubicInterpolationClass<DataType>::Init(const std::vector<DataType>& x, const std::vector<DataType>& y,                                                         std::vector<DataType> derivative)
{
    std::vector<DataType> h, delta;

    // compute h and delta
    CalculateDelta(x,y,h,delta);
    // check and if necessary modification of derivatives
    CheckAndModify(derivative,delta);

    ComputeHermiteCoefficients(x,y,h,delta,derivative);
}

#endif // DUNE_PM_INTERPOLATION_HH

