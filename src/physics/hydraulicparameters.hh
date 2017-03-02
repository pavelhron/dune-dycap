// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=4 sw=2 et sts=4:
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


/**
   \file
   This file contains classes describing hydraulic parameters for
   for water transport.

   \author Olaf Ippisch, Christian Engwer
*/

#ifndef DUNE_PM_HYDRPAR_HH
#define DUNE_PM_HYDRPAR_HH

#include <vector>
#include <assert.h>
#include <cmath>
#include <limits>
#include <iostream>
#include "physical_chemistry.hh"
#include "interpolation.h"
#include <dune/common/exceptions.hh>

namespace Dune {
    namespace PM {

        /**
           \brief Interface class for hydraulic parameters

           \todo add generic regularization for Sl
        */
        template <typename FT>
        class HydrParamBase
        {
            typedef FT field_type;
        protected:
            field_type resSatL_;
            field_type resSatG_;
            field_type availableWater_;
            field_type entryPress_;
            field_type delta_;
            const unsigned int maxIter_;

        public:
            HydrParamBase() :
                resSatL_(0.),
                resSatG_(0.),
                availableWater_(1.),
                entryPress_(std::numeric_limits<field_type>::epsilon()),
                delta_(1e-5),
                maxIter_(200)
            {}
            virtual ~HydrParamBase()
            {}
            /**
               set entry pressure
               \param entryPress entry pressure [Pa]
            */
            void SetEntryPressure(field_type entryPress)
            {
                assert(entryPress>=0.);
                entryPress_= entryPress;
                    // * PhysicalChemistry<field_type>::DENSITY_WATER
                    // * PhysicalChemistry<field_type>::GRAVITY;
            }
            void SetResSatL(field_type resSatL)
            {
                assert((resSatL>=0.)&&(resSatL<=1.));
                resSatL_=resSatL;
                availableWater_=1.-resSatL_-resSatG_;
            }
            void SetResSatG(field_type resSatG)
            {
                assert((resSatG>=0.)&&(resSatG<=1.));
                resSatG_=resSatG;
                availableWater_=1.-resSatL_-resSatG_;
            }
            void Setdelta(field_type delta)
            {
                assert((delta>=0.)&&(delta<=1.));
                delta_=delta;
            }
            field_type GetEntryPressure()
            {
                return(entryPress_);
            }
            field_type GetResSatL()
            {
                return(resSatL_);
            }
            field_type GetResSatG()
            {
                return(resSatG_);
            }
            field_type GetAvailableWater()
            {
                return(availableWater_);
            }
            field_type Getdelta()
            {
                return(delta_);
            }
            /**
               computes saturation
            */
            virtual field_type EffSl(field_type pC)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffSl");
            }
            /**
               computes derivative of saturation
            */
            virtual field_type EffDSl(field_type pC)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffDSl");
            }
            /**
               returns the unsaturated hydraulic conductivity
            */
            virtual field_type EffKRelLpC(field_type pC)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffKRelLpC");
            }
            /**
               returns the derivative of the unsaturated hydraulic
               conductivity for the given capillary pressure
            */
            virtual field_type EffDKRelLpC(field_type pC)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffDKRelLpC");
            }
            /**
               returns the unsaturated hydraulic conductivity for the given
               volumetrical gas saturation
            */
            virtual field_type EffKRelG(field_type sG)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffKRelG");
            }
            /**
               returns the derivative of the unsaturated hydraulic
               conductivity for the given volumetrical gas saturation
            */
            virtual field_type EffDKRelG(field_type sG)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffDKRelG");
            }
            /**
               returns the unsaturated hydraulic liquid conductivity for the
               given volumetrical liquid saturation
            */
            virtual field_type EffKRelL(field_type sL)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffKRelL");
            }
            /**
               returns the derivative of unsaturated hydraulic
               liquid conductivity
               for the given volumetrical liquid saturation
            */
            virtual field_type EffDKRelL(field_type sL)
            {
                DUNE_THROW(NotImplemented, "You have to write your own EffDKRelL");
            }
            /**
               returns the matrix potential for the given
               volumetrical water content
            */
            virtual field_type EffPc(field_type sL)
            {
                if (sL>=1.)
                    return(0.);
                unsigned int iter=0;
                double xLo=entryPress_;
                double xHi=1e10;
                do
                {
                    ++iter;
                    double middle=0.5*(xLo+xHi);
                    double sLMiddle=Sl(middle);
                    if (sLMiddle<sL)
                        xHi=middle;
                    else if (std::abs(sL-sLMiddle)<std::numeric_limits<field_type>::epsilon())
                        xLo=xHi=middle;
                    else
                        xLo=middle;
                } while ((xHi-xLo)>1e-10 && iter < maxIter_);
                if ((xHi-xLo)>1e-10)
                {
                    DUNE_THROW(MathError, "root solver not converged!");
                }

                return(0.5*(xLo+xHi));
            }
            /**
               returns the matrix potential for the given
               volumetrical water content
            */
            virtual field_type EffDPc(field_type sL)
            {
                return(1./EffDSl(EffPc(sL)));
            }
            /**
               computes saturation
            */
            field_type Sl(field_type pC)
            {
                field_type Sl0 = 1.0-delta_/(1+delta_*(entryPress_-pC));
                field_type Sl1 = (1.0-delta_)*EffSl(pC);

                if (pC<=entryPress_)
                    return Sl0 * availableWater_ + resSatL_;

                return Sl1 * availableWater_ + resSatL_;
            }
            /**
               computes derivative of saturation
            */
            field_type DSl(field_type pC)
            {
                if (pC<=std::numeric_limits<field_type>::epsilon())
                    return(0.);
                if (pC<=entryPress_)
                    return -delta_*delta_/((1+delta_*(entryPress_-pC))*(1+delta_*(entryPress_-pC)));
                return(EffDSl(pC)*availableWater_);
            }
            /**
               returns the unsaturated hydraulic conductivity for the given
               volumetrical water content
            */
            field_type KRelLpC(field_type pC)
            {
                if (pC<=entryPress_)
                    return(1.);
                field_type temp = EffKRelLpC(pC);
                if (temp > std::numeric_limits<field_type>::epsilon())
                    return (temp);
                else
                    return (std::numeric_limits<field_type>::epsilon());
            }
            /**
               returns the derivative of the unsaturated hydraulic
               conductivity for the given capillary pressure
            */
            field_type DKRelLpC(field_type pC)
            {
                if (pC<=entryPress_)
                    return(0.);
                return (EffDKRelLpC(pC));
            }
            /**
               returns the unsaturated hydraulic conductivity for the given
               volumetrical gas saturation
            */
            field_type KRelG(field_type sG)
            {
                if (sG >= (1.0-resSatL_))
                    return (1.0);
                else if (sG<=resSatG_)
                    return (0.0);
                sG = (sG-resSatG_) / availableWater_;
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return(EffKRelG(sG));
            }
            /**
               returns the derivative of the unsaturated hydraulic
               conductivity for the given volumetrical gas saturation
            */
            field_type DKRelG(field_type sG)
            {
                return(EffDKRelG(sG));
            }
            /**
               returns the unsaturated hydraulic liquid conductivity for the
               given volumetrical liquid saturation
            */
            field_type KRelL(field_type sL)
            {
                if (sL >= 1.0 - resSatG_)
                    {
                        // std::cout << "sl > 1." << std::endl;
                        return (1.0);
                    }
                else if (sL <= resSatL_)
                    {
                        //std::cout << "sl < ressat" << std::endl;
                    return (0.0);
                    }
                sL = (sL-resSatL_) / availableWater_;
                field_type temp = EffKRelL(sL);
                if (temp > std::numeric_limits<field_type>::epsilon())
                    {
                        return (temp);
                    }
                else
                    {
                        //std::cout << "ret eps" << std::endl;
                    return (std::numeric_limits<field_type>::epsilon());
                    }
            }
            /**
               returns the derivative of unsaturated hydraulic liquid
               conductivity for the given volumetrical liquid saturation
            */
            field_type DKRelL(field_type sL)
            {
                if (sL >= 1.0 - resSatG_)
                    return (0.0);
                else if (sL <= resSatL_)
                    return (0.0);
                sL = (sL-resSatL_) / availableWater_;
                return(EffDKRelL(sL) / availableWater_);
            }
            /**
               returns the matrix potential for the given volumetrical water
               content
            */
            field_type Pc(field_type sL)
            {
                sL = (sL-resSatL_)/availableWater_;
                if (sL <= std::numeric_limits<field_type>::epsilon())
                    return(std::numeric_limits<field_type>::max());
                else if (sL >= 1.)
                    return(entryPress_);
                else
                    return(EffPc(sL));
            }
            /**
               returns the matrix potential for the given volumetrical water
               content
            */
            field_type DPc(field_type sL)
            {
                sL = (sL-resSatL_)/availableWater_;
                if (sL <= std::numeric_limits<field_type>::epsilon())
                    return(0.);
                else if (sL >= 1.)
                    return(0.);
                else
                    return(EffDPc(sL));
            }
        };


        template<typename FT>
        class VanGenuchtenParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            field_type alpha_;
            field_type vgN_;
            field_type vgM_;
            field_type tau_;

        public:
            VanGenuchtenParam() :
                alpha_(-1.), // [1/Pa]
                vgN_(-1.),   // [.]
                vgM_(-1.),   // [.]
                tau_(0.5)
            {}

            ~VanGenuchtenParam() {}

            void CheckValidity()
            {
                // Check if all parameters are set to a reasonable value
                if (alpha_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter alpha>0 for the van Genuchten/Mualem model. (Use functionSetAlpha)");
                if (vgN_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter n>0 for the van Genuchten/Mualem model. (Use functionSetN)");
                if (vgM_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter n>0 for the van Genuchten/Mualem model. (Use functionSetM)");
            }

            void SetEntryPressure(field_type entryPress)
            {
                DUNE_THROW(Dune::Exception,
                           "Setting EntryPressure for the van Genuchten/Mualem model is not valid");
            }

            void SetAlpha(field_type alpha)
            {
                assert(alpha>0);
                alpha_ = alpha;
            }

            void SetN(field_type vgN)
            {
                assert(vgN>1.);
                vgN_=vgN;
                if (vgN_<2.)
                    std::cerr << "BIG FAT WARNING: Using van Genuchten/Mualem parameterisation with n = "
                              << vgN_ << " < 2!" << std::endl
                              << "The van Genuchten/Mualem model is not valid in this range and can lead to numerical instability and unphysical results!!!" << std::endl;
                if (vgM_<0.)
                    vgM_=1.-1./vgN;
            }

            void SetM(field_type vgM)
            {
                assert((vgM>0.)&&(vgM<1.));
                vgM_=vgM;
            }

            void SetTau(field_type tau)
            {
                tau_=tau;
            }

            field_type EffSl(field_type pC)
            {
                ;
                return(std::pow(1.0 + std::pow(alpha_ * pC, vgN_), -vgM_));
            }

            field_type EffDSl(field_type pC)
            {
                ;
                return(-alpha_ * vgN_ * vgM_ *
                    std::pow(alpha_*pC,vgN_-1.) * std::pow(1.0+std::pow(alpha_*pC, vgN_),-vgM_-1.));
            }

            field_type EffKRelLpC(field_type pC)
            {
                ;
                pC*=alpha_;
                field_type temp = 1.+std::pow(pC, vgN_);
                return(std::pow(1.-std::pow(pC,vgN_-1.)*std::pow(temp,-vgM_),2.) * std::pow(temp,-vgM_*tau_));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                ;
                return(EffDSl(pC)*EffKRelL(EffSl(pC)));
            }

            field_type EffKRelG(field_type sG)
            {
                ;
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                // return sG;
                return (std::pow(sG,tau_) * std::pow(1.0-std::pow(1.0-sG,(1.0/vgM_)),2.0*vgM_));
            }

            field_type EffKRelL(field_type sL)
            {
                ;
                  if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                  //return pow(sL,2.);
                  return(std::pow(sL,tau_)*std::pow(1.0-std::pow(1.0-std::pow(sL,1.0/vgM_),vgM_),2.0));
            }

            field_type EffDKRelL(field_type sL)
            {
                ;
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return ((tau_*std::pow(sL, tau_-1.0) *
                        std::pow(1.0-std::pow(1.0-std::pow(sL,1.0/vgM_),vgM_),2.0) +
                        std::pow(sL,tau_)*2.0*(1.-std::pow(1.0-std::pow(sL,1.0/vgM_),vgM_)) *
                        std::pow(1.0-std::pow(sL,1.0/vgM_),vgM_-1.0)*
                        std::pow(sL,(1.0/vgM_)-1.0)) / availableWater_);
            }

            field_type EffPc(field_type sL)
            {
                ;
                return(std::pow(std::pow(sL,-1.0/vgM_)-1.0,1.0/vgN_) / alpha_);
            }
        };


        template<typename FT>
        class ModVanGenuchtenParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        protected:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            field_type alpha_;
            field_type vgN_;
            field_type vgM_;
            field_type tau_;
            field_type entrySat_;
            field_type entryCondTermL_;
            field_type entryCondTermG_;

        public:
            ModVanGenuchtenParam() :
                alpha_(-1.),
                vgN_(-1.),
                vgM_(-1.),
                tau_(0.5),
                entrySat_(0.),
                entryCondTermL_(1.),
                entryCondTermG_(1.)
            {}

            void CheckValidity()
            {
                // Check if all parameters are set to a reasonable value
                if (alpha_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter alpha>0 for the modified van Genuchten/Mualem model. (Use functionSetAlpha)");
                if (vgN_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter n>0 for the modified van Genuchten/Mualem model. (Use functionSetN)");
                if (vgM_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter n>0 for the modified van Genuchten/Mualem model. (Use functionSetM)");

                // Check if necessary functions were called
                if (entrySat_ == 0)
                    DUNE_THROW(Dune::Exception,
                               "entrySat_ of the modified van Genuchten/Mualem model is 0. Did you forget to call Init()?");
            }

            void SetAlpha(field_type alpha)
            {
                assert(alpha>0);
                alpha_= alpha;
            }

            void SetN(field_type vgN)
            {
                assert(vgN>1.);
                vgN_=vgN;
                if (vgM_<0.)
                    vgM_=1.-1./vgN;
            }

            void SetM(field_type vgM)
            {
                assert((vgM>0.)&&(vgM<1.));
                vgM_=vgM;
            }

            void SetTau(field_type tau)
            {
                tau_=tau;
            }

            void Init()
            {
                entrySat_ = std::pow(1.0 + std::pow(alpha_ * entryPress_, vgN_), -vgM_);
                entryCondTermL_ =
                    1.0 - std::pow(1.0 - std::pow(entrySat_, (1.0 / vgM_)), vgM_);
                entryCondTermG_ = std::pow((1.0 - std::pow(entrySat_, (1.0 / vgM_))), vgM_);
            }

            field_type EffSl(field_type pC)
            {
                ;
                return(std::pow(1.0 + std::pow(alpha_ * pC, vgN_), -vgM_) / entrySat_);
            }

            field_type EffDSl(field_type pC)
            {
                ;
                return(-alpha_ * vgN_ * vgM_ *
                    std::pow(alpha_*pC,vgN_-1.) * std::pow(1.0+std::pow(alpha_*pC, vgN_),-vgM_-1.) / entrySat_);
            }

            field_type EffKRelLpC(field_type pC)
            {
                ;
                pC*=alpha_;
                field_type temp = 1.+std::pow(pC, vgN_);
                return(std::pow(std::pow(temp,-vgM_)/entrySat_,tau_)) * std::pow((1.-std::pow(pC,vgN_-1.)*std::pow(temp,-vgM_)) / entryCondTermL_,2.);
            }

            field_type EffDKRelLpC(field_type pC)
            {
                ;
                return(EffDSl(pC)*EffKRelL(EffSl(pC)));
            }

            field_type EffKRelG(field_type sG)
            {
                ;
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return (std::pow(sG,tau_) *
                    std::pow((std::pow((1.0 - std::pow((1.0 - sG)*entrySat_,(1.0/vgM_))),vgM_)
                            -entryCondTermG_)/(1.-entryCondTermG_),2.));
            }

            field_type EffKRelL(field_type sL)
            {
                ;
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return(std::pow(sL,tau_)*std::pow((1.-std::pow(1.-std::pow(sL*entrySat_,1./vgM_),vgM_)) / entryCondTermL_,2.));
            }

            field_type EffDKRelL(field_type sL)
            {
                ;
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return ((tau_*std::pow(sL, tau_-1.) *
                        std::pow((1.-std::pow(1.-std::pow(sL,1./vgM_),vgM_)) / entryCondTermL_,2.) +
                        std::pow(sL,tau_)*2.*(1.-std::pow(1.-std::pow(sL,1./vgM_),vgM_)) / entryCondTermL_ *
                        std::pow(1.-std::pow(sL,1./vgM_),vgM_-1.) / entryCondTermL_ *
                        std::pow(sL,(1./vgM_)-1.)) / availableWater_);
            }

            field_type EffPc(field_type sL)
            {
                ;
                return(std::pow(std::pow(sL * entrySat_,-1./vgM_)-1.,1./vgN_) / alpha_);
            }
        };


        template<typename FT>
        class MultVanGenuchtenParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        protected:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            std::vector<field_type> alpha_;
            std::vector<field_type> vgN_;
            std::vector<field_type> vgM_;
            std::vector<field_type> fraction_;
            field_type tau_;
            field_type totPoreIntegral_;

            void SetOnlyAlpha(const std::vector<field_type>& alpha)
            {
                alpha_=alpha;
                // for (typename std::vector<field_type>::iterator i=alpha_.begin();i!=alpha_.end();++i)
                //     (*i) /= (PhysicalChemistry<field_type>::DENSITY_WATER
                //         * PhysicalChemistry<field_type>::GRAVITY);
            }

        public:
            MultVanGenuchtenParam() :
                tau_(0.5),
                totPoreIntegral_(1.)
            {}

            ~MultVanGenuchtenParam() {}

            void SetAlpha(const std::vector<field_type>& alpha)
            {
                SetOnlyAlpha(alpha);
                if (fraction_.size()==alpha_.size())
                {
                    totPoreIntegral_=0.;
                    for (size_t i=0;i<fraction_.size();++i)
                        totPoreIntegral_ += alpha_[i] * fraction_[i];
                }
            }

            void SetN(const std::vector<field_type>& vgN)
            {
                vgN_=vgN;
                for(size_t i=0;i<vgN_.size();++i)
                    if (vgN_[i]<2.)
                        std::cerr << "BIG FAT WARNING: Using van Genuchten/Mualem parameterisation with n = " << vgN_[i] << " < 2!" << std::endl << "The van Genuchten/Mualem model is not valid in this range and can lead to numerical instability and unphysical results!!!" << std::endl;
                if (vgM_.size()==0)
                {
                    vgM_=vgN;
                    for (typename std::vector<field_type>::iterator i=vgM_.begin();i!=vgM_.end();++i)
                        (*i)=1.0 - 1.0/(*i);
                }
            }

            void SetM(const std::vector<field_type>& vgM)
            {
                vgM_=vgM;
            }

            void SetFraction(const std::vector<field_type>& fraction)
            {
                fraction_=fraction;
                fraction_.push_back(0.0);
                field_type norm=0.;
                for (size_t i=0;i<fraction_.size()-1;++i)
                    norm+=fraction_[i];
                if (norm>=1.)
                {
                    fraction_.back()=0.;
                    for (size_t i=0;i<fraction_.size()-1;++i)
                        fraction_[i]/=norm;
                }
                else
                    fraction_.back()=1.-norm;
                if (fraction_.size()==alpha_.size())
                {
                    totPoreIntegral_=0.;
                    for (size_t i=0;i<fraction_.size();++i)
                        totPoreIntegral_ += alpha_[i] * fraction_[i];
                }
            }

            void SetTau(field_type tau)
            {
                tau_=tau;
            }

            field_type EffSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                    result+=fraction_[i]*std::pow(1.0 + std::pow(alpha_[i] * pC, vgN_[i]), -vgM_[i]);
                return(result);
            }

            field_type EffDSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc = alpha_[i] * pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    result+=fraction_[i]*-alpha_[i] * vgN_[i] * vgM_[i] *
                        temp/alphaPc * std::pow(1.0 + temp, -vgM_[i]-1.);
                }
                return(result);
            }

            field_type EffKRelLpC(field_type pC)
            {
                field_type poreIntegral=0.;
                field_type sL=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc=alpha_[i]*pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    field_type temp2 = std::pow(1.+temp,-vgM_[i]);
                    sL+=fraction_[i]*temp2;
                    poreIntegral+=fraction_[i]*alpha_[i]*(1. - (temp/alphaPc)*temp2);
                }
                poreIntegral/=totPoreIntegral_;
                return(std::pow(sL,tau_)*poreIntegral*poreIntegral);
            }

            field_type EffDKRelLpC(field_type pC)
            {
                field_type poreIntegral=0.;
                field_type dPoreIntegral=0.;
                field_type sL=0.;
                field_type dSl=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc=alpha_[i]*pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    field_type dTemp = alpha_[i] * vgN_[i] * temp/alphaPc;
                    field_type temp2 = std::pow(1.+temp,-vgM_[i]);
                    field_type dTemp2 = dTemp * -vgM_[i] * temp2/(1.+temp);
                    sL+=fraction_[i]*temp2;
                    dSl+=fraction_[i]*dTemp2;
                    poreIntegral+=fraction_[i]*alpha_[i]*(1. - (temp/alphaPc)*temp2);
                    dPoreIntegral-=fraction_[i]*alpha_[i]*(((dTemp-temp/pC)/alphaPc)*temp2+(temp/alphaPc)*dTemp2);
                }
                poreIntegral/=totPoreIntegral_;
                dPoreIntegral/=totPoreIntegral_;
                return(tau_*std::pow(sL,tau_-1.)*dSl*poreIntegral*poreIntegral + std::pow(sL,tau_)*2.*poreIntegral*dPoreIntegral);
            }

            field_type EffKRelG(field_type sG)
            {
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return(0.);
            }

            field_type EffKRelL(field_type sL)
            {
                return(EffKRelLpC(EffPc(sL)));
            }

            field_type EffDKRelL(field_type sL)
            {
                return(EffDPc(sL)*EffKRelLpC(EffPc(sL)));
            }
        };


        template<typename FT>
        class ModMultVanGenuchtenParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        protected:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            std::vector<field_type> alpha_;
            std::vector<field_type> vgN_;
            std::vector<field_type> vgM_;
            std::vector<field_type> fraction_;
            field_type tau_;
            field_type totPoreIntegral_;
            field_type entrySat_;
            field_type entryCondTermL_;
            field_type entryCondTermG_;

        public:
            ModMultVanGenuchtenParam() :
                tau_(0.5),
                totPoreIntegral_(1.),
                entrySat_(0.),
                entryCondTermL_(1.),
                entryCondTermG_(1.)
            {}

            ~ModMultVanGenuchtenParam() {}

            void Init()
            {
                entrySat_=0.;
                for(size_t i=0;i<alpha_.size();++i)
                    entrySat_+=fraction_[i]*std::pow(1.0 + std::pow(alpha_[i] * entryPress_, vgN_[i]), -vgM_[i]);
                entryCondTermL_=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc=alpha_[i]*entryPress_;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    field_type temp2 = std::pow(1.+temp,-vgM_[i]);
                    entryCondTermL_+=fraction_[i]*alpha_[i]*(1. - (temp/alphaPc)*temp2);
                }
                entryCondTermL_/=totPoreIntegral_;
            }

            void SetAlpha(const std::vector<field_type>& alpha)
            {
                alpha_=alpha;
                // for (typename std::vector<field_type>::iterator i=alpha_.begin();i!=alpha_.end();++i)
                //     (*i) /= (PhysicalChemistry<field_type>::DENSITY_WATER
                //         * PhysicalChemistry<field_type>::GRAVITY);
                if (fraction_.size()==alpha_.size())
                {
                    totPoreIntegral_=0.;
                    for (size_t i=0;i<fraction_.size();++i)
                        totPoreIntegral_ += alpha_[i] * fraction_[i];
                }
            }

            void SetN(const std::vector<field_type>& vgN)
            {
                vgN_=vgN;
                if (vgM_.size()==0)
                {
                    vgM_=vgN;
                    for (typename std::vector<field_type>::iterator i=vgM_.begin();i!=vgM_.end();++i)
                        (*i)=1.-1./(*i);
                }
            }

            void SetM(const std::vector<field_type>& vgM)
            {
                vgM_=vgM;
            }

            void SetFraction(const std::vector<field_type>& fraction)
            {
                fraction_=fraction;
                fraction_.push_back(0.0);
                field_type norm=0.;
                for (size_t i=0;i<fraction_.size()-1;++i)
                    norm+=fraction_[i];
                if (norm>=1.)
                {
                    fraction_.back()=0.;
                    for (size_t i=0;i<fraction_.size()-1;++i)
                        fraction_[i]/=norm;
                }
                else
                    fraction_.back()=1.-norm;
                if (fraction_.size()==alpha_.size())
                {
                    totPoreIntegral_=0.;
                    for (size_t i=0;i<fraction_.size();++i)
                        totPoreIntegral_ += alpha_[i] * fraction_[i];
                }
            }

            void SetTau(field_type tau)
            {
                tau_=tau;
            }

            field_type EffSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                    result+=fraction_[i]*std::pow(1.0 + std::pow(alpha_[i] * pC, vgN_[i]), -vgM_[i]);
                return(result/entrySat_);
            }

            field_type EffDSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc = alpha_[i] * pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    result+=fraction_[i]*-alpha_[i] * vgN_[i] * vgM_[i] *
                        temp/alphaPc * std::pow(1.0 + temp, -vgM_[i]-1.);
                }
                return(result/entrySat_);
            }

            field_type EffKRelLpC(field_type pC)
            {
                field_type poreIntegral=0.;
                field_type sL=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc=alpha_[i]*pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    field_type temp2 = std::pow(1.+temp,-vgM_[i]);
                    sL+=fraction_[i]*temp2;
                    poreIntegral+=fraction_[i]*alpha_[i]*(1. - (temp/alphaPc)*temp2);
                }
                poreIntegral/=totPoreIntegral_*entryCondTermL_;
                return(std::pow(sL/entrySat_,tau_)*poreIntegral*poreIntegral);
            }

            field_type EffDKRelLpC(field_type pC)
            {
                field_type poreIntegral=0.;
                field_type dPoreIntegral=0.;
                field_type sL=0.;
                field_type dSl=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc=alpha_[i]*pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    field_type dTemp = alpha_[i] * vgN_[i] * temp/alphaPc;
                    field_type temp2 = std::pow(1.+temp,-vgM_[i]);
                    field_type dTemp2 = dTemp * -vgM_[i] * temp2/(1.+temp);
                    sL+=fraction_[i]*temp2;
                    dSl+=fraction_[i]*dTemp2;
                    poreIntegral+=fraction_[i]*alpha_[i]*(1. - (temp/alphaPc)*temp2);
                    dPoreIntegral-=fraction_[i]*alpha_[i]*(((dTemp-temp/pC)/alphaPc)*temp2+(temp/alphaPc)*dTemp2);
                }
                poreIntegral/=totPoreIntegral_*entryCondTermL_;
                dPoreIntegral/=totPoreIntegral_*entryCondTermL_;
                return(tau_*std::pow(sL/entrySat_,tau_-1.)*dSl/entrySat_*poreIntegral*poreIntegral + std::pow(sL/entrySat_,tau_)*2.*poreIntegral*dPoreIntegral);
            }

            field_type EffKRelG(field_type sG)
            {
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return(0.);
            }

            field_type EffKRelL(field_type sL)
            {
                return(EffKRelLpC(EffPc(sL)));
            }

            field_type EffDKRelL(field_type sL)
            {
                return(EffDPc(sL)*EffKRelLpC(EffPc(sL)));
            }

        };


        template<typename FT>
        class RegModVanGenuchtenParam : public ModVanGenuchtenParam<FT>
        {
            typedef FT field_type;
        private:
            using ModVanGenuchtenParam<FT>::resSatL_;
            using ModVanGenuchtenParam<FT>::resSatG_;
            using ModVanGenuchtenParam<FT>::availableWater_;
            using ModVanGenuchtenParam<FT>::entryPress_;
            using ModVanGenuchtenParam<FT>::delta_;

            struct RegParType
            {
                field_type x0;
                field_type x1;
                std::vector<field_type> coeff;
            };

            field_type interval_;
            RegParType regPar_;

            void CheckValidity()
            {
                // Check if necessary functions were called
                if (regPar_.x1== 0)
                    DUNE_THROW(Dune::Exception,
                               "regPar_.x1 of the regularized modified van Genuchten/Mualem model is 0. Did you forget to call InitRegularisation()?");
            }

            void CompSmoothParam(std::vector<field_type> &regPar,field_type x0,field_type x1,field_type g0,field_type dg0,field_type g1,field_type dg1)
            {
                field_type t1,t2,t3,t4,t5,t6,t7,t8;
                field_type xa=(x0+x1)/2.;
                t1 = x0*x0;
                t2 = t1*xa;
                t3 = g0*xa;
                t4 = x1*x1;
                t5 = xa*x0;
                t6 = t4*xa;
                t7 = t4*x0;
                t8 = g1*xa;
                regPar[0] = (-x1*dg1-dg0*x1+xa*dg1+2.0*x0*dg0+2.0*g1-2.0*g0-xa*dg0)/(-xa+x0)/(-x1+x0)/2.0;
                regPar[1] = -(-x1*x0*dg1+xa*x0*dg1+t1*dg0+2.0*x0*g1-2.0*g0*x0-dg0*xa*x1)/(-xa+x0)/(-x1+x0);
                regPar[2] = (-t1*x1*dg1+t1*dg0*x1+t2*dg1+2.0*t1*g1+t2*dg0-2.0*x0*dg0*xa*x1+2.0*t3*x1-2.0*t3*x0-2.0*g0*x0*x1)/(-xa+x0)/(-x1+x0)/2.0;
                regPar[3] = (-xa*dg0+xa*dg1+x0*dg0+x0*dg1-2.0*g0-2.0*x1*dg1+2.0*g1)/(xa*x1-xa*x0+x0*x1-t4)/2.0;
                regPar[4] = -(-dg0*xa*x1+dg0*x0*x1-2.0*x1*g0-t4*dg1+2.0*x1*g1+t5*dg1)/(xa*x1-t5+x0*x1-t4);
                regPar[5] = (-t6*dg0-t6*dg1+t7*dg0-t7*dg1-2.0*t4*g0+2.0*x1*dg1*t5+2.0*t8*x1-2.0*t8*x0+2.0*g1*x0*x1)/(xa*x1-t5+x0*x1-t4)/2.0;
            }

            void CompSmooth()
            {
                field_type g0,g1,dg0,dg1;
                field_type x0,x1;
                field_type smoothIntervall;

                smoothIntervall=Pc(1.-resSatG_-interval_)-entryPress_;
                x0 = entryPress_;
                x1 = entryPress_+smoothIntervall;
                g0 = 1.;
                dg0 = DSl(x0);
                g1 = ModVanGenuchtenParam<FT>::EffSl(x1);
                dg1 = ModVanGenuchtenParam<FT>::EffDSl(x1);
                CompSmoothParam(regPar_.coeff,x0,x1,g0,dg0,g1,dg1);
                if ((2.*regPar_.coeff[0]*(x0+smoothIntervall*0.5)+regPar_.coeff[1])<0.)
                {
                    regPar_.x1 = entryPress_+smoothIntervall;
                    regPar_.x0 = entryPress_+0.5*smoothIntervall;
                }
                else
                {
                    std::cerr << "Deactivating smoothing pressure/saturation curve" << std::endl;
                    for (unsigned int i=0;i<6;i++)
                        std::cerr << regPar_.coeff[i] << "  ";
                    std::cerr << std::endl;
                    regPar_.x1 = 0.;
                }
            }



        public:
            RegModVanGenuchtenParam() :
                ModVanGenuchtenParam<FT>(),
                interval_(0.001)
            {
                regPar_.x1=0.;
                regPar_.coeff.resize(6);
            }

            ~RegModVanGenuchtenParam() {}

            void InitRegularisation(field_type interval=0.001)
            {
                interval_=interval;
                CompSmooth();
            }
            field_type EffSl(field_type pC)
            {
                ;
                if (pC <= regPar_.x1)
                {
                    if (pC <= regPar_.x0)
                        return(regPar_.coeff[0]*pC*pC+regPar_.coeff[1]*pC+regPar_.coeff[2]);
                    else
                        return(regPar_.coeff[3]*pC*pC+regPar_.coeff[4]*pC+regPar_.coeff[5]);
                }
                return(ModVanGenuchtenParam<FT>::EffSl(pC));
            }

            field_type EffDSl(field_type pC)
            {
                ;
                if (pC <= regPar_.x1)
                {
                    if (pC <= regPar_.x0)
                        return(2.*regPar_.coeff[0]*pC+regPar_.coeff[1]);
                    else
                        return(2.*regPar_.coeff[3]*pC+regPar_.coeff[4]);
                }
                return(ModVanGenuchtenParam<FT>::EffDSl(pC));
            }
        };


        template<typename FT>
        class MultVanGenuchtenMPermParam : public MultVanGenuchtenParam<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            using MultVanGenuchtenParam<FT>::alpha_;
            using MultVanGenuchtenParam<FT>::fraction_;
            using MultVanGenuchtenParam<FT>::tau_;
            using MultVanGenuchtenParam<FT>::totPoreIntegral_;
            field_type vgN_;
            field_type vgM_;
            field_type vgMPerm_;

        public:
            MultVanGenuchtenMPermParam() :
                vgN_(-1.),
                vgM_(-1.),
                vgMPerm_(-1.)
            {}

            ~MultVanGenuchtenMPermParam() {}

            void SetAlpha(const std::vector<field_type>& alpha)
            {
                SetOnlyAlpha(alpha);
                // alpha_=alpha;
                // for (typename std::vector<field_type>::iterator i=alpha_.begin();i!=alpha_.end();++i)
                //     (*i) /= (PhysicalChemistry<field_type>::DENSITY_WATER *
                //         PhysicalChemistry<field_type>::GRAVITY);
            }

            void SetN(field_type vgN)
            {
                vgN_=vgN;
                if (vgN_<2.)
                    std::cerr << "BIG FAT WARNING: Using van Genuchten/Mualem parameterisation with n = " << vgN_ << " < 2!" << std::endl << "The van Genuchten/Mualem model is not valid in this range and can lead to numerical instability and unphysical results!!!" << std::endl;
                if (vgM_<0.)
                    vgM_=1.-1./vgN;
                if (vgMPerm_<0.)
                    vgMPerm_=vgM_;
            }

            void SetM(field_type vgM)
            {
                vgM_=vgM;
                if (vgMPerm_<0.)
                    vgMPerm_=vgM_;
            }

            void SetMPerm(field_type vgM)
            {
                vgMPerm_=vgM;
            }

            void SetFraction(const std::vector<field_type>& fraction)
            {
                fraction_=fraction;
                fraction_.push_back(0.0);
                field_type norm=0.;
                for (size_t i=0;i<fraction_.size()-1;++i)
                    norm+=fraction_[i];
                if (norm>=1.)
                {
                    fraction_.back()=0.;
                    for (size_t i=0;i<fraction_.size()-1;++i)
                        fraction_[i]/=norm;
                }
                else
                    fraction_.back()=1.-norm;
            }

            field_type EffSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                    result+=fraction_[i]*std::pow(1.0 + std::pow(alpha_[i] * pC, vgN_), -vgM_);
                return(result);
            }

            field_type EffDSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc = alpha_[i] * pC;
                    field_type temp = std::pow(alphaPc,vgN_);
                    result+=fraction_[i]*-alpha_[i] * vgN_ * vgM_ *
                        temp/alphaPc * std::pow(1.0 + temp, -vgM_-1.);
                }
                return(result);
            }

            field_type EffKRelLpC(field_type pC)
            {
                return(EffKRelL(EffSl(pC)));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                return(EffDSl(pC)*EffDKRelL(EffSl(pC)));
            }

            field_type EffKRelL(field_type sL)
            {
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return(std::pow(sL,tau_)*std::pow(1.-std::pow(1.-std::pow(sL,1./vgMPerm_),vgMPerm_),2.));
            }

            field_type EffDKRelL(field_type sL)
            {
                if (sL < std::numeric_limits<field_type>::epsilon())
                    return(0.);
                return(tau_ * std::pow(sL,tau_-1.) * std::pow(1.-std::pow(1.-std::pow(sL,1./vgMPerm_),vgMPerm_),2.) +
                    std::pow(sL,tau_) * 2. * std::pow(1.-std::pow(1.-std::pow(sL,1./vgMPerm_),vgMPerm_),2.) *
                    vgMPerm_ * std::pow(1.-std::pow(sL,1./vgMPerm_),vgMPerm_-1.) * std::pow(sL,1./vgMPerm_-1.)/vgMPerm_);
            }
        };


        template<typename FT>
        class MultipleVanGenuchtenSplineParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        protected:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            std::vector<field_type> alpha_;
            std::vector<field_type> vgN_;
            std::vector<field_type> vgM_;
            std::vector<field_type> fraction_;

        private:
            CubicInterpolationClass<field_type> krLPcTab_;

        public:
            void ReadInterpolationValuesKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationValues(filename);
            }

            MultipleVanGenuchtenSplineParam()
            {
            }

            void SetAlpha(const std::vector<field_type>& alpha)
            {
                alpha_=alpha;
                // for (typename std::vector<field_type>::iterator i=alpha_.begin();i!=alpha_.end();++i)
                //     (*i) /= (PhysicalChemistry<field_type>::DENSITY_WATER *
                //         PhysicalChemistry<field_type>::GRAVITY);
            }

            void SetN(const std::vector<field_type>& vgN)
            {
                vgN_=vgN;
                if (vgM_.size()==0)
                {
                    vgM_=vgN;
                    for (typename std::vector<field_type>::iterator i=vgM_.begin();i!=vgM_.end();++i)
                        (*i)=1.-1./(*i);
                }
            }

            void SetM(const std::vector<field_type>& vgM)
            {
                vgM_=vgM;
            }

            void SetFraction(const std::vector<field_type>& fraction)
            {
                fraction_=fraction;
                fraction_.push_back(0.0);
                field_type norm=0.;
                for (size_t i=0;i<fraction_.size()-1;++i)
                    norm+=fraction_[i];
                if (norm>=1.)
                {
                    fraction_.back()=0.;
                    for (size_t i=0;i<fraction_.size()-1;++i)
                        fraction_[i]/=norm;
                }
                else
                    fraction_.back()=1.-norm;
            }

            field_type EffSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                    result+=fraction_[i]*std::pow(1.0 + std::pow(alpha_[i] * pC, vgN_[i]), -vgM_[i]);
                return(result);
            }

            field_type EffDSl(field_type pC)
            {
                field_type result=0.;
                for(size_t i=0;i<alpha_.size();++i)
                {
                    field_type alphaPc = alpha_[i] * pC;
                    field_type temp = std::pow(alphaPc,vgN_[i]);
                    result+=fraction_[i]*-alpha_[i] * vgN_[i] * vgM_[i] *
                        temp/alphaPc * std::pow(1.0 + temp, -vgM_[i]-1.);
                }
                return(result);
            }

            field_type EffKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetValue(pC));
                else
                    return(krLPcTab_.GetValue(krLPcTab_.GetMaxX()));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetDerivative(pC));
                else
                    return(0.);
            }

            field_type EffKRelL(field_type sL)
            {
                return(EffKRelLpC(EffPc(sL)));
            }

            field_type EffDKRelL(field_type sL)
            {
                return(EffDPc(sL)*EffDKRelLpC(EffPc(sL)));
            }

            field_type EffKRelG(field_type sG)
            {
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return(0.);
            }
        };


        template<typename FT>
        class MultVanGenuchtenSplineParam : public MultVanGenuchtenParam<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            using MultVanGenuchtenParam<FT>::alpha_;
            using MultVanGenuchtenParam<FT>::vgN_;
            using MultVanGenuchtenParam<FT>::vgM_;
            using MultVanGenuchtenParam<FT>::fraction_;
            using MultVanGenuchtenParam<FT>::tau_;
            using MultVanGenuchtenParam<FT>::totPoreIntegral_;
            CubicInterpolationClass<field_type> krLPcTab_;

        public:
            void ReadInterpolationValuesKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationValues(filename);
            }

            void SetAlpha(const std::vector<field_type>& alpha)
            {
                alpha_=alpha;
                // for (typename std::vector<field_type>::iterator i=alpha_.begin();i!=alpha_.end();++i)
                //     (*i) /= (PhysicalChemistry<field_type>::DENSITY_WATER
                //         * PhysicalChemistry<field_type>::GRAVITY);
            }

            void SetFraction(const std::vector<field_type>& fraction)
            {
                fraction_=fraction;
                fraction_.push_back(0.0);
                field_type norm=0.;
                for (size_t i=0;i<fraction_.size()-1;++i)
                    norm+=fraction_[i];
                if (norm>=1.)
                {
                    fraction_.back()=0.;
                    for (size_t i=0;i<fraction_.size()-1;++i)
                        fraction_[i]/=norm;
                }
                else
                    fraction_.back()=1.-norm;
            }

            field_type EffKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetValue(pC));
                else
                    return(krLPcTab_.GetValue(krLPcTab_.GetMaxX()));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetDerivative(pC));
                else
                    return(0.);
            }

            field_type EffKRelL(field_type sL)
            {
                return(EffKRelLpC(EffPc(sL)));
            }

            field_type EffDKRelL(field_type sL)
            {
                return(EffDPc(sL)*EffDKRelLpC(EffPc(sL)));
            }
        };


        template<typename FT>
        class BrooksCoreyParam : public HydrParamBase<FT>
        {
            typedef FT field_type;
        protected:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            field_type lambda_;
            field_type tau_;

        public:
            BrooksCoreyParam() :
                lambda_(-1.),
                tau_(0.5)
            {}

            ~BrooksCoreyParam() {}

            void CheckValidity()
            {
                // Check if all parameters are set to a reasonable value
                if (lambda_ < 0)
                    DUNE_THROW(Dune::Exception,
                               "Please set parameter lambda>0 for the Brooks Corey model. (Use functionSetLambda)");
            }

            void SetLambda(field_type lambda)
            {
                assert(lambda>0.);
                lambda_=lambda;
            }

            void SetTau(field_type tau)
            {
                tau_=tau;
            }

            field_type EffSl(field_type pC)
            {
                ;
                return(std::pow(entryPress_/pC,lambda_));
            }

            field_type EffDSl(field_type pC)
            {
                ;
                return((-lambda_/entryPress_)*std::pow(entryPress_/pC,lambda_+1.));
            }

            field_type EffKRelLpC(field_type pC)
            {
                ;
                return(EffKRelL(EffSl(pC)));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                ;
                return(EffDSl(pC)*EffKRelL(EffSl(pC)));
            }

            field_type EffKRelG(field_type sG)
            {
                ;
                if (sG < std::numeric_limits<field_type>::epsilon())
                    sG = std::numeric_limits<field_type>::epsilon();
                return (std::pow(sG,tau_) * std::pow(1.0-std::pow(1.0-sG,((1.0+lambda_)/lambda_)),2.));
            }

            field_type EffKRelL(field_type sL)
            {
                ;
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return(std::pow(sL,2.0 + tau_ + (2.0/lambda_)));
            }

            field_type EffDKRelL(field_type sL)
            {
                ;
                if (sL < std::numeric_limits<field_type>::epsilon())
                    sL = std::numeric_limits<field_type>::epsilon();
                return ((2.0 + tau_ + (2.0/lambda_)) *
                    std::pow(sL,1.0 + tau_ + (2.0/lambda_)));
            }

            field_type EffPc(field_type sL)
            {
                ;
                return(entryPress_*std::pow(sL,-1./lambda_));
            }

            field_type EffDPc(field_type sL)
            {
                ;
                return(-entryPress_/lambda_*std::pow(sL,-1./lambda_-1.));
            }
        };


        template<typename FT>
        class RegBrooksCoreyParam : public BrooksCoreyParam<FT>
        {
            typedef FT field_type;
        private:
            // import base class members
            using BrooksCoreyParam<FT>::resSatL_;
            using BrooksCoreyParam<FT>::resSatG_;
            using BrooksCoreyParam<FT>::availableWater_;
            using BrooksCoreyParam<FT>::entryPress_;
            using BrooksCoreyParam<FT>::delta_;
            using BrooksCoreyParam<FT>::lambda_;
            using BrooksCoreyParam<FT>::tau_;

            struct RegParType
            {
                field_type x0;
                field_type x1;
                std::vector<field_type> coeff;
            };

            field_type interval_;
            RegParType regPar_;

            void CheckValidity()
            {
                // Check if necessary functions were called
                if (regPar_.x1== 0)
                    DUNE_THROW(Dune::Exception,
                               "regPar_.x1 of the regularized Brooks Corey model is 0. Did you forget to call InitRegularisation()?");
            }

            void CompSmoothParam(std::vector<field_type> &regPar,field_type x0,field_type x1,field_type g0,field_type dg0,field_type g1,field_type dg1)
            {
                field_type t1,t2,t3,t4,t5,t6,t7,t8;
                field_type xa=(x0+x1)/2.;
                t1 = x0*x0;
                t2 = t1*xa;
                t3 = g0*xa;
                t4 = x1*x1;
                t5 = xa*x0;
                t6 = t4*xa;
                t7 = t4*x0;
                t8 = g1*xa;
                regPar[0] = (-x1*dg1-dg0*x1+xa*dg1+2.0*x0*dg0+2.0*g1-2.0*g0-xa*dg0)/(-xa+x0)/(-x1+x0)/2.0;
                regPar[1] = -(-x1*x0*dg1+xa*x0*dg1+t1*dg0+2.0*x0*g1-2.0*g0*x0-dg0*xa*x1)/(-xa+x0)/(-x1+x0);
                regPar[2] = (-t1*x1*dg1+t1*dg0*x1+t2*dg1+2.0*t1*g1+t2*dg0-2.0*x0*dg0*xa*x1+2.0*t3*x1-2.0*t3*x0-2.0*g0*x0*x1)/(-xa+x0)/(-x1+x0)/2.0;
                regPar[3] = (-xa*dg0+xa*dg1+x0*dg0+x0*dg1-2.0*g0-2.0*x1*dg1+2.0*g1)/(xa*x1-xa*x0+x0*x1-t4)/2.0;
                regPar[4] = -(-dg0*xa*x1+dg0*x0*x1-2.0*x1*g0-t4*dg1+2.0*x1*g1+t5*dg1)/(xa*x1-t5+x0*x1-t4);
                regPar[5] = (-t6*dg0-t6*dg1+t7*dg0-t7*dg1-2.0*t4*g0+2.0*x1*dg1*t5+2.0*t8*x1-2.0*t8*x0+2.0*g1*x0*x1)/(xa*x1-t5+x0*x1-t4)/2.0;
            }

            void CompSmooth()
            {
                field_type g0,g1,dg0,dg1;
                field_type x0,x1;
                field_type smoothIntervall;

                smoothIntervall=this->Pc(1.-resSatG_-interval_)-entryPress_;
                x0 = entryPress_;
                x1 = entryPress_+smoothIntervall;
                g0 = 1.;
                dg0 = this->DSl(x0);
                g1 = BrooksCoreyParam<FT>::EffSl(x1);
                dg1 = BrooksCoreyParam<FT>::EffDSl(x1);
                CompSmoothParam(regPar_.coeff,x0,x1,g0,dg0,g1,dg1);
                if ((2.*regPar_.coeff[0]*(x0+smoothIntervall*0.5)+regPar_.coeff[1])<0.)
                {
                    regPar_.x1 = entryPress_+smoothIntervall;
                    regPar_.x0 = entryPress_+0.5*smoothIntervall;
                }
                else
                {
                    std::cerr << "Deactivating smoothing pressure/saturation curve" << std::endl;
                    for (unsigned int i=0;i<6;i++)
                        std::cerr << regPar_.coeff[i] << "  ";
                    std::cerr << std::endl;
                    regPar_.x1 = 0.;
                }
            }

        public:
            RegBrooksCoreyParam() : interval_(0.001)
            {
                regPar_.x1=0.;
                regPar_.coeff.resize(6);
            }

            ~RegBrooksCoreyParam() {}

            void InitRegularisation(field_type interval=0.001)
            {
                interval_=interval;
                CompSmooth();
            }

            field_type EffSl(field_type pC)
            {
                ;
                if (pC <= regPar_.x1)
                    {
                        if (pC <= regPar_.x0)
                            return(regPar_.coeff[0]*pC*pC+regPar_.coeff[1]*pC+regPar_.coeff[2]);
                        else
                            return(regPar_.coeff[3]*pC*pC+regPar_.coeff[4]*pC+regPar_.coeff[5]);
                    }
                return(std::pow(entryPress_/pC,lambda_));
            }

            field_type EffDSl(field_type pC)
            {
                ;
                if (pC <= regPar_.x1)
                {
                    if (pC <= regPar_.x0)
                        return(2.*regPar_.coeff[0]*pC+regPar_.coeff[1]);
                    else
                        return(2.*regPar_.coeff[3]*pC+regPar_.coeff[4]);
                }
                return((-lambda_/entryPress_)*std::pow(entryPress_/pC,lambda_+1.));
            }
        };

        template<typename FT>
        class HydrLinearInterpolation : public HydrParamBase<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            HydrParamBase<FT> *hydrParamObject_;
            field_type maxPc_;
            field_type minSl_;
            LinearInterpolationClass<field_type> sLTab_;
            LinearInterpolationClass<field_type> krLTab_;
            LinearInterpolationClass<field_type> krGTab_;
            LinearInterpolationClass<field_type> krLPcTab_;
            LinearInterpolationClass<field_type> pCTab_;
            bool initset_;

            void CheckValidity()
            {
                // Check if necessary functions were called
                if (initset_ == false)
                    DUNE_THROW(Dune::Exception,
                               "initset_ of the HydrLinearInterpolation class is false. Did you forget to call Init()?");
            }

            void ReadInterpolationValuesSl(std::string filename)
            {
                sLTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrL(std::string filename)
            {
                krLTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrG(std::string filename)
            {
                krGTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesPc(std::string filename)
            {
                pCTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationTableSl(std::string filename)
            {
                sLTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrL(std::string filename)
            {
                krLTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrG(std::string filename)
            {
                krGTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTablePc(std::string filename)
            {
                pCTab_.ReadInterpolationTable(filename);
            }
            void WriteInterpolationTableSl(std::string filename)
            {
                sLTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrL(std::string filename)
            {
                krLTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrG(std::string filename)
            {
                krGTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTablePc(std::string filename)
            {
                pCTab_.WriteInterpolationTable(filename);
            }
            void InterpolateSl(int numPoints)
            {
                std::vector<field_type> y(numPoints);
                field_type minX=hydrParamObject_->GetEntryPressure();
                field_type maxX=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type interval=(maxX-minX)/(numPoints-1.);
                y[0]=1.;
                for (int i=1;i<numPoints;i++)
                    y[i]=hydrParamObject_->EffSl(minX+i*interval);
                sLTab_.Init(minX,maxX,y);
            }
            void InterpolatePc(int numPoints)
            {
                std::vector<field_type> y(numPoints);
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->EffSl(pc);
                field_type maxX=1.0;
                field_type interval=(maxX-minX)/(numPoints-1.);
                for (int i=0;i<numPoints-1;i++)
                    y[i]=hydrParamObject_->EffPc(minX+i*interval);
                y[numPoints-1]=hydrParamObject_->GetEntryPressure();
                pCTab_.Init(minX,maxX,y);
            }
            void InterpolateKrLPc(int numPoints)
            {
                std::vector<field_type> y(numPoints);
                field_type minX=hydrParamObject_->GetEntryPressure();
                field_type maxX=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type interval=(maxX-minX)/(numPoints-1.);
                y[0]=1.;
                for (int i=1;i<numPoints;i++)
                    y[i]=hydrParamObject_->EffKRelLpC(minX+i*interval);
                krLPcTab_.Init(minX,maxX,y);
            }
            void InterpolateKrL(int numPoints)
            {
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->EffSl(pc);
                field_type maxX=1.0;
                numPoints=std::min(numPoints,int((maxX-minX)/minSl_)+1);
                std::vector<field_type> y(numPoints);
                field_type interval=(maxX-minX)/(numPoints-1.);
                for (int i=0;i<numPoints-1;i++)
                    y[i]=hydrParamObject_->EffKRelL(minX+i*interval);
                y[numPoints-1]=1.;
                krLTab_.Init(minX,maxX,y);
            }
            void Interpolate(int numPoints, int interpolationSwitch)
            {
                entryPress_=hydrParamObject_->GetEntryPressure();
                delta_=hydrParamObject_->Getdelta();
                resSatL_=hydrParamObject_->GetResSatL();
                resSatG_=hydrParamObject_->GetResSatG();
                availableWater_=hydrParamObject_->GetAvailableWater();
                std::vector<field_type> pcT,krlT,krlPcT,slT;
                field_type minX=1.-int(((1.-minSl_)/minSl_)+0.5)*minSl_;
                field_type maxX=1.0;
                int numIntervalls=int((maxX-minX)/minSl_+0.5);
                field_type interval=maxPc_;
                pcT.push_back(hydrParamObject_->EffPc(minX));
                if (interpolationSwitch&KRL)
                    krlT.push_back(hydrParamObject_->EffKRelL(minX));
                field_type newPc;
                for (int i=1;i<numIntervalls;++i)
                {
                    double sL=minX+i*minSl_;
                    newPc=hydrParamObject_->EffPc(sL);
                    interval=std::min(interval,pcT.back()-newPc);
                    pcT.push_back(newPc);
                    if (interpolationSwitch&KRL)
                        krlT.push_back(hydrParamObject_->EffKRelL(sL));
                }
                newPc=hydrParamObject_->GetEntryPressure();
                interval=std::min(interval,pcT.back()-newPc);
                pcT.push_back(newPc);
                pCTab_.Init(minX,maxX,pcT);
                if (interpolationSwitch&KRL)
                {
                    krlT.push_back(1.);
                    krLTab_.Init(minX,maxX,krlT);
                }
                minX=pcT.back();
                maxX=std::min(maxPc_,pcT.front());
                if (((numPoints-1)*interval+minX)<maxX)
                    maxX=(numPoints-1)*interval+minX;
                else
                {
                    numPoints=int((maxX-minX)/interval + 0.5)+1;
                    interval=(maxX-minX)/(numPoints-1.);
                }
#ifdef DEBUG_INTERPOLATION_
                std::std::cout << left << scientific << setw(10) << setprecision(5) << "interpolate with resolution " << minSl_ << " using an interval of " << interval << " Pa times " << numPoints << " points resulting in a upper limit for capillary pressure of " << maxX << " Pa at a water saturation of " << hydrParamObject_->EffSl(maxX) << std::std::endl;
#endif
                slT.push_back(1.);
                if (interpolationSwitch&KRLPC)
                    krlPcT.push_back(1.);
                for (int i=1;i<numPoints;i++)
                {
                    slT.push_back(hydrParamObject_->EffSl(minX+i*interval));
                    if (interpolationSwitch&KRLPC)
                        krlPcT.push_back(hydrParamObject_->EffKRelLpC(minX+i*interval));
                }
                sLTab_.Init(minX,maxX,slT);
                if (interpolationSwitch&KRLPC)
                    krLPcTab_.Init(minX,maxX,krlPcT);
            }
            void InterpolateKrG(int numPoints)
            {
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=0.;
                field_type maxX=1.-hydrParamObject_->EffSl(pc);
                numPoints=std::min(numPoints,int((maxX-minX)/minSl_)+1);
                std::vector<field_type> y(numPoints);
                field_type interval=(maxX-minX)/(numPoints-1.);
                for (int i=0;i<numPoints-1;i++)
                    y[i]=hydrParamObject_->EffKRelG(minX+i*interval);
                y[numPoints-1]=1.;
                krGTab_.Init(minX,maxX,y);
            }

        public:
            enum interpolationSwitchType
            {
                SL=1,KRL=2,KRG=4,KRLPC=8,PC=16
            };

            HydrLinearInterpolation() :
                hydrParamObject_(0),
                initset_(false)
            {}

            HydrLinearInterpolation(HydrParamBase<FT> *hydrParamObject) :
                hydrParamObject_(hydrParamObject),
                maxPc_(1e5),
                minSl_(0.001),
                initset_(false)
            {}

            ~HydrLinearInterpolation()
            {
                if (hydrParamObject_!=0)
                    delete hydrParamObject_;
            }

            void SetMaxPc(field_type value)
            {
                maxPc_=value;
            }

            void SetMinSl(field_type value)
            {
                minSl_=value;
            }

            void Init(int numPoints=10000)
            {
                Interpolate(numPoints,SL|KRL|KRLPC|PC);
                InterpolateKrG(numPoints);
                initset_ = true;
            }

            void Init(int interpolationSwitch, int numPoints)
            {
                Interpolate(numPoints,interpolationSwitch);
                if (interpolationSwitch&KRG)
                    InterpolateKrG(numPoints);
                initset_ = true;
            }

            field_type EffSl(field_type pC)
            {
                ;
                if (pC<=sLTab_.GetMaxX())
                    return(sLTab_.GetValue(pC));
                else
                    return(hydrParamObject_->EffSl(pC));
            }

            field_type EffDSl(field_type pC)
            {
                ;
                if (pC<=sLTab_.GetMaxX())
                    return(sLTab_.GetDerivative(pC));
                else
                    return(hydrParamObject_->EffDSl(pC));
            }

            field_type EffKRelLpC(field_type pC)
            {
                ;
                if (pC<=krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetValue(pC));
                else
                    return(hydrParamObject_->EffKRelLpC(pC));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                ;
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetDerivative(pC));
                else
                    return(hydrParamObject_->EffDKRelLpC(pC));
            }

            field_type EffKRelG(field_type sG)
            {
                ;
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetValue(sG));
                else
                    return(hydrParamObject_->EffKRelG(sG));
            }

            field_type EffDKRelG(field_type sG)
            {
                ;
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetDerivative(sG));
                else
                    return(hydrParamObject_->EffDKRelG(sG));
            }

            field_type EffKRelL(field_type sL)
            {
                ;
                if (sL>=krLTab_.GetMinX())
                    return(krLTab_.GetValue(sL));
                else
                    return(hydrParamObject_->EffKRelL(sL));
            }

            field_type EffDKRelL(field_type sL)
            {
                ;
                if (sL>=krLTab_.GetMinX())
                    return(krLTab_.GetDerivative(sL));
                else
                    return(hydrParamObject_->EffDKRelL(sL));
            }

            field_type EffPc(field_type sL)
            {
                ;
                if (sL>=pCTab_.GetMinX())
                    return(pCTab_.GetValue(sL));
                else
                    return(hydrParamObject_->EffPc(sL));
            }

            field_type EffDPc(field_type sL)
            {
                ;
                if (sL>=pCTab_.GetMinX())
                    return(pCTab_.GetDerivative(sL));
                else
                    return(hydrParamObject_->EffDPc(sL));
            }
        };


        template<typename FT>
        class HydrCubicInterpolation : public HydrParamBase<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            HydrParamBase<FT> *hydrParamObject_;
            field_type maxPc_;
            field_type minSl_;
            CubicInterpolationClass<field_type> sLTab_;
            CubicInterpolationClass<field_type> krLTab_;
            CubicInterpolationClass<field_type> krGTab_;
            CubicInterpolationClass<field_type> krLPcTab_;
            CubicInterpolationClass<field_type> pCTab_;

            void InterpolateSl(int numPoints)
            {
                std::vector<field_type> x(numPoints),y(numPoints);
                field_type maxX=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->GetEntryPressure();
                field_type interval=(maxX-minX)/(numPoints-1.);
                x[0]=hydrParamObject_->GetEntryPressure();
                y[0]=1.;
                for (int i=1;i<(numPoints-1);i++)
                {
                    x[i]=x[0]+i*interval;
                    y[i]=hydrParamObject_->EffSl(x[i]);
                }
                sLTab_.Init(x,y,hydrParamObject_->EffDSl(x.front()+std::numeric_limits<field_type>::epsilon()),hydrParamObject_->EffDSl(x.back()));
                pCTab_.Init(y,x,hydrParamObject_->DPc(y.front()*availableWater_+resSatL_),hydrParamObject_->DPc(y.back()*availableWater_+resSatL_));
            }
            void InterpolateKrL(int numPoints)
            {
                std::vector<field_type> x(numPoints),y(numPoints);
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->EffSl(pc);
                field_type maxX=1.0;
                field_type interval=(maxX-minX)/(numPoints-1.);
                for (int i=0;i<numPoints-1;i++)
                {
                    x[i]=minX+interval*i;
                    y[i]=hydrParamObject_->EffKRelL(minX+i*interval);
                }
                x[numPoints-1]=1.;
                y[numPoints-1]=1.;
                krLTab_.Init(x,y,hydrParamObject_->EffDKRelL(x.front()+std::numeric_limits<field_type>::epsilon()),hydrParamObject_->EffDKRelL(x.back()));
            }
            void InterpolateKrG(int numPoints)
            {
                std::vector<field_type> x(numPoints),y(numPoints);
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=0.;
                field_type maxX=1.-hydrParamObject_->EffSl(pc);
                field_type interval=(maxX-minX)/(numPoints-1.);
                for (int i=0;i<numPoints-1;i++)
                {
                    x[i]=minX+interval*i;
                    y[i]=hydrParamObject_->EffKRelG(minX+i*interval);
                }
                x[numPoints-1]=1.;
                y[numPoints-1]=1.;
                krGTab_.Init(x,y,hydrParamObject_->EffDKRelG(x.front()+std::numeric_limits<field_type>::epsilon()),hydrParamObject_->EffDKRelG(x.back()));
            }
            void InterpolateKrLPc(int numPoints)
            {
                std::vector<field_type> x(numPoints),y(numPoints);
                field_type maxX=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->GetEntryPressure();
                field_type interval=(maxX-minX)/(numPoints-1.);
                x[0]=hydrParamObject_->GetEntryPressure();
                y[0]=1.;
                for (int i=1;i<numPoints;i++)
                {
                    x[i]=x[0]+i*interval;
                    y[i]=hydrParamObject_->EffKRelLpC(x[i]);
                }
                krLPcTab_.Init(x,y,hydrParamObject_->EffDKRelLpC(x.front()+std::numeric_limits<field_type>::epsilon()),hydrParamObject_->EffDKRelLpC(x.back()));
            }
            void InterpolatePc(int numPoints)
            {
                std::vector<field_type> x(numPoints),y(numPoints);
                field_type pc=std::min(maxPc_,hydrParamObject_->EffPc(minSl_));
                field_type minX=hydrParamObject_->EffSl(pc);
                field_type maxX=1.0;
                field_type interval=(maxX-minX)/(numPoints-1.);
                x[0]=minX;
                y[0]=hydrParamObject_->GetEntryPressure();
                for (int i=1;i<numPoints;i++)
                {
                    x[i]=minX+interval*i;
                    y[i]=hydrParamObject_->EffPc(minX+i*interval);
                }
                pCTab_.Init(x,y,hydrParamObject_->DPc(x.front()*availableWater_+resSatL_),hydrParamObject_->DPc(x.back()*availableWater_+resSatL_));
            }

        public:
            enum interpolationSwitchType
            {
                SL=1,KRL=2,KRG=4,KRLPC=8,PC=16
            };

            HydrCubicInterpolation() :
                hydrParamObject_(0),
                maxPc_(1e5)
            {}

            HydrCubicInterpolation(HydrParamBase<FT> *hydrParamObject) :
                hydrParamObject_(hydrParamObject),
                maxPc_(1e5),
                minSl_(0.001)
            {}

            ~HydrCubicInterpolation()
            {
                if (hydrParamObject_!=0)
                    delete hydrParamObject_;
            }

            //! \todo doc me!
            void SetMaxPc(field_type value)
            {
                maxPc_=value;
            }

            //! \todo doc me!
            void SetMinSl(field_type value)
            {
                minSl_=value;
            }

            //! \todo doc me!
            void ReadInterpolationTableSl(std::string filename)
            {
                sLTab_.ReadInterpolationTable(filename);
            }
            //! \todo doc me!
            void ReadInterpolationTableKrL(std::string filename)
            {
                krLTab_.ReadInterpolationTable(filename);
            }
            //! \todo doc me!
            void ReadInterpolationTableKrG(std::string filename)
            {
                krGTab_.ReadInterpolationTable(filename);
            }
            //! \todo doc me!
            void ReadInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationTable(filename);
            }
            //! \todo doc me!
            void ReadInterpolationTablePc(std::string filename)
            {
                pCTab_.ReadInterpolationTable(filename);
            }
            //! \todo doc me!
            void WriteInterpolationTableSl(std::string filename)
            {
                sLTab_.WriteInterpolationTable(filename);
            }
            //! \todo doc me!
            void WriteInterpolationTableKrL(std::string filename)
            {
                krLTab_.WriteInterpolationTable(filename);
            }
            //! \todo doc me!
            void WriteInterpolationTableKrG(std::string filename)
            {
                krGTab_.WriteInterpolationTable(filename);
            }
            //! \todo doc me!
            void WriteInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.WriteInterpolationTable(filename);
            }
            //! \todo doc me!
            void WriteInterpolationTablePc(std::string filename)
            {
                pCTab_.WriteInterpolationTable(filename);
            }
            //! \todo doc me!
            void ReadInterpolationValuesSl(std::string filename)
            {
                sLTab_.ReadInterpolationValues(filename);
            }
            //! \todo doc me!
            void ReadInterpolationValuesKrL(std::string filename)
            {
                krLTab_.ReadInterpolationValues(filename);
            }
            //! \todo doc me!
            void ReadInterpolationValuesKrG(std::string filename)
            {
                krGTab_.ReadInterpolationValues(filename);
            }
            //! \todo doc me!
            void ReadInterpolationValuesKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationValues(filename);
            }
            //! \todo doc me!
            void ReadInterpolationValuesPc(std::string filename)
            {
                pCTab_.ReadInterpolationValues(filename);
            }
            //! \todo doc me!
            typedef int INT;
            //! \todo doc me!
            void SetInterpolationValuesSl(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        sLTab_.Init(x,y);
                        pCTab_.Init(x,y);
                        break;
                    case 1:
                        sLTab_.Init(x,y,firstDerivative);
                        if (firstDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(x,y,1./firstDerivative);
                        else
                            pCTab_.Init(x,y,1e100);
                        break;
                    case 2:
                        sLTab_.InitLast(x,y,lastDerivative);
                        if (lastDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.InitLast(x,y,1./lastDerivative);
                        else
                            pCTab_.InitLast(x,y,1e100);
                        break;
                    case 3:
                        sLTab_.Init(x,y,firstDerivative,lastDerivative);
                        if ((firstDerivative>std::numeric_limits<field_type>::epsilon())&&(lastDerivative>std::numeric_limits<field_type>::epsilon()))
                            pCTab_.Init(x,y,1./firstDerivative,1./lastDerivative);
                        if (firstDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(x,y,1./firstDerivative,1e100);
                        else if (lastDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(x,y,1e100,1./lastDerivative);
                        else
                            pCTab_.Init(x,y,1e100,1e100);
                        break;
                }
            }
            void SetInterpolationValuesKrl(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krLTab_.Init(x,y);
                        break;
                    case 1:
                        krLTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krLTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krLTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }
            void SetInterpolationValuesKrLPc(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krLPcTab_.Init(x,y);
                        break;
                    case 1:
                        krLPcTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krLPcTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krLPcTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }
            void SetInterpolationValuesKrG(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krGTab_.Init(x,y);
                        break;
                    case 1:
                        krGTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krGTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krGTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }

            void Init(int numPoints=100)
            {
                entryPress_=hydrParamObject_->GetEntryPressure();
                delta_=hydrParamObject_->Getdelta();
                resSatL_=hydrParamObject_->GetResSatL();
                resSatG_=hydrParamObject_->GetResSatG();
                availableWater_=hydrParamObject_->GetAvailableWater();
                InterpolateSl(numPoints);
                InterpolateKrL(numPoints);
                InterpolateKrG(numPoints);
                InterpolateKrLPc(numPoints);
                InterpolatePc(numPoints);
            }

            void Init(int interpolationSwitch, int numPoints)
            {
                if (interpolationSwitch&SL)
                    InterpolateSl(numPoints);
                if (interpolationSwitch&KRL)
                    InterpolateKrL(numPoints);
                if (interpolationSwitch&KRLPC)
                    InterpolateKrLPc(numPoints);
                if (interpolationSwitch&PC)
                    InterpolatePc(numPoints);
                if (interpolationSwitch&KRG)
                    InterpolateKrG(numPoints);
            }

            field_type EffSl(field_type pC)
            {
                if (pC<=sLTab_.GetMaxX())
                    return(sLTab_.GetValue(pC));
                else
                    return(hydrParamObject_->EffSl(pC));
            }

            field_type EffDSl(field_type pC)
            {
                if (pC<=sLTab_.GetMaxX())
                    return(sLTab_.GetDerivative(pC));
                else
                    return(hydrParamObject_->EffDSl(pC));
            }

            field_type EffKRelLpC(field_type pC)
            {
                if (pC<=krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetValue(pC));
                else
                    return(hydrParamObject_->EffKRelLpC(pC));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetDerivative(pC));
                else
                    return(hydrParamObject_->EffDKRelLpC(pC));
            }

            field_type EffKRelG(field_type sG)
            {
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetValue(sG));
                else
                    return(hydrParamObject_->EffKRelG(sG));
            }

            field_type EffDKRelG(field_type sG)
            {
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetDerivative(sG));
                else
                    return(hydrParamObject_->EffDKRelG(sG));
            }

            field_type EffKRelL(field_type sL)
            {
                if (sL>=krLTab_.GetMinX())
                    return(krLTab_.GetValue(sL));
                else
                    return(hydrParamObject_->EffKRelL(sL));
            }

            field_type EffDKRelL(field_type sL)
            {
                if (sL>=krLTab_.GetMinX())
                    return(krLTab_.GetDerivative(sL));
                else
                    return(hydrParamObject_->EffDKRelL(sL));
            }

            field_type EffPc(field_type sL)
            {
                if (sL>=pCTab_.GetMinX())
                    return(pCTab_.GetValue(sL));
                else
                    return(hydrParamObject_->EffPc(sL));
            }

            field_type EffDPc(field_type sL)
            {
                if (sL>=pCTab_.GetMinX())
                    return(pCTab_.GetDerivative(sL));
                else
                    return(hydrParamObject_->EffDPc(sL));
            }
        };

        template<typename FT>
        class HydrSpline : public HydrParamBase<FT>
        {
            typedef FT field_type;
        private:
            using HydrParamBase<FT>::resSatL_;
            using HydrParamBase<FT>::resSatG_;
            using HydrParamBase<FT>::availableWater_;
            using HydrParamBase<FT>::entryPress_;
            using HydrParamBase<FT>::delta_;
            CubicInterpolationClass<field_type> sLTab_;
            CubicInterpolationClass<field_type> krLTab_;
            CubicInterpolationClass<field_type> krGTab_;
            CubicInterpolationClass<field_type> krLPcTab_;
            CubicInterpolationClass<field_type> pCTab_;

        public:
            void ReadInterpolationTableSl(std::string filename)
            {
                sLTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrL(std::string filename)
            {
                krLTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrG(std::string filename)
            {
                krGTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTablePc(std::string filename)
            {
                pCTab_.ReadInterpolationTable(filename);
            }
            void WriteInterpolationTableSl(std::string filename)
            {
                sLTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrL(std::string filename)
            {
                krLTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrG(std::string filename)
            {
                krGTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableKrLPc(std::string filename)
            {
                krLPcTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTablePc(std::string filename)
            {
                pCTab_.WriteInterpolationTable(filename);
            }
            void ReadInterpolationValuesSl(std::string filename)
            {
                sLTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrL(std::string filename)
            {
                krLTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrG(std::string filename)
            {
                krGTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesKrLPc(std::string filename)
            {
                krLPcTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesPc(std::string filename)
            {
                pCTab_.ReadInterpolationValues(filename);
            }
            typedef int INT;
            void SetInterpolationValuesSl(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        sLTab_.Init(x,y);
                        pCTab_.Init(y,x);
                        break;
                    case 1:
                        sLTab_.Init(x,y,firstDerivative);
                        if (firstDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(y,x,1./firstDerivative);
                        else
                            pCTab_.Init(y,x,1e100);
                        break;
                    case 2:
                        sLTab_.InitLast(x,y,lastDerivative);
                        if (lastDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.InitLast(y,x,1./lastDerivative);
                        else
                            pCTab_.InitLast(y,x,1e100);
                        break;
                    case 3:
                        sLTab_.Init(x,y,firstDerivative,lastDerivative);
                        if ((firstDerivative>std::numeric_limits<field_type>::epsilon())&&(lastDerivative>std::numeric_limits<field_type>::epsilon()))
                            pCTab_.Init(y,x,1./firstDerivative,1./lastDerivative);
                        if (firstDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(y,x,1./firstDerivative,1e100);
                        else if (lastDerivative>std::numeric_limits<field_type>::epsilon())
                            pCTab_.Init(y,x,1e100,1./lastDerivative);
                        else
                            pCTab_.Init(y,x,1e100,1e100);
                        break;
                }
            }
            void SetInterpolationValuesKrl(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krLTab_.Init(x,y);
                        break;
                    case 1:
                        krLTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krLTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krLTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }
            void SetInterpolationValuesKrLPc(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krLPcTab_.Init(x,y);
                        break;
                    case 1:
                        krLPcTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krLPcTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krLPcTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }
            void SetInterpolationValuesKrG(const std::vector<field_type> &x, const std::vector<field_type> &y, const INT derTp, const field_type firstDerivative=0., const field_type lastDerivative=0.)
            {
                switch(derTp)
                {
                    case 0:
                        krGTab_.Init(x,y);
                        break;
                    case 1:
                        krGTab_.Init(x,y,firstDerivative);
                        break;
                    case 2:
                        krGTab_.InitLast(x,y,lastDerivative);
                        break;
                    case 3:
                        krGTab_.Init(x,y,firstDerivative,lastDerivative);
                        break;
                }
            }

            field_type EffSl(field_type pC)
            {
                if (pC<sLTab_.GetMaxX())
                    return(sLTab_.GetValue(pC));
                else
                    return(sLTab_.GetValue(sLTab_.GetMaxX()));
            }

            field_type EffDSl(field_type pC)
            {
                if (pC<sLTab_.GetMaxX())
                    return(sLTab_.GetDerivative(pC));
                else
                    return(0.);
            }

            field_type EffKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetValue(pC));
                else
                    return(krLPcTab_.GetValue(krLPcTab_.GetMaxX()));
            }

            field_type EffDKRelLpC(field_type pC)
            {
                if (pC<krLPcTab_.GetMaxX())
                    return(krLPcTab_.GetDerivative(pC));
                else
                    return(0.);
            }

            field_type EffKRelG(field_type sG)
            {
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetValue(sG));
                else
                    return(krGTab_.GetValue(krGTab_.GetMinX()));
            }

            field_type EffDKRelG(field_type sG)
            {
                if (sG>krGTab_.GetMinX())
                    return(krGTab_.GetDerivative(sG));
                else
                    return(0.);
            }

            field_type EffKRelL(field_type sL)
            {
                if (sL>krLTab_.GetMinX())
                    return(krLTab_.GetValue(sL));
                else
                    return(krLTab_.GetValue(krLTab_.GetMinX()));
            }

            field_type EffDKRelL(field_type sL)
            {
                if (sL>krLTab_.GetMinX())
                    return(krLTab_.GetDerivative(sL));
                else
                    return(0.);
            }

            field_type EffPc(field_type sL)
            {
                if (sL>pCTab_.GetMinX())
                    return(pCTab_.GetValue(sL));
                else
                    return(pCTab_.GetValue(pCTab_.GetMinX()));
            }

            field_type EffDPc(field_type sL)
            {
                if (sL>pCTab_.GetMinX())
                    return(pCTab_.GetDerivative(sL));
                else
                    return(0.);
            }
        };

    } // end namespace PM
} // end namespace Dune

#endif // DUNE_PM_HYDRPAR_HH
