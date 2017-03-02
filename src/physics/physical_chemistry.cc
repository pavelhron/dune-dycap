#include <cmath>
#include <limits>
#include "physical_chemistry.hh"
#include "interpolation.h"
#include<dune/common/exceptions.hh>


namespace Dune {
    namespace PM {

        /**
           \brief compute the viscosity of water at given temperature.

           It uses a function given in P. W. Atkins: Physikalische
           Chemie, 2. korr. Nachdruck der 1. Auflage, VCH
           Verlagsgesellschaft, 1990 S. 872, Tabelle 24-5.

           \param temp temperature
           \return viscosity
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::ViscosityWater(FLOAT temp)
        {
            FLOAT dum;
            if (temp<273.15)
                temp=273.15;
            dum = temp - 293.15;
            return (1.0019e-3*std::exp(-M_LN10*(1.37023*dum+8.36e-4*dum*dum)/(temp-164.15)));  /* [kg/m/s] */
        }


        /**
           \brief compute the viscosity of air at given temperature.

           The values are linear interpolated from values given in
           P. W. Atkins: Physikalische Chemie, 2. korr. Nachdruck der
           1. Auflage, VCH Verlagsgesellschaft, 1990, S. 873, Tabelle
           26-3.

           \param temp temperature
           \return viscosity
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::ViscosityAir(FLOAT temp)
        {
            return (1.73e-5+0.0045e-5*(temp-273.0));
        }


        /**
           \brief compute the surface tension of water at given temperature
           relative to the surface tension at 20 °C.

           It uses a function given ???

           \param temp  temperature
           \return relative surface tension
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::RelSurfaceTensionWaterAir(FLOAT temp)
        {
            return((116.19 - temp * 0.1477)/72.89); // surface tension is 72.89 mN/m at 20° C
        }


        /**
           \brief compute the molar density of water vapor at a given
           temperature, capillary pressure and liquid phase
           composition

           \param pC    capillary pressure
           \param temp  temperature
           \param XWL   molar fraction of water
           \return molar density
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::MolarDensityVapor (FLOAT pC, FLOAT temp, FLOAT XWL)
        {
            return (XWL*(610.78/(IDEAL_GAS_CONSTANT*temp)))*std::exp((17.08085*(temp-273.15))/(temp-38.975)+(-pC*MOLAR_VOLUME_WATER)/(IDEAL_GAS_CONSTANT*temp));  // [mol/m^3]
        }


        /**
           \brief compute the molar density of water vapor at a given
           temperature, and liquid phase composition.

           The gas pressure is the sum of the water vapor pressure and
           the partial pressure of dry air.

           \param pL    liquid phase pressure
           \param pA    partial pressure of air
           \param temp  temperature
           \param XWL   molar fraction of water
           \return molar density
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::MolarDensityVaporVT (FLOAT pL, FLOAT pA, FLOAT temp, FLOAT XWL)
        {
            FLOAT pVap = XWL*(610.78* std::exp((17.08085*(temp-273.15))/(temp-38.975)+((pL-pA)*MOLAR_VOLUME_WATER)/(IDEAL_GAS_CONSTANT*temp)));  // [mol/m^3]
            FLOAT a = MOLAR_VOLUME_WATER/(IDEAL_GAS_CONSTANT*temp);
            FLOAT b = pVap;
            FLOAT delta;
            do
                {
                    delta = (b*std::exp(-a*pVap)-pVap)/(a*b*std::exp(-a*pVap)+1.0);
                    pVap=pVap+delta;
                } while (delta>VAPOR_LIMIT);
            return(pVap/(IDEAL_GAS_CONSTANT*temp));
        }


        /**
           \brief compute mass density of water at a given
           temperature.

           The polynom is interpolated to values given by Kuchling:
           Handbuch der Physik, 13. korr. Auflage, Thun Verlag, 1991,
           S. 596, Tabelle 15 for temperatures above 0° C and a
           function given by Hare and Sorensen, J. Chem. Phys. 87:
           4840-4845 below 0° C.

           \param temp temperature
           \return density
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::DensityWater(FLOAT temp)
        {
            FLOAT result;
            if (temp<273.15)
                {
                    result = -1.475e-2*(temp-273.15);
                    result = 999.8+(result + 1.584e-2)*(temp-273.15);
                }
            else
                {
                    result = -3.62471e-7*temp;
                    result = (result+4.65374e-04)*temp;
                    result = (result-0.227669)*temp;
                    result = (result+49.8233)*temp-3089.35;
                }
            return(result);
        }


        /**
           \brief compute molar density of water at a given
           temperature.

           The polynom is interpolated to values given by Kuchling:
           Handbuch der Physik, 13. korr. Auflage, Thun Verlag, 1991,
           S. 596, Tabelle 15 for temperatures above 0° C and a
           function given by Hare and Sorensen, J. Chem. Phys. 87:
           4840-4845 below 0° C. The coefficients for the mass density
           are divided by the molar mass.

           \param temp - temperature
           \return density
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::MolarDensityWater(FLOAT temp)
        {
            return(PhysicalChemistry<FLOAT>::DensityWater(temp)/MOLAR_MASS_WATER);
        }


        /**
           \brief compute the capillary pressure between an ice phase
           and a liquid water phase at a given temperature and ice
           pressure.

           The function is given in Spaans, E. J. A. and J. M. Baker:
           "The Soil Freezing Characteristic: Its Measurement and
           Similarity to the Soil Moisture Characteristic", SSSAJ 60:
           13-19 (1996) with a modification for the ice pressure by
           the author of this program.

           \param temp temperature
           \param pIce ice pressure
           \return freezing pressure
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::FreezingPressure(FLOAT temp,FLOAT pIce)
        {
            FLOAT result=(712380.*log(temp/273.15)-5545.*(temp-273.15)+3.14*(temp*temp-273.15*273.15)); /* [J/kg] */
            result*=PhysicalChemistry<FLOAT>::DensityWater(temp); /* [J/m3 = Pa] */
            return (result-(PhysicalChemistry<FLOAT>::DensityWater(temp)/DENSITY_ICE-1.)*(pIce-STANDARD_AIR_PRESSURE));
        }


        /**
           \brief compute the temperature at which water with a given
           capillary pressure at a given gas pressure starts freezing.

           \param pC capillary pressure
           \param pG gas pressure
           \return freezing pressure

           The function is given in Spaans, E. J. A. and J. M. Baker:
           "The Soil Freezing Characteristic: Its Measurement and
           Similarity to the Soil Moisture Characteristic", SSSAJ 60:
           13-19 (1996) with a modification for the ice pressure by
           the author of this program.
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::FreezingTemperature(FLOAT pC,FLOAT pG)
        {
            FLOAT temp=273.15;
            FLOAT dx;
            FLOAT density;
            do
                {
                    density=PhysicalChemistry<FLOAT>::DensityWater(temp);
                    dx=(712380.*log(temp/273.15)-5545.*(temp-273.15)+3.14*(temp*temp-273.15*273.15)); /* [J/kg] */
                    dx*=density; /* [J/m3 = Pa] */
                    dx=pC-(dx-(density/DENSITY_ICE-1.)*(pG-STANDARD_AIR_PRESSURE));
                    dx/=density*((712380./temp)-5545.+6.28*temp*temp);
                    temp-=dx;
                }while(fabs(dx)>1e-13);
            return temp;
        }

        /**
           \brief compute the temperature dependend latent heat of the
           ice/water phase change

           The function is given in Spaans, E. J. A. and J. M. Baker:
           "The Soil Freezing Characteristic: Its Measurement and
           Similarity to the Soil Moisture Characteristic", SSSAJ 60:
           13-19 (1996) and the values are multiplied by the molar
           mass of water (0.018015 kg/mol)

           \param temp temperature
           \return latent heat
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::LatentHeatIceWater(FLOAT temp)
        {
            return ((-12834.+(99.89-0.113*temp)*temp)); /* [J/mol]  */
        }


        /**
           \brief compute the molar fraction of air dissolved in water
           at a given temperature and partial pressure of air

           The values are taken from the Handbook of chemistry and
           Physics, 76th Edition, CRC Press, 1995, page 6-5,
           6-6. Solubility is computed for 21% Oxygen and 79%
           Nitrogen.

           \param temp temperature
           \param pA   partial pressure of air
           \return molar fraction
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::SolubilityAir(FLOAT temp,FLOAT pA)
        {
            if (pA<std::numeric_limits<FLOAT>::epsilon())
                return(0.);

            FLOAT dum=temp/100.;
            FLOAT Xoxygen=std::exp(-66.7354+87.4755/dum+24.4526*log(dum));
            FLOAT Xnitrogen=std::exp(-67.3877+86.3213/dum+24.7981*log(dum));

            return (0.21*Xoxygen+0.79*Xnitrogen)*(pA/101325.);
        }


        /**
           This function computes tortuosity according to a function
           given by Millington and Quirk.

           \param saturation phase saturation
           \param porosity porosity
           \return tortuosity
        */

        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::TortuosityMQ(FLOAT saturation,FLOAT porosity)
        {
            return saturation*saturation*porosity*porosity/std::pow(porosity,2./3.);
        }


        /**
           This function computes the square root of the temperature
           dependend dielectricity constant of water

           \param temp temperature
           \return tortuosity
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::CompSqrtDielConstW(FLOAT temp)
        {
            FLOAT result=(-2.8e-8)*(temp-298.15);
            result=(result+1.19e-5)*(temp-298.15);
            result=(result-4.579e-3)*(temp-298.15);
            result=(result+1.0)*78.54;
            return sqrt(result);
        }

        /**
           computes enthalpy of pure water,
           confer IAPWS: "Revised Release on the IAPWS Industrial Formulation
           1997 for the Thermodynamic Properties of Water and Steam",
           http://www.iapws.org/relguide/IF97-Rev.pdf
           valid for temp < 623.15 K and press < 100 MPa

           \param temp temperature [K]
           \param press pressure [Pa]
           \return enthalpy [J/kg]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::EnthalpyWater(FLOAT temp, FLOAT press)
        {
            FLOAT spec_R = 461.526; // specific gas constant of water [(J/kg*K)]
            FLOAT tau = 1386/temp; // inverse reduced temperature [1/K]
            FLOAT pi = press/16530000; // reduced pressure [Pa]
            // partial derivative of Gibbs free energy to normalized temperature [dimensionless]
            FLOAT gamma_pi = 0;
            const int I[34] = {
                0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32};
            const int J[34] = {
                -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17,
                -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41};
            const FLOAT n[34] = {
                0.14632971213167, -0.84548187169114, -0.37563603672040e1,
                0.33855169168385e1, -0.95791963387872, 0.15772038513228,
                -0.16616417199501e-1, 0.81214629983568e-3, 0.28319080123804e-3,
                -0.60706301565874e-3, -0.18990068218419e-1, -0.32529748770505e-1,
                -0.21841717175414e-1, -0.52838357969930e-4, -0.47184321073267e-3,
                -0.30001780793026e-3, 0.47661393906987e-4, -0.44141845330846e-5,
                -0.72694996297594e-15,-0.31679644845054e-4, -0.28270797985312e-5,
                -0.85205128120103e-9, -0.22425281908000e-5, -0.65171222895601e-6,
                -0.14341729937924e-12,-0.40516996860117e-6, -0.12734301741641e-8,
                -0.17424871230634e-9, -0.68762131295531e-18, 0.14478307828521e-19,
                0.26335781662795e-22,-0.11947622640071e-22, 0.18228094581404e-23,
                -0.93537087292458e-25};

            for (int i=0; i<34; i++){
                gamma_pi += n[i] * pow(7.1-pi,I[i]) * pow(tau-1.222,J[i]-1) * J[i];
            }

            return tau*gamma_pi*spec_R*temp;
        }

        /**
           computes mass density of a brine-CO2 solution,
           confer Batzle & Wang (1992) (salt contribution)
           confer "Density of Aqueous Solutions of CO2" from Garcia (A (19)) (co2 contribution)

           \param t temperature [K]
           \param pl liquid phase pressure [Pa]
           \param sal salinity [mol/kg]
           \param x_l_g mass fraction of CO2 in water [kg/kg]
           \return density [kg/m³]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::DensityBrine(FLOAT t, FLOAT pl, FLOAT sal, FLOAT x_l_g)
        {
            FLOAT tc = t - 273.15; // t in °C
            FLOAT pmpa = pl/1000000; // pl in MPa
            // salinity as mass fraction [dimensionless]
            FLOAT x_l_b = MOLAR_MASS_NACL * sal;
            //mass fraction to mole fraction
            FLOAT c_l_g = ((x_l_g/MOLAR_MASS_CO2) /
                           (x_l_g/MOLAR_MASS_CO2 + (1-x_l_g)/MOLAR_MASS_H2O));
            FLOAT c_l_l = 1 - c_l_g;
            FLOAT rho_water = PhysicalChemistry<FLOAT>::DensityWater(t);

            FLOAT tmp = 80.0 - 3*tc - 3300*x_l_b - 13*pmpa + 47*pmpa*x_l_b;
            FLOAT brinedens = rho_water +
                1000*x_l_b*(0.668 + 0.44*x_l_b + 1.0E-6*(300*pmpa -2400*pmpa*x_l_b + tc*tmp));

            FLOAT M_T = MOLAR_MASS_H2O*c_l_l + MOLAR_MASS_CO2*c_l_g;
            FLOAT V_phi = (37.51-9.585e-2*tc+8.74e-4*pow(tc,2)-
                           5.044e-7*pow(tc,3))/1.0e6;
            FLOAT rho_inv = (c_l_g*V_phi/M_T
                             + MOLAR_MASS_H2O*c_l_l/(rho_water*M_T));
            FLOAT co2_contribution = 1/rho_inv - rho_water;
            return brinedens + co2_contribution; // [kg/m^3]
        }

        /**
           equation for reduced volume of CO2 at a given
           temperature and pressure.
           The Equation of State by Duan is used.

           see "An equation of state for the CH4-CO2-H20 system:
           I. Pure systems from 0 to 1000°C and 0 to 8000 bar",
           Geochimica ey Cosmochimica Acta Vol.56: 2605-2617

        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::ReducedVolumeCO2(FLOAT x, FLOAT temp, FLOAT press)
        {
            FLOAT a1 =  0.0899288497;
            FLOAT a2 = -0.494783127;
            FLOAT a3 =  0.0477922245;
            FLOAT a4 =  0.0103808883;
            FLOAT a5 = -0.0282516861;
            FLOAT a6 =  0.0949887563;
            FLOAT a7 =  0.00052060088;
            FLOAT a8 = -0.000293540971;
            FLOAT a9 = -0.00177265112;
            FLOAT a10 = -0.0000251101973;
            FLOAT a11 =  0.0000893353441;
            FLOAT a12 =  0.0000788998563;
            FLOAT a13 = -0.0166727022;
            FLOAT a14 =  1.398;
            FLOAT a15 =  0.0296;

            FLOAT TR = temp/CRIT_TEMP_CO2;
            FLOAT PR = press/CRIT_PRESS_CO2;
            FLOAT A = a1 + a2/(TR*TR) + a3/(TR*TR*TR);
            FLOAT B = a4 + a5/(TR*TR) + a6/(TR*TR*TR);
            FLOAT C = a7 + a8/(TR*TR) + a9/(TR*TR*TR);
            FLOAT D = a10 + a11/(TR*TR) + a12/(TR*TR*TR);
            return 1 - PR*x/TR + A/x + B/(x*x) + C/(x*x*x*x) + D/(x*x*x*x*x)
                + (a13/(TR*TR*TR*x*x))*(a14 + a15/(x*x))*exp(-a15/(x*x));
        }

        /**
           \brief compute mass density of CO2 at a given
           temperature and pressure.

           To compute the density the Equation of State by Duan is used:
           "An equation of state for the CH4-CO2-H20 system:
           I. Pure systems from 0 to 1000°C and 0 to 8000 bar",
           Geochimica ey Cosmochimica Acta Vol.56: 2605-2617

           \param temp temperature [K]
           \param press pressure [Pa]
           \return density [kg/m³]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::DensityCO2(FLOAT temp, FLOAT press)
        {
            if (temp < 273.15)
                temp = 273.15;
            if (press < 0.)
                press = 0.;

            FLOAT result;
            FLOAT xLo=0.1;
            FLOAT xHi=2000;
            FLOAT mid=0;
            unsigned int iter=0;
            do{
                ++iter;
                mid=0.5*(xLo+xHi);
                FLOAT fMid=PhysicalChemistry<FLOAT>::ReducedVolumeCO2(mid,temp,press);
                FLOAT fLo=PhysicalChemistry<FLOAT>::ReducedVolumeCO2(xLo,temp,press);
                FLOAT fHi=PhysicalChemistry<FLOAT>::ReducedVolumeCO2(xHi,temp,press);
                if (fMid*fLo<0)
                    xHi=mid;
                else if (fMid*fHi<0)
                    xLo=mid;
                else if (abs(fMid)<1e-6)
                    xLo=xHi=mid;
                else{
                    std::cout << "dens co2 at pascal pg=" << press << std::endl;
                    DUNE_THROW(MathError, "method of nested intervals did not converge! (Density CO2)");}
            } while ((xHi-xLo)>1e-10 && iter < 100);
            mid = (0.5*(xLo+xHi));
            result = mid*IDEAL_GAS_CONSTANT*CRIT_TEMP_CO2/CRIT_PRESS_CO2;  /* molar volume [m³/mol] */
            return (MOLAR_MASS_CO2/result);  /* [kg/m³] */
        }

        /**
           \brief compute viscosity of CO2 at a given
           temperature and pressure.

           see "The Viscosity of Carbon Dioxide" by Fenghour and Wakeham and Vesovic,
           J. Chem. Phys. 27: 31-44.
           The "viscosity in the critical region"-term is currently not implemented
           (has very little influence on result).

           \param temp temperature [K]
           \param press pressure [Pa]
           \return viscosity [kg/m/s]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::ViscosityCO2(FLOAT temp, FLOAT press)
        {
            FLOAT b0 =  0.235156;
            FLOAT b1 = -0.491266;
            FLOAT b2 =  5.211155e-2;
            FLOAT b3 =  5.347906e-2;
            FLOAT b4 = -1.537102e-2;
            FLOAT d11 =  0.4071119e-2;
            FLOAT d21 =  0.7198037e-4;
            FLOAT d64 =  0.2411697e-16;
            FLOAT d81 =  0.2971072e-22;
            FLOAT d82 = -0.1627888e-22;

            FLOAT dens = PhysicalChemistry<FLOAT>::DensityCO2(temp,press); /* [kg/m³] */
            FLOAT rT = temp/251.196;
            FLOAT viscosity_zero_density =
                1.00697*pow(temp,0.5) / exp(b0 + b1*log(rT) + b2*log(rT)*log(rT)
                                            + b3*log(rT)*log(rT)*log(rT)
                                            + b4*log(rT)*log(rT)*log(rT)*log(rT));
            FLOAT excess_viscosity =
                d11*dens + d21*dens*dens + d64*pow(dens,6)/pow(rT,3)
                + d81*pow(dens,8) + d82*pow(dens,8)/rT;
            return (viscosity_zero_density + excess_viscosity)/1000000; /* [kg/m/s] dynamic viscosity */
        }

        /**

           The Equation of State by Redlich-Kwong is used.
           see "CO2-H2O mixtures in the geological sequestration of CO2.
           I. Assessment and calculation of mutual solubilities
           from 12 to 100°C and up to 600 bar"
           by Spycher and Pruess and Ennis King,
           Geochimica et Cosmochimica Acta Vol.67: 3015-3031

           \param x x
           \param temp temperature [K]
           \param press pressure [Pa]
           \return molar volume [cm³/mol]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::MolarVolumeCO2(FLOAT x, FLOAT temp, FLOAT press)
        {
            FLOAT a_co2 = 7.54e7 - 4.13e4*temp;
            FLOAT b_co2 = 27.8;

            return x*x*x - x*x*10*IDEAL_GAS_CONSTANT*temp/press
                - x*10*IDEAL_GAS_CONSTANT*temp*b_co2/press
                + x*a_co2/(press*pow(temp,0.5))
                - x*b_co2*b_co2 - a_co2*b_co2/(press*pow(temp,0.5));
        }

        /**
           \brief solubility of water/co2 in co2/brine phase at a given
           temperature, pressure and salinity.

           see "CO2-H2O mixtures in the geological sequestration of CO2.
           I. Assessment and calculation of mutual solubilities
           from 12 to 100°C and up to 600 bar"
           by Spycher and Pruess and Ennis King,
           Geochimica et Cosmochimica Acta Vol.67: 3015-3031
           and
           "CO2-H2O Mixtures in the Geological Sequestration of CO2.
           II. Partitioning in Chloride Brines at 12-100°C and up to 600 bar"
           by Spycher and Pruess

           \param temp temperature [K]
           \param press pressure [Pa]
           \param salinity salinity [mol/kg]
           \param massfracwater massfraction water in co2 [kg/kg]
           \param massfracco2 massfraction co2 in brine [kg/kg]

        */
        template<typename FLOAT>
        void PhysicalChemistry<FLOAT>::BrineCO2Solubility(FLOAT temp, FLOAT press_pa,
                                                          FLOAT salinity,
                                                          FLOAT& massfracwater,
                                                          FLOAT& massfracco2)
        { //negative values for pressure < 7500 Pa!
            FLOAT press = press_pa/100000; //use bar for the computations
            FLOAT deltaP = press-1; // delta P = P - P_0, P_0 = 1 bar
            FLOAT v_h2o = 18.1; // average partial molar volume of H2O [cm^3/mol]
            FLOAT v_co2 = 32.6; // average partial molar volume of CO2 [cm^3/mol]
            FLOAT Tc = temp-273.15; // temparature in degrees
            FLOAT w = 1000/(MOLAR_MASS_H2O*1000); //1kg water are w moles
            FLOAT R = 10*IDEAL_GAS_CONSTANT; //gas constant [cm^3*bar/(K mol)]

            //intermolecular attraction and repulsion coefficients
            FLOAT a_co2 = 7.54e7 - 4.13e4*temp;
            FLOAT b_co2 = 27.8;
            FLOAT b_h2o = 18.18;
            FLOAT a_h2oco2 = 7.89e7;

            // molar volume of CO2
            FLOAT xLo=0;
            FLOAT xHi=1000000;
            FLOAT mid;
            unsigned int iter=0;
            do{
                ++iter;
                mid=0.5*(xLo+xHi);
                FLOAT fMid=PhysicalChemistry<FLOAT>::MolarVolumeCO2(mid,temp,press);
                FLOAT fLo=PhysicalChemistry<FLOAT>::MolarVolumeCO2(xLo,temp,press);
                FLOAT fHi=PhysicalChemistry<FLOAT>::MolarVolumeCO2(xHi,temp,press);
                if (fMid*fLo<0)
                    xHi=mid;
                else if (fMid*fHi<0)
                    xLo=mid;
                else if (abs(fMid)<1e-6)
                    xLo=xHi=mid;
                else
                    DUNE_THROW(MathError, "method of nested intervals did not converge!");
            } while ((xHi-xLo)>1e-10 && iter < 10000);
            mid = (0.5*(xLo+xHi));
            FLOAT v = mid; // molar volume of CO2

            //compute molar volume for the case of liquid CO2 (v<=96) (see figure B1)
            xLo=0;
            xHi=97;
            iter = 0;
            do{
                ++iter;
                mid=0.5*(xLo+xHi);
                FLOAT fMid=PhysicalChemistry<FLOAT>::MolarVolumeCO2(mid,temp,press);
                FLOAT fLo=PhysicalChemistry<FLOAT>::MolarVolumeCO2(xLo,temp,press);
                FLOAT fHi=PhysicalChemistry<FLOAT>::MolarVolumeCO2(xHi,temp,press);
                if (fMid*fLo<0)
                    xHi=mid;
                else if (fMid*fHi<0)
                    xLo=mid;
                else if (abs(fMid)<1e-6)
                    xLo=xHi=mid;
                else
                    xLo=xHi=-1; //molar volume is too big, v is already ok!
            } while ((xHi-xLo)>1e-10 && iter < 100);
            FLOAT v_l = (0.5*(xLo+xHi));

            if (v_l>0 && v_l<v)
                {
                    //determine which volume should be taken (liquid or gas)
                    FLOAT w1 = press*(v-v_l);
                    FLOAT w2 = R*temp*log((v-b_co2)/(v_l-b_co2)) +
                        a_co2*log((v+b_co2)*v_l/(v*(v_l+b_co2)))/
                        (pow(temp,0.5)*b_co2);
                    if ((w2-w1) < 0)
                        v=v_l;
                }

            // equilibrium constant H2O
            FLOAT k0_h2o =
                pow(10,-2.209 +  3.097e-2*Tc - 1.098e-4*Tc*Tc + 2.048e-7*Tc*Tc*Tc);
            // fugacity coefficient of H2O in waterphase
            FLOAT phi_h2o =
                exp(log(v/(v-b_co2)) + b_h2o/(v-b_co2)
                    - 2*a_h2oco2/(R*pow(temp,1.5)*b_co2) * log((v+b_co2)/v)
                    + (a_co2*b_h2o/(R*pow(temp,1.5)*b_co2*b_co2))
                    * (log((v+b_co2)/v) - b_co2/(v+b_co2)) - log(press*v/(R*temp)));
            // equilibrium constant CO2
            FLOAT k0_co2;
            if(Tc<31 && v<94) // different equilibrium constant for liquid co2
                k0_co2 = pow(10,1.169 +  1.368e-2*Tc - 5.380e-5*Tc*Tc);
            else
                k0_co2 = pow(10,1.189 +  1.304e-2*Tc - 5.446e-5*Tc*Tc);
            // fugacity coefficient of CO2
            FLOAT phi_co2 =
                exp(log(v/(v-b_co2)) + b_co2/(v-b_co2)
                    - (2*a_co2/(R*pow(temp,1.5)*b_co2)) * log((v+b_co2)/v)
                    + (a_co2*b_co2/(R*pow(temp,1.5)*b_co2*b_co2))
                    * (log((v+b_co2)/v) - b_co2/(v+b_co2)) - log(press*v/(R*temp)));
            //activity coefficient
            FLOAT nu =
                0.000336389723 - 0.0000198298980*temp + 0.00212220830*press/temp
                - 0.00524873303*press/(630 - temp);
            FLOAT lambda =
                - 0.411370585 + 0.000607632013*temp + 97.5347708/temp
                - 0.0237622469*press/temp + 0.0170656236*press/(630 - temp)
                + 0.0000141335834*temp*log(press);
            FLOAT act_coeff =
                exp(2*lambda*salinity + nu*salinity*salinity);

            FLOAT A = (k0_h2o / (phi_h2o*press)) * exp( deltaP*v_h2o / (R*temp));
            FLOAT B = (phi_co2*press / (w*k0_co2))
                * exp( -deltaP*v_co2 / (R*temp));

            FLOAT molh2o = (1.-B) / (1./A - B);
            FLOAT molco2 = B*(1.-molh2o);
            FLOAT mco2 = molco2*w / (1-molco2); //molality of aqueous CO2
            if (salinity > 0){
                mco2 = mco2/act_coeff;
                FLOAT molsalt = 2*salinity / (w+2*salinity+mco2); // mole fraction of NaCl
                molco2 = mco2 / (molsalt+w+mco2);
                molh2o = (1.-B-molsalt) / (1./A - B);
            }

            //use g/mol as unit for molar mass instead of kg/mol
            massfracco2 =
                (molco2*MOLAR_MASS_CO2*1000)/
                (molco2*MOLAR_MASS_CO2*1000 + (1-molco2)*MOLAR_MASS_H2O*1000); /* [kg/kg] */
            massfracwater =
                (molh2o*MOLAR_MASS_H2O*1000)/
                (molh2o*MOLAR_MASS_H2O*1000 + (1-molh2o)*MOLAR_MASS_CO2*1000);  /* [kg/kg] */
        }

        /**
           \brief equation for solubility of CO2 in water at a given
           temperature, pressure and salinity.

           see BrineCO2Solubility

           \param temp temperature [K]
           \param press pressure [Pa]
           \param salinity salinity [mol/kg]
           \return massfracco2 massfraction co2 in brine [kg/kg]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::SolubilityCO2InWater(FLOAT temp, FLOAT press,
                                                             FLOAT salinity)
        {
            if (press > 60000000) //only valid for p<=600 bar
                press = 60000000;
            FLOAT massfracco2inbrine = 0;
            FLOAT massfracwaterinco2 = 0;
            PhysicalChemistry::BrineCO2Solubility(temp,press,salinity,
                                                  massfracwaterinco2,massfracco2inbrine);
            if (massfracco2inbrine > 0)
                return massfracco2inbrine;
            else
                return 0.;
        }

        /**
           \brief equation for solubility of water in CO2 at a given
           temperature, pressure and salinity.

           see BrineCO2Solubility

           \param temp temperature [K]
           \param press pressure [Pa]
           \param salinity salinity [mol/kg]
           \return massfracwater massfraction water in co2 [kg/kg]
        */
        template<typename FLOAT>
        FLOAT PhysicalChemistry<FLOAT>::SolubilityWaterInCO2(FLOAT temp, FLOAT press,
                                                             FLOAT salinity)
        {
            if (press==0)
                return 0.;
            if (press > 60000000) //only valid for p<=600 bar
                press = 60000000;

            FLOAT massfracco2inbrine = 0;
            FLOAT massfracwaterinco2 = 0;
            PhysicalChemistry::BrineCO2Solubility(temp,press,salinity,
                                                  massfracwaterinco2,
                                                  massfracco2inbrine);
            if (massfracwaterinco2 > 0)
                return massfracwaterinco2;
            else
                return 0.;
        }

        template<typename FLOAT>
        class InterpolateCO2
        {
            typedef FLOAT FT;
        private:
            FT temp;
            FT sal;
            LinearInterpolationClass<FT> XlgTab_;
            LinearInterpolationClass<FT> XgwTab_;
            LinearInterpolationClass<FT> xlgTab_;
            LinearInterpolationClass<FT> xgwTab_;
            LinearInterpolationClass<FT> DensityTab_;
            LinearInterpolationClass<FT> ViscosityTab_;

            void ReadInterpolationValuesXlg(std::string filename)
            {
                XlgTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesXgw(std::string filename)
            {
                XgwTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesxlg(std::string filename)
            {
                xlgTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesxgw(std::string filename)
            {
                xgwTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesDensity(std::string filename)
            {
                DensityTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationValuesViscosity(std::string filename)
            {
                ViscosityTab_.ReadInterpolationValues(filename);
            }
            void ReadInterpolationTableXlg(std::string filename)
            {
                XlgTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableXgw(std::string filename)
            {
                XgwTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTablexlg(std::string filename)
            {
                xlgTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTablexgw(std::string filename)
            {
                xgwTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableDensity(std::string filename)
            {
                DensityTab_.ReadInterpolationTable(filename);
            }
            void ReadInterpolationTableViscosity(std::string filename)
            {
                ViscosityTab_.ReadInterpolationTable(filename);
            }
            void WriteInterpolationTableXlg(std::string filename)
            {
                XlgTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableXgw(std::string filename)
            {
                XgwTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTablexlg(std::string filename)
            {
                xlgTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTablexgw(std::string filename)
            {
                xgwTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableDensity(std::string filename)
            {
                DensityTab_.WriteInterpolationTable(filename);
            }
            void WriteInterpolationTableViscosity(std::string filename)
            {
                ViscosityTab_.WriteInterpolationTable(filename);
            }
            void InterpolateXlg(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=0.;
                FT maxX=200000000;
                FT interval=(maxX-minX)/(numPoints-1.);
                y[0]=0.;
                for (int i=1;i<numPoints;i++)
                    y[i]=PhysicalChemistry<FT>::
                        SolubilityCO2InWater(temp, minX+i*interval,sal);
                XlgTab_.Init(minX,maxX,y);
            }
            void InterpolateXgw(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=100000;
                FT maxX=200000000;
                FT interval=(maxX-minX)/(numPoints-1.);
                y[0]=PhysicalChemistry<FT>::
                    SolubilityWaterInCO2(temp, minX,sal);
                for (int i=1;i<numPoints;i++)
                    y[i]=PhysicalChemistry<FT>::
                        SolubilityWaterInCO2(temp, minX+i*interval,sal);
                XgwTab_.Init(minX,maxX,y);
            }
            void Interpolatexlg(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=0.;
                FT maxX=200000000;
                FT interval=(maxX-minX)/(numPoints-1.);
                y[0]=0.;
                for (int i=1;i<numPoints;i++){
                    FT tmp = PhysicalChemistry<FT>::
                        SolubilityCO2InWater(temp, minX+i*interval,sal);
                    //mass fraction to mole fraction
                    y[i]=(tmp*PhysicalChemistry<FT>::MOLAR_MASS_H2O)/
                        (tmp*PhysicalChemistry<FT>::MOLAR_MASS_H2O + (1-tmp)*PhysicalChemistry<FT>::MOLAR_MASS_CO2);}
                xlgTab_.Init(minX,maxX,y);
            }
            void Interpolatexgw(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=100000;
                FT maxX=200000000;
                FT interval=(maxX-minX)/(numPoints-1.);
                FT tmp = PhysicalChemistry<FT>::
                    SolubilityWaterInCO2(temp, minX,sal);
                y[0]=(tmp*PhysicalChemistry<FT>::MOLAR_MASS_CO2)/
                    (tmp*PhysicalChemistry<FT>::MOLAR_MASS_CO2 + (1-tmp)*PhysicalChemistry<FT>::MOLAR_MASS_H2O);;
                for (int i=1;i<numPoints;i++){
                    tmp = PhysicalChemistry<FT>::
                        SolubilityWaterInCO2(temp, minX+i*interval,sal);
                    y[i]=(tmp*PhysicalChemistry<FT>::MOLAR_MASS_CO2)/
                        (tmp*PhysicalChemistry<FT>::MOLAR_MASS_CO2 + (1-tmp)*PhysicalChemistry<FT>::MOLAR_MASS_H2O);}
                xgwTab_.Init(minX,maxX,y);
            }
            void InterpolateDensity(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=0.;
                FT maxX=100000000; //max 1000 bar
                FT interval=(maxX-minX)/(numPoints-1.);
                y[0]=0.;
                for (int i=1;i<numPoints;i++)
                    y[i]=PhysicalChemistry<FT>::DensityCO2(temp,minX+i*interval);
                DensityTab_.Init(minX,maxX,y);
            }
            void InterpolateViscosity(int numPoints)
            {
                std::vector<FT> y(numPoints);
                FT minX=0.;
                FT maxX=100000000; //max 1000 bar
                FT interval=(maxX-minX)/(numPoints-1.);
                y[0]=PhysicalChemistry<FT>::ViscosityCO2(temp,minX);
                for (int i=1;i<numPoints;i++)
                    y[i]=PhysicalChemistry<FT>::ViscosityCO2(temp,minX+i*interval);
                ViscosityTab_.Init(minX,maxX,y);
            }

        public:
            InterpolateCO2(FT temp_,FT sal_) :
                temp(temp_),
                sal(sal_)
            {}

            void Init(int numPoints=10000)
            {
                InterpolateXlg(numPoints);
                InterpolateXgw(numPoints);
                Interpolatexlg(numPoints);
                Interpolatexgw(numPoints);
                InterpolateDensity(numPoints);
                InterpolateViscosity(numPoints);
            }

            FT Xlg(FT pg)
            {
                if (pg > 60000000)
                    return(XlgTab_.GetValue(60000000));
                else
                    return(XlgTab_.GetValue(pg));
            }

            FT Xgw(FT pg)
            {
                if (pg > 60000000)
                    return(XgwTab_.GetValue(60000000));
                else if ( pg < 100000 )
                    return (XgwTab_.GetValue(100000));
                else
                    return(XgwTab_.GetValue(pg));
            }

            FT xlg(FT pg)
            {
                if (pg > 60000000)
                    return(xlgTab_.GetValue(60000000));
                else
                    return(xlgTab_.GetValue(pg));
            }

            FT xgw(FT pg)
            {
                if (pg > 60000000)
                    return(xgwTab_.GetValue(60000000));
                else if ( pg < 100000 )
                    return (xgwTab_.GetValue(100000));
                else
                    return(xgwTab_.GetValue(pg));
            }

            FT density(FT pg)
            {
                if (pg > 100000000)
                    return(DensityTab_.GetValue(100000000));
                else
                    return(DensityTab_.GetValue(pg));
            }

            FT viscosity(FT pg)
            {
                if (pg > 100000000)
                    return(ViscosityTab_.GetValue(100000000));
                else
                    return(ViscosityTab_.GetValue(pg));
            }

        };

    } // end namespace PM
} // end namespace Dune
