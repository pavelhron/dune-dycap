// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set ts=4 sw=2 et sts=4:
#ifndef PHYSICAL_CHEMISTRY_H
#define PHYSICAL_CHEMISTRY_H

namespace Dune {
    namespace PM {

        /**
           \brief handle physical and chemical constants and relations
        */
        template<typename FLOAT>
        class PhysicalChemistry
        {
        public:
            static FLOAT ViscosityWater(FLOAT temp);
            static FLOAT ViscosityAir(FLOAT temp);
            static FLOAT RelSurfaceTensionWaterAir(FLOAT temp);
            static FLOAT DensityWater(FLOAT temp);
            static FLOAT MolarDensityWater(FLOAT temp);
            static FLOAT FreezingPressure(FLOAT temp,FLOAT pIce);
            static FLOAT FreezingTemperature(FLOAT pC,FLOAT pG);
            static FLOAT LatentHeatIceWater(FLOAT temp);
            static FLOAT MolarDensityVapor (FLOAT pC, FLOAT temp, FLOAT XWL);
            static FLOAT MolarDensityVaporVT (FLOAT pL, FLOAT pA, FLOAT temp, FLOAT XWL);
            static FLOAT TortuosityMQ(FLOAT saturation,FLOAT porosity);
            static FLOAT SolubilityAir(FLOAT temp,FLOAT pA);
            static FLOAT CompSqrtDielConstW(FLOAT temp);
            static FLOAT EnthalpyWater(FLOAT temp, FLOAT press);
            static FLOAT DensityBrine(FLOAT temp, FLOAT pl, FLOAT sal, FLOAT x_l_g);
            static FLOAT DensityCO2(FLOAT temp, FLOAT press);
            static FLOAT ViscosityCO2(FLOAT temp, FLOAT press);
            static FLOAT SolubilityCO2InWater(FLOAT temp, FLOAT press, FLOAT salinity);
            static FLOAT SolubilityWaterInCO2(FLOAT temp, FLOAT press, FLOAT salinity);

            static const FLOAT IDEAL_GAS_CONSTANT;     /* [J/mol/K] Physikalische Blätter 3/2000       */
            static const FLOAT MOLAR_MASS_WATER;
            static const FLOAT MOLAR_VOLUME_WATER;
            static const FLOAT MOLAR_DENSITY_ICE;
            static const FLOAT MOLAR_MASS_AIR;
            static const FLOAT MOLAR_MASS_CO2;
            static const FLOAT MOLAR_MASS_H2O;
            static const FLOAT DENSITY_ICE;
            static const FLOAT DENSITY_WATER;
            static const FLOAT DIFF_COEFF_H2O_WATER;
            static const FLOAT DIFF_COEFF_AIR_WATER;
            static const FLOAT DIFF_COEFF_H2O_AIR;
            static const FLOAT MOLAR_LATHEAT_WATVAP;
            static const FLOAT GRAVITY;
            static const FLOAT STANDARD_AIR_PRESSURE;
            static const FLOAT SQRT_DIELECTRIC_CONST_AIR;
            static const FLOAT SQRT_DIELECTRIC_CONST_MATRIX;
            static const FLOAT SQRT_DIELECTRIC_CONST_ICE;
            static const FLOAT MOLAR_HEAT_CAP_WATER;
            static const FLOAT MOLAR_HEAT_CAP_ICE;
            static const FLOAT MOLAR_HEAT_CAP_AIR;
            static const FLOAT MOLAR_HEAT_CAP_VAPOR;
            static const FLOAT CRIT_PRESS_CO2;
            static const FLOAT CRIT_TEMP_CO2;
            static const FLOAT MOLAR_MASS_NACL;
            static const FLOAT MOLAR_MASS_NA;
            static const FLOAT MOLAR_MASS_CL;
        private:
            static FLOAT ReducedVolumeCO2(FLOAT x, FLOAT temp, FLOAT press);
            static FLOAT MolarVolumeCO2(FLOAT x, FLOAT temp, FLOAT press);
            static void BrineCO2Solubility(FLOAT temp, FLOAT press, FLOAT salinity,
                                           FLOAT& massfracwaterinco2,
                                           FLOAT& massfracco2inbrine);
            static const FLOAT VAPOR_LIMIT;
        };

        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::IDEAL_GAS_CONSTANT = 8.314472;     /* [J/mol/K] Physikalische Blätter 3/2000       */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_WATER = 0.018015;       /* [kg/mol] Kuchling: Handbuch der Physik       */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_VOLUME_WATER = 1.8047e-05;   /* [m3/mol] 20°C Kuchling: Handbuch der Physik  */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_DENSITY_ICE = 50902.0;       /* [mol/m3] 0°C Kuchling: Handbuch der Physik   */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_AIR = 0.026;            /* [kg/mol] 21% oxygen, 79% nitrogen data from Kuchling: Handbuch der Physik    */
        /* => according to Wikipedia: 0.0289644 kg/mol (http://de.wikipedia.org/wiki/Universelle_Gaskonstante) corresponds to specific gas constant */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_CO2 = 0.044;            /* [kg/mol] 27.3% carbon 72.7% oxygen http://www.webqc.org/molecular-weight-of-CO2.html    */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_H2O = 0.018;            /* [kg/mol] 11.2% hydrogem 88.8% oxygen http://www.webqc.org/molecular-weight-of-H2O.html    */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::DENSITY_ICE = 917.0;               /* [kg/m3] 0°C Kuchling: Handbuch der Physik    */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::DENSITY_WATER = 998.227;           /* [kg/m3] 0°C Kuchling: Handbuch der Physik    */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::DIFF_COEFF_H2O_WATER = 2.26e-9;    /* [m2/s] Atkins: Physikalische Chemie          */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::DIFF_COEFF_AIR_WATER = 2e-9;       /* [m2/s] estimated                             */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::DIFF_COEFF_H2O_AIR = 0.242e-4;     /* [m2/s] Handbook of Chemistry and Physics     */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_LATHEAT_WATVAP = 40656.;     /* [J/mol] 100°C Atikins: Physikalische Chemie  */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::GRAVITY = 9.81;                    /* [m/s2]   */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::STANDARD_AIR_PRESSURE = 101300.0;  /* [Pa]        */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::SQRT_DIELECTRIC_CONST_AIR = 1.0;
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::SQRT_DIELECTRIC_CONST_MATRIX = 1.96;
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::SQRT_DIELECTRIC_CONST_ICE = 1.96;
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_HEAT_CAP_WATER = 75.7528;    /* [J/mol/K]  =  4.19753e6 J/m3/K bei 20°C       */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_HEAT_CAP_ICE = 37.1;         /* [J/mol/K]  =  1.89e6 J/m3/K bei 0°C            */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_HEAT_CAP_AIR = 28.98;        /* [J/mol/K]  =  1005 J/kg/K bei 20°C Kuchling: Handbuch der Physik */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_HEAT_CAP_VAPOR = 28.98;      /* [J/mol/K]  =  1005 J/kg/K bei 20°C Kuchling: Handbuch der Physik  = Luft) */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::VAPOR_LIMIT = 1.0E-13;             /* convergence limit for vapor pressure computation */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::CRIT_PRESS_CO2 = 7382500;          /* [Pa] critical pressure for co2 from Duan: Equation of state for CH4, C02, and H20 */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::CRIT_TEMP_CO2 = 304.2;             /* [K] critical temperature for co2 from Duan: Equation of state for CH4, C02, and H20 */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_NACL = 0.058;        /* [kg/mol] 39.3% sodium 60.7% chlorine http://www.webqc.org/molecular-weight-of-NaCl.html  */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_NA =  0.023;        /* [kg/mol] http://www.webqc.org/molecular-weight-of-Na.html  */
        template<typename FLOAT>
        const FLOAT PhysicalChemistry<FLOAT>::MOLAR_MASS_CL = 0.035;           /* [kg/mol] http://www.webqc.org/molecular-weight-of-Cl.html  */

    } //end namespace PM
} // end namespace Dune

#include "physical_chemistry.cc"

#endif
