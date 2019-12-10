//  Author: Robert Scheller, Melissa Lucash
//  Modified: Chihiro Haga (chihiro.haga@ge.see.eng.osaka-u.ac.jp)

using System.Collections.Generic;
using System.IO;
using System;
using System.Globalization;
using Landis.Core;
using Landis.SpatialModeling;
using Landis.Library.Climate;

namespace Landis.Extension.Succession.NECN
{

    public enum WaterType { Linear, Ratio }


    public class SoilWater
    {
        private static double Precipitation;
        private static double H2Oinputs;
        private static double tave;
        private static double tmax;
        private static double tmin;
        private static double pet;
        private static int daysInMonth;
        private static int beginGrowing;
        private static int endGrowing;
        private static double holdingTank = 0.0;


        public static void Run(int year, int month, double liveBiomass, Site site, out double baseFlow, out double stormFlow, out double AET)
        {

            //Originally from h2olos.f of CENTURY model
            //Water Submodel for Century - written by Bill Parton
            //     Updated from Fortran 4 - rm 2/92
            //     Rewritten by Bill Pulliam - 9/94
            //     Rewritten by Melissa Lucash- 11/2014

            PlugIn.ModelCore.UI.WriteLine("month={0}.", Main.Month);
        
            //...Initialize Local Variables
            double addToSoil = 0.0;
            double bareSoilEvap = 0.0;
            baseFlow = 0.0;
            double relativeWaterContent = 0.0;
            double snow = 0.0;
            stormFlow = 0.0;
            double actualET = 0.0;
            double remainingPET = 0.0;
            double availableWaterMax = 0.0;  //amount of water available after precipitation and snowmelt (over-estimate of available water)
            //double availableWaterMin = 0.0;   //amount of water available after stormflow (runoff) evaporation and transpiration, but before baseflow/leaching (under-estimate of available water)
            double availableWater = 0.0;     //amount of water deemed available to the trees, which will be the average between the max and min
            double priorWaterAvail = SiteVars.AvailableWater[site];
            double waterFull = 0.0;

            //...Calculate external inputs
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];

            double litterBiomass = (SiteVars.SurfaceStructural[site].Carbon + SiteVars.SurfaceMetabolic[site].Carbon) * 2.0;
            double deadBiomass = SiteVars.SurfaceDeadWood[site].Carbon / 0.47;
            double soilWaterContent = SiteVars.SoilWaterContent[site];
            double liquidSnowpack = SiteVars.LiquidSnowPack[site];

            H2Oinputs = ClimateRegionData.AnnualWeather[ecoregion].MonthlyPrecip[month]; //rain + irract in cm;
            Precipitation = ClimateRegionData.AnnualWeather[ecoregion].MonthlyPrecip[month]; //rain + irract in cm;
            //PlugIn.ModelCore.UI.WriteLine("SoilWater. WaterInputs={0:0.00}, .", H2Oinputs);
            tave = ClimateRegionData.AnnualWeather[ecoregion].MonthlyTemp[month];
            //PlugIn.ModelCore.UI.WriteLine("SoilWater. AvgTemp={0:0.00}, .", tave);
            tmax = ClimateRegionData.AnnualWeather[ecoregion].MonthlyMaxTemp[month];
            tmin = ClimateRegionData.AnnualWeather[ecoregion].MonthlyMinTemp[month];
            pet = ClimateRegionData.AnnualWeather[ecoregion].MonthlyPET[month];
            daysInMonth = AnnualClimate.DaysInMonth(month, year);
            beginGrowing = ClimateRegionData.AnnualWeather[ecoregion].BeginGrowing;
            endGrowing = ClimateRegionData.AnnualWeather[ecoregion].EndGrowing;

            double wiltingPoint = SiteVars.SoilWiltingPoint[site];
            double soilDepth = SiteVars.SoilDepth[site]; 
            double fieldCapacity = SiteVars.SoilFieldCapacity[site];
            double stormFlowFraction = SiteVars.SoilStormFlowFraction[site]; 
            double baseFlowFraction = SiteVars.SoilBaseFlowFraction[site];
            double drain = SiteVars.SoilDrain[site]; 
           
                      
            //...Calculating snow pack first. Occurs when mean monthly air temperature is equal to or below freezing,
            //     precipitation is in the form of snow.
            
            if (tmin <= 0.0) // Use tmin to dictate whether it snows or rains. 
            {
                snow = H2Oinputs; 
                H2Oinputs = 0.0;  
                liquidSnowpack += snow;  //only tracking liquidsnowpack (water equivalent) and not the actual amount of snow on the ground (i.e. not snowpack).
                //PlugIn.ModelCore.UI.WriteLine("Let it snow!! snow={0}, liquidsnowpack={1}.", snow, liquidSnowpack);
            }
            else
            {
                //PH: Accumulate precipitation and snowmelt before adding to soil so that interception and soil evaporation can come out first
                addToSoil += H2Oinputs;
                //soilWaterContent += H2Oinputs;
                //PlugIn.ModelCore.UI.WriteLine("Let it rain and add it to soil! rain={0}, soilWaterContent={1}.", H2Oinputs, soilWaterContent);
            }

           
            //...Then melt snow if there is snow on the ground and air temperature (tmax) is above minimum.            
            if (liquidSnowpack > 0.0 && tmax > 0.0)
            {
                //...Calculate the amount of snow to melt:
                //This relationship ultimately derives from http://www.nps.gov/yose/planyourvisit/climate.htm which described the relationship between snow melting and air temp.
                //Documentation for the regression equation is in spreadsheet called WaterCalcs.xls by M. Lucash
                double snowMeltFraction = Math.Max((tmax * 0.05) + 0.024, 0.0);//This equation assumes a linear increase in the fraction of snow that melts as a function of air temp.  

               if (snowMeltFraction > 1.0)
                    snowMeltFraction = 1.0;

               //PH: 
                addToSoil += liquidSnowpack * snowMeltFraction;  //PH: Add melted snow to addToSoil. Amount of liquidsnowpack that melts = liquidsnowpack multiplied by the fraction that melts.
                //addToSoil = liquidSnowpack * snowMeltFraction;  //Amount of liquidsnowpack that melts = liquidsnowpack multiplied by the fraction that melts.
              
                //Subtracted melted snow from snowpack and add it to the soil
               liquidSnowpack = liquidSnowpack - addToSoil;  
               //PH: Add to soil later. 
                //soilWaterContent += addToSoil;
            }
            
            //Calculate the max amout of water available to trees, an over-estimate of the water available to trees.  It only reflects precip and melting of precip.
            //PH: available water may need to be calculated differently with my proposed changes, but I moved this variable down for now so that soilEvaporation comes out first.
            //availableWaterMax = soilWaterContent;
            
            //...Evaporate water from the snow pack (rewritten by Pulliam 9/94)
            //...Coefficient 0.87 relates to heat of fusion for ice vs. liquid water
            if (liquidSnowpack > 0.0)
            {
                //...Calculate cm of snow that remaining pet energy can evaporate:
                double evaporatedSnow = pet * 0.87;

                //...Don't evaporate more snow than actually exists:
                if (evaporatedSnow > liquidSnowpack)
                    evaporatedSnow = liquidSnowpack;

                liquidSnowpack = liquidSnowpack - evaporatedSnow;

                //...Decrement remaining pet by energy used to evaporate snow:

                //PH: CENTURY code divides evaporatedSnow by 0.87 so it matches the PET used to melt snow.
                remainingPET = pet - evaporatedSnow / 0.87;
                //remainingPET = pet - evaporatedSnow;
                
                if (remainingPET < 0.0) 
                    remainingPET = 0.0;

                //Subtract evaporated snowfrom the soil water content
                //PH: Take evaporated snow out of snowmelt instead of soil...or is that double counting evaporation? Could this cause addToSoil to go negative?  
                //It is already taken out of the snowpack, and is used to decrement PET so that it won't affect AET
                addToSoil -= evaporatedSnow;
                if (addToSoil < 0.0)
                    addToSoil = 0.0;
                //soilWaterContent -= evaporatedSnow;
            }
            
            // ********************************************************
            //...Calculate bare soil water loss and interception  when air temperature is above freezing and no snow cover.
            //...Mofified 9/94 to allow interception when t < 0 but no snow cover, Pulliam
            //PH: I moved this up to remove intercepted precipitation and bare soil evaporation from accumulated precipitation and snowmelt before it goes to the soil.
            if (liquidSnowpack <= 0.0)
            {
                //...Calculate total canopy cover and litter, put cap on effects:
                double standingBiomass = liveBiomass + deadBiomass;

                if (standingBiomass > 800.0) standingBiomass = 800.0;
                if (litterBiomass > 400.0) litterBiomass = 400.0;

                //...canopy interception, fraction of  precip (canopyIntercept):
                double canopyIntercept = ((0.0003 * litterBiomass) + (0.0006 * standingBiomass)) * OtherData.WaterLossFactor1;

                //...Bare soil evaporation, fraction of precip (bareSoilEvap):
                bareSoilEvap = 0.5 * System.Math.Exp((-0.002 * litterBiomass) - (0.004 * standingBiomass)) * OtherData.WaterLossFactor2;
                
                //...Calculate total surface evaporation losses, maximum allowable is 0.4 * pet. -rm 6/94
                remainingPET = pet;
                double soilEvaporation = System.Math.Min(((bareSoilEvap + canopyIntercept) * H2Oinputs), (0.4 * remainingPET));
                PlugIn.ModelCore.UI.WriteLine("SWdebug soilEvaporation={0:0.00}, .", soilEvaporation); //Chihiro:

                //Subtract soil evaporation from soil water content
                //PH: Subtract soilEvaporation from addToSoil so it won't drive down soil water. 
                //PH: SoilEvaporation represents water that evaporates before reaching soil, so should not be subtracted from soil.
                addToSoil -= soilEvaporation;
                //soilWaterContent -= soilEvaporation;
            }

            PlugIn.ModelCore.UI.WriteLine("SWdebug addToSoil={0:0.00}, .", addToSoil); //Chihiro:

            //PH: Add liquid water to soil
            soilWaterContent += addToSoil;
            //Calculate the max amout of water available to trees, an over-estimate of the water available to trees.  It only reflects precip and melting of precip.
            //PH: Moved down so that soilEvaporation comes out first.
            availableWaterMax = soilWaterContent;
            // ********************************************************


            // Calculate actual evapotranspiration.  This equation is derived from the stand equation for calculating AET from PET
            //  Bergström, 1992

            // ********************************************************
            //PH: Moved up to take evapotranspiration out before excess drains away. This is different from the CENTURY approach, where evaporation is taken out of the add first, but is 
            //less complex because it does not require partitioning the evaporation if evapotranspiration exceeds addToSoil.

            double waterEmpty = wiltingPoint * soilDepth;
            waterFull = soilDepth * fieldCapacity;  //units of cm


            if (soilWaterContent > waterFull)
                actualET = remainingPET;
            else
            {
                actualET = Math.Max(remainingPET * ((soilWaterContent - waterEmpty) / (waterFull - waterEmpty)), 0.0);
            }

            if (actualET < 0.0)
                actualET = 0.0;
            AET = actualET;
            PlugIn.ModelCore.UI.WriteLine("SWdebug actualET={0:0.00}, .", actualET); //Chihiro:
            // ********************************************************

            //PlugIn.ModelCore.UI.WriteLine("AET {0} = ", AET);

            //Subtract transpiration from soil water content
            soilWaterContent -= actualET;
            PlugIn.ModelCore.UI.WriteLine("SWdebug soilWaterContent={0:0.00}, .", soilWaterContent); //Chihiro:

            //Allow excess water to run off during storm events (stormflow)
            double waterMovement = 0.0;            

            if (soilWaterContent > waterFull)
            {

                waterMovement = Math.Max((soilWaterContent - waterFull), 0.0); // How much water should move during a storm event, which is based on how much water the soil can hold.
                soilWaterContent = waterFull;
                
                //...Compute storm flow.
                stormFlow = waterMovement * stormFlowFraction;

                // ********************************************************
                //Subtract stormflow from soil water
                //PH: I don't see why this should come out of soil water. It should come out of excess water
                //soilWaterContent -= stormFlow;
                //PlugIn.ModelCore.UI.WriteLine("Water Runs Off. stormflow={0}.", stormFlow);
                // ********************************************************
            }

            // ********************************************************
            //PH: add new variable to track excess water and calclulate baseFlow
            //PH: Remove stormFlow from from excess water
            waterMovement -= stormFlow;
            holdingTank += waterMovement;
            // ********************************************************


            // ********************************************************
            //PH: moved up to take out soil evaporation and interception before water is added to soil. 
            //...Calculate bare soil water loss and interception  when air temperature is above freezing and no snow cover.
            //...Mofified 9/94 to allow interception when t < 0 but no snow cover, Pulliam
            //if (liquidSnowpack <= 0.0)
            //{
            //...Calculate total canopy cover and litter, put cap on effects:
            //     double standingBiomass = liveBiomass + deadBiomass;
            //
            //     if (standingBiomass > 800.0) standingBiomass = 800.0;
            //     if (litterBiomass > 400.0) litterBiomass = 400.0;

            //...canopy interception, fraction of  precip (canopyIntercept):
            //double canopyIntercept = ((0.0003 * litterBiomass) + (0.0006 * standingBiomass)) * OtherData.WaterLossFactor1;

            //...Bare soil evaporation, fraction of precip (bareSoilEvap):
            //bareSoilEvap = 0.5 * System.Math.Exp((-0.002 * litterBiomass) - (0.004 * standingBiomass)) * OtherData.WaterLossFactor2;

            //...Calculate total surface evaporation losses, maximum allowable is 0.4 * pet. -rm 6/94
            //remainingPET = pet;
            //double soilEvaporation = System.Math.Min(((bareSoilEvap + canopyIntercept) * H2Oinputs), (0.4 * remainingPET));

            //Subtract soil evaporation from soil water content
            //soilWaterContent -= soilEvaporation;
            //}

            // Calculate actual evapotranspiration.  This equation is derived from the stand equation for calculating AET from PET
            //  Bergström, 1992

            //double waterEmpty = wiltingPoint * soilDepth;

            //if (soilWaterContent > waterFull)
            //    actualET = remainingPET;
            //else
            //{
            //    actualET = Math.Max(remainingPET * ((soilWaterContent - waterEmpty) / (waterFull - waterEmpty)), 0.0);
            //}

            //if (actualET < 0.0)
            //    actualET = 0.0;
            //AET = actualET;

            //PlugIn.ModelCore.UI.WriteLine("AET {0} = ", AET);

            //Subtract transpiration from soil water content
            //soilWaterContent -= actualET;
            // ********************************************************


            // ********************************************************
            //Leaching occurs. Drain baseflow fraction from holding tank.
            //PH: Now baseflow comes from holding tank.
            baseFlow = holdingTank * baseFlowFraction;
            //baseFlow = soilWaterContent * baseFlowFraction;

            //Subtract baseflow from soil water
            //PH: Subtract from holding tank instead. To not deplete soil water but still allow estimation of baseFlow.
            holdingTank -= baseFlow;
            //soilWaterContent -= baseFlow;
            // ********************************************************


            //Calculate the amount of available water after all the evapotranspiration and leaching has taken place (minimum available water)           
            availableWater = Math.Max(soilWaterContent - waterEmpty, 0.0);

            //Calculate the final amount of available water to the trees, which is the average of the max and min          
            //PH: availableWater is affected by my changes, and soilWaterContent should be higher now.  Therefore, I propose calculating using soilWaterContent directly instead
            //availableWater = soilWaterContent - waterEmpty;
            
            // Here calculate available water at the midpoint of the month
            //availableWater = (availableWaterMax + availableWaterMin)/ 2.0;

            // Compute the ratio of precipitation to PET
            double ratioPrecipPET = 0.0;
            if (pet > 0.0) ratioPrecipPET = availableWater / pet;  //assumes that the ratio is the amount of incoming precip divided by PET.

            SiteVars.AnnualPPT_AET[site] += actualET; // RMS:  Currently using this to test AET by itself // Precipitation - actualET;
            SiteVars.AnnualClimaticWaterDeficit[site] += (pet - actualET) * 10.0;  // Convert to mm, the standard definition
            //PlugIn.ModelCore.UI.WriteLine("Month={0}, PET={1}, AET={2}.", month, pet, actualET);

            SiteVars.LiquidSnowPack[site] = liquidSnowpack;
            SiteVars.WaterMovement[site] = waterMovement;
            SiteVars.AvailableWater[site] = availableWater;  //available to plants for growth     
            SiteVars.SoilWaterContent[site] = soilWaterContent;
            SiteVars.SoilTemperature[site] = CalculateSoilTemp(tmin, tmax, liveBiomass, litterBiomass, month);
            SiteVars.DecayFactor[site] = CalculateDecayFactor((int)OtherData.WType, SiteVars.SoilTemperature[site], relativeWaterContent, ratioPrecipPET, month);
            SiteVars.AnaerobicEffect[site] = CalculateAnaerobicEffect(drain, ratioPrecipPET, pet, tave, site, month);
            // chihiro; add monthly variables
            SiteVars.MonthlyWaterMovement[site][month]      = waterMovement;
            SiteVars.MonthlyBaseFlow[site][month]           = baseFlow;
            SiteVars.MonthlyStormFlow[site][month]          = stormFlow;
            SiteVars.MonthlyPet[site][month]                = pet;
            SiteVars.MonthlyLiquidSnowPack[site][month]     = liquidSnowpack;
            SiteVars.MonthlyAvailableWater[site][month]     = availableWater;
            SiteVars.MonthlySoilWaterContent[site][month]   = soilWaterContent;
            SiteVars.MonthlySoilTemperature[site][month]    = SiteVars.SoilTemperature[site];
            SiteVars.MonthlyDecayFactor[site][month]        = SiteVars.DecayFactor[site];
            SiteVars.MonthlyAnaerobicEffect[site][month]    = SiteVars.AnaerobicEffect[site];
            if (month == 0)
                SiteVars.DryDays[site] = 0;
            else
                SiteVars.DryDays[site] += CalculateDryDays(month, beginGrowing, endGrowing, waterEmpty, availableWater, priorWaterAvail);
                        
            return;
        }

        private static int CalculateDryDays(int month, int beginGrowing, int endGrowing, double wiltingPoint, double waterAvail, double priorWaterAvail)
        {
            //PlugIn.ModelCore.UI.WriteLine("Month={0}, begin={1}, end={2}.", month, beginGrowing, endGrowing);
            int[] julianMidMonth = { 15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 };
            int dryDays = 0;
            int julianDay = julianMidMonth[month]; 
            int oldJulianDay = julianMidMonth[month-1]; 
            double dryDayInterp = 0.0;
            //PlugIn.ModelCore.UI.WriteLine("Month={0}, begin={1}, end={2}, wiltPt={3:0.0}, waterAvail={4:0.0}, priorWater={5:0.0}.", month, beginGrowing, endGrowing, wiltingPoint, waterAvail, priorWaterAvail);
            
            //Increment number of dry days, truncate at end of growing season
                if ((julianDay > beginGrowing) && (oldJulianDay < endGrowing)) 
                {
                    if ((priorWaterAvail >= wiltingPoint)  && (waterAvail >= wiltingPoint))
                        {
                        dryDayInterp += 0.0;  // NONE below wilting point
                    }
                    else if ((priorWaterAvail > wiltingPoint) && (waterAvail < wiltingPoint)) 
                    {
                        dryDayInterp = daysInMonth * (wiltingPoint - waterAvail) / 
                                        (priorWaterAvail - waterAvail);
                        if ((oldJulianDay < beginGrowing) && (julianDay > beginGrowing))
                            if ((julianDay - beginGrowing) < dryDayInterp)
                                dryDayInterp = julianDay - beginGrowing;
    
                        if ((oldJulianDay < endGrowing) && (julianDay > endGrowing))
                            dryDayInterp = endGrowing - julianDay + dryDayInterp;
    
                        if (dryDayInterp < 0.0)
                            dryDayInterp = 0.0;
    
                    } 
                    else if ((priorWaterAvail < wiltingPoint) && (waterAvail > wiltingPoint)) 
                    {
                        dryDayInterp = daysInMonth * (wiltingPoint - priorWaterAvail) / 
                                        (waterAvail - priorWaterAvail);
          
                        if ((oldJulianDay < beginGrowing) && (julianDay > beginGrowing))
                            dryDayInterp = oldJulianDay + dryDayInterp - beginGrowing;
    
                        if (dryDayInterp < 0.0)
                            dryDayInterp = 0.0;
    
                        if ((oldJulianDay < endGrowing) && (julianDay > endGrowing))
                            if ((endGrowing - oldJulianDay) < dryDayInterp)
                                dryDayInterp = endGrowing - oldJulianDay;
                    } 
                    else // ALL below wilting point
                    {
                        dryDayInterp = daysInMonth;
          
                        if ((oldJulianDay < beginGrowing) && (julianDay > beginGrowing))
                            dryDayInterp = julianDay - beginGrowing;
    
                        if ((oldJulianDay < endGrowing) && (julianDay > endGrowing))
                            dryDayInterp = endGrowing - oldJulianDay;
                    }
      
                    dryDays += (int) dryDayInterp;
                }
                return dryDays;
        }
        
        //---------------------------------------------------------------------------

        private static double CalculateDecayFactor(int idef, double soilTemp, double rwc, double ratioPrecipPET, int month)
        {
            // Decomposition factor relfecting the effects of soil temperature and moisture on decomposition
            // Irrigation is zero for natural forests
            double decayFactor = 0.0;   //represents defac in the original program defac.f
            double W_Decomp = 0.0;      //Water effect on decomposition

            //...where
            //      soilTemp;        //Soil temperature
            //      T_Decomp;     //Effect of soil temperature on decomposition
            //      W_Decomp;     //Effect of soil moisture on decompostion
            //      rwcf[10];     //Initial relative water content for 10 soil layers
            //      avh2o;        //Water available to plants for growth in soil profile
            //      precipitation;       //Precipitation of current month
            //      irract;       //Actual amount of irrigation per month (cm H2O/month)
            //      pet;          //Monthly potential evapotranspiration in centimeters (cm)

            //Option selection for wfunc depending on idef
            //      idef = 0;     // for linear option
            //      idef = 1;     // for ratio option


            if (idef == 0)
            {
                if (rwc > 13.0)
                    W_Decomp = 1.0;
                else
                    W_Decomp = 1.0 / (1.0 + 4.0 * System.Math.Exp(-6.0 * rwc));
            }
            else if (idef == 1)
            {
                if (ratioPrecipPET > 9)
                    W_Decomp = 1.0;
                else
                    W_Decomp = 1.0 / (1.0 + 30.0 * System.Math.Exp(-8.5 * ratioPrecipPET));
            }

            double tempModifier = T_Decomp(soilTemp);

            decayFactor = tempModifier * W_Decomp;

            //defac must >= 0.0
            if (decayFactor < 0.0) decayFactor = 0.0;
            if (decayFactor > 1.0) decayFactor = 1.0;

            //if (soilTemp < 0 && decayFactor > 0.01)
            //{
            //    PlugIn.ModelCore.UI.WriteLine("Yr={0},Mo={1}, PET={2:0.00}, MinT={3:0.0}, MaxT={4:0.0}, AveT={5:0.0}, H20={6:0.0}.", Century.Year, month, pet, tmin, tmax, tave, H2Oinputs);
            //    PlugIn.ModelCore.UI.WriteLine("Yr={0},Mo={1}, DecayFactor={2:0.00}, tempFactor={3:0.00}, waterFactor={4:0.00}, ratioPrecipPET={5:0.000}, soilT={6:0.0}.", Century.Year, month, decayFactor, tempModifier, W_Decomp, ratioPrecipPET, soilTemp);
            //}

            return decayFactor;   //Combination of water and temperature effects on decomposition
        }

        //---------------------------------------------------------------------------
        private static double T_Decomp(double soilTemp)
        {
            //Originally from tcalc.f
            //This function computes the effect of temperature on decomposition.
            //It is an exponential function.  Older versions of Century used a density function.
            //Created 10/95 - rm


            double Teff0 = OtherData.TemperatureEffectIntercept;
            double Teff1 = OtherData.TemperatureEffectSlope;
            double Teff2 = OtherData.TemperatureEffectExponent;

            double r = Teff0 + (Teff1 * System.Math.Exp(Teff2 * soilTemp));

            return r;
        }
        //---------------------------------------------------------------------------
        private static double CalculateAnaerobicEffect(double drain, double ratioPrecipPET, double pet, double tave, Site site, int month)
        {

            //Originally from anerob.f of Century

            //...This function calculates the impact of soil anerobic conditions
            //     on decomposition.  It returns a multiplier 'anerob' whose value
            //     is 0-1.

            //...Declaration explanations:
            //     aneref[1] - ratio RAIN/PET with maximum impact
            //     aneref[2] - ratio RAIN/PET with minimum impact
            //     aneref[3] - minimum impact
            //     drain     - percentage of excess water lost by drainage
            //     newrat    - local var calculated new (RAIN+IRRACT+AVH2O[3])/PET ratio
            //     pet       - potential evapotranspiration
            //     rprpet    - actual (RAIN+IRRACT+AVH2O[3])/PET ratio

            double aneref1 = OtherData.RatioPrecipPETMaximum;  //This value is 1.5
            double aneref2 = OtherData.RatioPrecipPETMinimum;   //This value is 3.0
            double aneref3 = OtherData.AnerobicEffectMinimum;   //This value is 0.3

            double anerob = 1.0;

            //...Determine if RAIN/PET ratio is GREATER than the ratio with
            //     maximum impact.

            if ((ratioPrecipPET > aneref1) && (tave > 2.0))
            {
                double xh2o = (ratioPrecipPET - aneref1) * pet * (1.0 - drain);

                if (xh2o > 0)
                {
                    double newrat = aneref1 + (xh2o / pet);
                    double slope = (1.0 - aneref3) / (aneref1 - aneref2);
                    anerob = 1.0 + slope * (newrat - aneref1);
                    //PlugIn.ModelCore.UI.WriteLine("If higher threshold. newrat={0:0.0}, slope={1:0.00}, anerob={2:0.00}", newrat, slope, anerob);      
                }

                if (anerob < aneref3)
                    anerob = aneref3;
                //PlugIn.ModelCore.UI.WriteLine("Lower than threshold. Anaerobic={0}", anerob);      
                SiteVars.MonthlyXh2o[site][Main.Month] = xh2o; // Chihiro
            }
            //PlugIn.ModelCore.UI.WriteLine("ratioPrecipPET={0:0.0}, tave={1:0.00}, pet={2:0.00}, AnaerobicFactor={3:0.00}, Drainage={4:0.00}", ratioPrecipPET, tave, pet, anerob, drain);         
            //PlugIn.ModelCore.UI.WriteLine("Anaerobic Effect = {0:0.00}.", anerob);

            return anerob;
        }
        //---------------------------------------------------------------------------
        private static double CalculateSoilTemp(double tmin, double tmax, double liveBiomass, double litterBiomass, int month)
        {
            // ----------- Calculate Soil Temperature -----------
            double bio = liveBiomass + (OtherData.EffectLitterSoilT * litterBiomass);
            bio = Math.Min(bio, 600.0);

            //...Maximum temperature
            double maxSoilTemp = tmax + (25.4 / (1.0 + 18.0 * Math.Exp(-0.20 * tmax))) * (Math.Exp(OtherData.EffectBiomassMaxSurfT * bio) - 0.13);

            //...Minimum temperature
            double minSoilTemp = tmin + OtherData.EffectBiomassMinSurfT * bio - 1.78;

            //...Average surface temperature
            //...Note: soil temperature used to calculate potential production does not
            //         take into account the effect of snow (AKM)
            double soilTemp = (maxSoilTemp + minSoilTemp) / 2.0;

            //PlugIn.ModelCore.UI.WriteLine("Month={0}, Soil Temperature = {1}.", month+1, soilTemp);

            return soilTemp;
        }
        //--------------------------------------------------------------------------
        public static void Leach(Site site, double baseFlow, double stormFlow)
        {
           
            //  double minlch, double frlech[3], double stream[8], double basef, double stormf)
            //Originally from leach.f of CENTURY model
            //...This routine computes the leaching of inorganic nitrogen (potential for use with phosphorus, and sulfur)
            //...Written 2/92 -rm. Revised on 12/11 by ML
            // ML left out leaching intensity factor.  Cap on MAX leaching (MINLECH/OMLECH3) is poorly defined in CENTURY manual. Added a NO3frac factor to account 
            //for the fact that only NO3 (not NH4) is leached from soils.  

            //...Called From:   SIMSOM

            //...amtlea:    amount leached
            //...linten:    leaching intensity
            //...strm:      storm flow
            //...base:      base flow

            //Outputs:
            //minerl and stream are recomputed
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];
            double waterMove = SiteVars.WaterMovement[site];

            double amtNLeached = 0.0;

            //PlugIn.ModelCore.UI.WriteLine("WaterMove={0:0}, ", waterMove);         
           
         //...waterMove > 0. indicates a saturated water flow out of layer lyr
            if (waterMove > 0.0 && SiteVars.MineralN[site] > 0.0)
            {
                double textureEffect = OtherData.MineralLeachIntercept + OtherData.MineralLeachSlope * SiteVars.SoilPercentSand[site];//ClimateRegionData.PercentSand[ecoregion];
                //double leachIntensity = (1.0 - (OtherData.OMLeachWater - waterMove) / OtherData.OMLeachWater);
                //amtNLeached = textureEffect * SiteVars.MineralN[site] * OtherData.NfracLeachWater * OtherData.NO3frac;
                amtNLeached = textureEffect * SiteVars.MineralN[site] *  OtherData.NO3frac;
                
                //PlugIn.ModelCore.UI.WriteLine("amtNLeach={0:0.0}, textureEffect={1:0.0}, waterMove={2:0.0}, MineralN={3:0.00}", amtNLeached, textureEffect, waterMove, SiteVars.MineralN[site]);      
            }        
            
            double totalNleached = (baseFlow * amtNLeached) + (stormFlow * amtNLeached);
                        
            SiteVars.MineralN[site] -= totalNleached;
            //PlugIn.ModelCore.UI.WriteLine("AfterSoilWaterLeaching. totalNLeach={0:0.0}, MineralN={1:0.00}", totalNleached, SiteVars.MineralN[site]);         

            SiteVars.Stream[site].Nitrogen += totalNleached;
            SiteVars.MonthlyStreamN[site][Main.Month] += totalNleached;
            //PlugIn.ModelCore.UI.WriteLine("AfterSoilWaterLeaching. totalNLeach={0:0.0}, MineralN={1:0.00}", totalNleached, SiteVars.MineralN[site]);        

            // Chihiro; to input into water model
            SiteVars.MonthlyAmtNLeached[site][Main.Month] = amtNLeached;
            SiteVars.MonthlyStreamNitrate[site][Main.Month] = totalNleached;
            return;
        }

    }
}

