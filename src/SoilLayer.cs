//  Author: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Utilities;
using System;

namespace Landis.Extension.Succession.NECN
{
    /// <summary>
    /// </summary>
    public class SoilLayer 
    {
        
        public static void Decompose(ActiveSite site)
        {

            // IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];
            
            //---------------------------------------------------------------------
            // Surface SOM1 decomposes to SOM1 with CO2 lost to respiration.
            
            //double som1c_surface = SiteVars.SOM1surface[site].Carbon;  
                                
            //if (som1c_surface > 0.0000001)
            //{
            //    // Determine C/N ratios for flows to som2
            //    double radds1 = OtherData.SurfaceActivePoolCNIntercept 
            //        + OtherData.SurfaceActivePoolCNSlope * ((som1c_surface / SiteVars.SOM1surface[site].Nitrogen) 
            //        - OtherData.MinCNSurfMicrobes);
                
            //    double ratioCNtoSOM2 = (som1c_surface / SiteVars.SOM1surface[site].Nitrogen) + radds1;
            //    ratioCNtoSOM2 = System.Math.Max(ratioCNtoSOM2, OtherData.SurfaceActivePoolCNMinimum);
        
            //    //Compute total C flow out of surface microbes.
            //    double totalCflow = som1c_surface
            //        * SiteVars.DecayFactor[site]
            //        * PlugIn.Parameters.DecayRateSurf //ClimateRegionData.DecayRateSurf[ecoregion]
            //        * OtherData.MonthAdjust
            //        * OtherData.LitterParameters[(int)LayerType.Surface].DecayRateMicrobes;
                    
            //    // If decomposition can occur, schedule flows associated with respiration
            //    // and decomposition
            //    //if (SiteVars.SOM1surface[site].DecomposePossible(ratioCNtoSOM2, SiteVars.MineralN[site]))
            //    //{   
                    
            //        //CO2 loss - Compute and schedule respiration flows.
            //        double co2loss = totalCflow * OtherData.P1CO2_Surface;
            //        double netCFlow = totalCflow - co2loss;
            //        SiteVars.SOM1surface[site].Respiration(co2loss, site);

            //        // Decompose Surface SOM1 to SOM2
            //        SiteVars.SOM1surface[site].TransferCarbon(SiteVars.SOM1soil[site], netCFlow);
            //        SiteVars.SOM1surface[site].TransferNitrogen(SiteVars.SOM1soil[site], netCFlow, som1c_surface, ratioCNtoSOM2, site);
            //        //PlugIn.ModelCore.UI.WriteLine(".  MineralN={0:0.00}.", SiteVars.MineralN[site]);

            //    //}
            //}


            //---------------------------------------------------------------------
            // Soil SOM1 decomposes to SOM2 and SOM3 with CO2 loss and leaching
            
            double som1c_soil = SiteVars.SOM1soil[site].Carbon;
            //PlugIn.ModelCore.UI.WriteLine("SOM1soil[site].Carbon={0:0.00}", som1c_soil);
            //PlugIn.ModelCore.UI.WriteLine("SiteVars.MineralN = {0:0.00} - pre SOM1.", SiteVars.MineralN[site]);
          
            if (som1c_soil > 0.0000001)
            {

                double r = -0.008314; // gas constant
                double SoilT = 10.0;  //TBD
                double SoilMoisture = 5.0; // TBD
                double SOC = SiteVars.SOM1soil[site].Carbon;
                double SON = SiteVars.SOM1soil[site].Nitrogen;
                double DOC = 0.0020;
                double DON = 0.0011;
                double microbial_C = 1.9703;
                double microbial_N = 0.1970;
                double ec = 0.0339;  // WHAT IS EC??

                // seasonal DOC input
                double A1 = 0.0005; //       #seasonal amplitude
                double w1 = 2 * Math.PI / 4559;
                double dref = 0.0005;     //#reference input
                double t = 1.0;
                double DOC_input = dref + A1 * Math.Sin(w1 * t - Math.PI / 2);  // # mg C cm-3 soil 
                                                                                //double A2 = 0;            //#daily amplitude
                                                                                //double w2 = 2 * Math.PI;
                                                                                // double t = 1:4559; = MINUTES??


                // fitted2009
                double ea_dep = 64.34135015;    //p[1]  #activation energy of SOM depolymerization
                double ea_upt = 60.26134704;    //p[2]  #activation energy of DOC uptake
                double a_dep = 1.05598E+11;     //p[3]   #pre-exponential constant for SOM depolymerization
                double a_upt = 1.07944E+11;     //p[4]  #pre-exponential constant for uptake
                double frac = 0.000378981;      //p[5]  #fraction of unprotected SOM, Magill et al. 2000
                                                // p[6] 0.000466501;  
                                                // p[7] 0.000486085;
                double cnl = 48.81375915;       //p[8] #C:N of litter
                double cns = 28.089739;         //p[9] #C:N of soil
                double cnm = 9.803885526;         //p[10] #C:N of microbial biomass
                double cne = 2.857178317;       //p[11] #C:N of enzymes
                double km_dep = 0.002467667;    //p[12] #half-saturation constant for SOM depolymerization
                double km_upt = 0.289955018;    //p[13] #half-saturation constant for DOC uptake
                double r_ecloss = 0.00098887;  //p[14]                     #enzyme turnover rate
                double r_death = 0.000150496;    //p[15]                      #microbial turnover rate
                double cue = 0.299595251;       //p[16]                          #carbon use efficiency
                double a = 0.50853192;         //p[17]                            #proportion of enzyme pool acting on SOC
                double pconst = 0.481045149;     //p[18]                       #proportion of assimilated C allocated to enzyme production
                double qconst = 0.491011779;    //p[19]                       #proportion of assimilated N allocated to enzyme production
                double mic_to_som = 0.493381459;//p[20]                   #fraction of dead microbial biomass allocated to SOM
                double km_o2 = 0.115473482;     //p[21]                        #Michaelis constant for O2
                double dgas = 1.631550734;      //p[22]                         #diffusion coefficient for O2 in air
                double dliq = 3.135405232;      //p[23]                         #diffusion coefficient for unprotected SOM and DOM in liquid
                double o2airfrac = 0.202587072; //p[24]                    #volume fraction of O2 in air
                double bd = 0.75743956;        //p[25]                           #bulk density
                double pd = 2.50156948;          //p[26]                           #particle density
                double soilMoistureA = -1.92593874;            //p[27]
                double soilMoistureB = 9.774504229;       //p[28]                     
                double sat = 0.492526228;           //p[29]#saturation level

                // Calculated internally or otherwise
                double litter_c = 20;
                double litter_n = -litter_c / cnl;    //#litter N input to SOC
                double doc_input = -99.99;            //#litter C input to DOC
                double don_input = -doc_input / cns;  //litter N input to DOC

                double porosity = 1 - bd / pd;                                              //calculate porosity
                double soilm = soilMoistureA + soilMoistureB * SoilMoisture;                //calculate soil moisture scalar
                soilm = (soilm > sat) ? sat : soilm;                                        //set upper bound on soil moisture (saturation)
                soilm = (soilm < 0.1) ? 0.1 : soilm;                                        //set lower bound on soil moisture
                double o2 = dgas * o2airfrac * Math.Pow((porosity - soilm), (4.0 / 3.0));    //calculate oxygen concentration
                double sol_soc = dliq * Math.Pow(soilm, 3) * frac * SOC;                     //calculate unprotected SOC
                double sol_son = dliq * Math.Pow(soilm, 3) * frac * SON;                     //calculate unprotected SON
                double vmax_dep = a_dep * Math.Exp(-ea_dep / (r * (SoilT + 273)));          //calculate maximum depolymerization rate
                double vmax_upt = a_upt * Math.Exp(-ea_upt / (r * (SoilT + 273)));          //calculate maximum depolymerization rate

                double upt_c = microbial_C * vmax_upt * DOC / (km_upt + DOC) * o2 / (km_o2 + o2); //calculate DOC uptake
                double cmin = upt_c * (1 - cue);                                            //calculate initial C mineralization
                double upt_n = microbial_N * vmax_upt * DON / (km_upt + DON) * o2 / (km_o2 + o2); //calculate DON uptake
                double death_c = r_death * Math.Pow(microbial_C, 2);                         //calculate density-dependent microbial C turnover
                double death_n = r_death * Math.Pow(microbial_N, 2);                         //calculate density-dependent microbial N turnover

                double enz_c = pconst * cue * upt_c;                                        //calculate potential enzyme C production
                double enz_n = qconst * upt_n;                                              //calculate potential enzyme N production
                double eprod = (enz_c / cne >= enz_n) ? enz_n : (enz_c / cne);              //calculate actual enzyme based on Liebig's Law
                double growth_c = (1 - pconst) * (upt_c * cue) + enz_c - cne * eprod;       //calculate potential microbial biomass C growth
                double growth_n = (1 - qconst) * upt_n + enz_n - eprod;                     //calculate potential microbial biomass N growth
                double growth = (growth_c / cnm >= growth_n) ? growth_n : (growth_c / cnm); //calculate actual microbial biomass growth based on Liebig's Law of the minimum (Schimel & Weintraub 2003 SBB)

                double overflow = growth_c - cnm * growth;                                  //calculate overflow metabolism of C
                double nmin = growth_n - growth;                                            //calculate N mineralization

                double dmic_c = cnm * growth - death_c;                                     //calculate change in microbial C pool
                double dmic_n = growth - death_n;                                           //calculate change in microbial N pool

                double eloss = r_ecloss * ec;                                               //calculate enzyme turnover
                double dec = eprod - eloss;                                                 //calculate change in enzyme pool

                double decom_c = vmax_dep * a * ec * sol_soc / (km_dep + sol_soc + ec);     //calculate depolymerization of SOC using ECA kinetics (Tang 2015 GMD)
                double decom_n = vmax_dep * (1 - a) * ec * sol_son / (km_dep + sol_son + ec); //calculate depolymerization of SON using ECA kinetics 


                double dsoc = litter_c + death_c * mic_to_som - decom_c;                    //calculate change in SOC pool
                double dson = litter_n + death_n * mic_to_som - decom_n;                    //calculate change in SON pool
                double ddoc = doc_input + decom_c + death_c * (1 - mic_to_som) + cne * eloss - upt_c; //calculate change in DOC pool
                double ddon = don_input + decom_n + death_n * (1 - mic_to_som) + eloss - upt_n; //calculate change in DON pool
                double dcout = cmin + overflow;                                             //calculate C efflux

                double co2loss = dcout;
                //Determine C/N ratios for flows to som2
                //double ratioCNtoSOM2  = Layer.BelowgroundDecompositionRatio(site,
                //                            OtherData.MinCNenterSOM2, 
                //                            OtherData.MaxCNenterSOM2,
                //                            OtherData.MinContentN_SOM2);

                //Compute total C flow out of soil microbes.
                //Added impact of soil anaerobic conditions -rm 12/91
                //double textureEffect = OtherData.TextureEffectIntercept + OtherData.TextureEffectSlope * SiteVars.SoilPercentSand[site];
                
                //double anerb = SiteVars.AnaerobicEffect[site];

                ////PlugIn.ModelCore.UI.WriteLine("SiteVars.DecayFactor = {0:0.00}, SoilDecayRateMicrobes = {1:0.00}, texture = {2:0.00}, anerb = {3:0.00}, MonthAdjust = {4:0.00}.",
                //double totalCflow = som1c_soil 
                //            * SiteVars.DecayFactor[site]
                //            * OtherData.LitterParameters[(int) LayerType.Soil].DecayRateMicrobes
                //            * PlugIn.Parameters.DecayRateSOM1 
                //             * textureEffect
                //             * anerb
                //             * OtherData.MonthAdjust;

                // First determine if decomposition can occur:
                //if (SiteVars.SOM1soil[site].DecomposePossible(ratioCNtoSOM2, SiteVars.MineralN[site]))
                //{   
                    //CO2 Loss - Compute and schedule respiration flows
                    //double P1CO2_Soil = OtherData.P1CO2_Soil_Intercept + OtherData.P1CO2_Soil_Slope * SiteVars.SoilPercentSand[site];

                    //double co2loss = totalCflow * P1CO2_Soil;
                    //double netCFlow = totalCflow - co2loss;


                SiteVars.SOM1soil[site].Respiration(co2loss, site);
 
                    // Decompose Soil SOM1 to SOM3
                    // The fraction of totalCflow that goes to SOM3 is a function of clay content.
                    //double clayEffect = OtherData.PS1S3_Intercept + (OtherData.PS1S3_Slope * SiteVars.SoilPercentClay[site]);//ClimateRegionData.PercentClay[ecoregion]);
                    //double cFlowS1S3 = netCFlow * clayEffect * (1.0 + OtherData.AnaerobicImpactSlope * (1.0 - anerb));

                    ////Compute and schedule C & N flows and update mineralization accumulators
                    //double ratioCNto3 = Layer.BelowgroundDecompositionRatio(site,
                    //                        OtherData.MinCNenterSOM3, 
                    //                        OtherData.MaxCNenterSOM3,
                    //                        OtherData.MinContentN_SOM3);
                     
                    ////Partition and schedule C and N flows 
                    //SiteVars.SOM1soil[site].TransferCarbon(SiteVars.SOM3[site], cFlowS1S3);
                    //SiteVars.SOM1soil[site].TransferNitrogen(SiteVars.SOM3[site], cFlowS1S3, som1c_soil, ratioCNto3, site);
                    //PlugIn.ModelCore.UI.WriteLine("AfterSOM1.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
                     
                    // Leaching of Organics
                    // This only occurs when the water flow out of water layer 2
                    // exceeds a critical value.  Use the same C/N ratios as for the flow to SOM3.

                    double cLeached = 0.0;  // Carbon leached to a stream
                    
                    if(SiteVars.WaterMovement[site] > 0.0)  //Volume of water moving-ML.  Used to be an index of water movement that indicates saturation (amov)
                    {

                        double leachTextureEffect = OtherData.OMLeachIntercept + OtherData.OMLeachSlope * SiteVars.SoilPercentSand[site];

                        double indexWaterMovement = SiteVars.WaterMovement[site] / (SiteVars.SoilDepth[site] * SiteVars.SoilFieldCapacity[site]);

                    //cLeached = netCFlow * leachTextureEffect * indexWaterMovement;
                    cLeached = co2loss * leachTextureEffect * indexWaterMovement;

                    //Partition and schedule C flows 
                    SiteVars.SOM1soil[site].TransferCarbon(SiteVars.Stream[site], cLeached);

                        // Compute and schedule N flows and update mineralization accumulators
                        // Need to use the ratio for som1 for organic leaching
                        double ratioCN_SOM1soil = som1c_soil / SiteVars.SOM1soil[site].Nitrogen;
                        double orgflow = cLeached / ratioCN_SOM1soil;

                        SiteVars.SOM1soil[site].Nitrogen -= orgflow; 
                        SiteVars.Stream[site].Nitrogen += orgflow;

                        SiteVars.MonthlyStreamN[site][Main.Month] += orgflow;


                    }

                    // C & N movement from SOM1 to SOM2.
                    // SOM2 gets what's left of totalCflow.
                    //double cFlowS1S2 = netCFlow - cFlowS1S3 - cLeached;

                    ////Partition and schedule C and N flows 
                    //SiteVars.SOM1soil[site].TransferCarbon(SiteVars.SOM2[site], cFlowS1S2);
                    //SiteVars.SOM1soil[site].TransferNitrogen(SiteVars.SOM2[site], cFlowS1S2, som1c_soil, ratioCNtoSOM2, site);
                    //PlugIn.ModelCore.UI.WriteLine("PartitionCN.  MineralN={0:0.00}.", SiteVars.MineralN[site]);

                //}  
            } 


            //---------------------------------------------------------------------
            //**********SOM2 decomposes to soil SOM1 and SOM3 with CO2 loss**********

            //double som2c = SiteVars.SOM2[site].Carbon;
          
            //if (som2c > 0.0000001)
            //{
            //    // Determine C/N ratios for flows to SOM1
            //    double ratioCNto1 = Layer.BelowgroundDecompositionRatio(site,
            //                            OtherData.MinCNenterSOM1, 
            //                            OtherData.MaxCNenterSOM1,
            //                            OtherData.MinContentN_SOM1);
                
            //    double anerb = SiteVars.AnaerobicEffect[site];  

            //    // Compute total C flow out of SOM2C
            //    double totalCflow = som2c 
            //                    * SiteVars.DecayFactor[site] 
            //                    * PlugIn.Parameters.DecayRateSOM2 
            //                    * anerb //impact of soil anaerobic conditions
            //                    * OtherData.MonthAdjust;
            //    //PlugIn.ModelCore.UI.WriteLine("som2c={0:0.00}, decayFactor={1:0.00}, decayRateSOM2={2:0.00}, anerb={3:0.00}, monthAdj={4:0.00}", som2c, SiteVars.DecayFactor[site], ClimateRegionData.DecayRateSOM2[ecoregion], anerb, OtherData.MonthAdjust);

            //    // If SOM2 can decompose to SOM1, it will also go to SOM3.
            //    // If it can't go to SOM1, it can't decompose at all.

            //    if (SiteVars.SOM2[site].DecomposePossible(ratioCNto1, SiteVars.MineralN[site]))
            //        //PlugIn.ModelCore.UI.WriteLine("DecomposePoss.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
            //    {
                
            //        //CO2 loss - Compute and schedule respiration flows
            //        double co2loss = totalCflow * OtherData.FractionSOM2toCO2;
            //        double netCFlow = totalCflow - co2loss;
            //        SiteVars.SOM2[site].Respiration(co2loss, site);
            //        //PlugIn.ModelCore.UI.WriteLine("AfterTransferto.  MineralN={0:0.00}.", SiteVars.MineralN[site]);

            //        // -----------------------------------------------
            //        // Decompose SOM2 to SOM3, SOM3 gets what's left of totalCflow.
            //        double clayEffect = OtherData.PS2S3_Intercept + OtherData.PS2S3_Slope * SiteVars.SoilPercentClay[site];//ClimateRegionData.PercentClay[ecoregion];
            //        double cFlowS2S3 = netCFlow * clayEffect * (1.0 + OtherData.AnaerobicImpactSlope * (1.0 - anerb));

            //        //Compute and schedule C and N flows and update mineralization accumulators
            //        double ratioCNto3 = Layer.BelowgroundDecompositionRatio(site,
            //                                OtherData.MinCNenterSOM3, 
            //                                OtherData.MaxCNenterSOM3,
            //                                OtherData.MinContentN_SOM3);
            //        //PlugIn.ModelCore.UI.WriteLine("TransferSOM2.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
                    
            //        //Partition and schedule C and N flows 
            //        SiteVars.SOM2[site].TransferCarbon(SiteVars.SOM3[site], cFlowS2S3);
            //        SiteVars.SOM2[site].TransferNitrogen(SiteVars.SOM3[site], cFlowS2S3, som2c, ratioCNto3, site);
                   
            //        // -----------------------------------------------
            //        // Decompose SOM2 to SOM1
            //        double cFlowS2S1 = netCFlow - cFlowS2S3;

            //        // Compute and schedule N and C flows and update mineralization accumulators
            //        ratioCNto1 = Layer.BelowgroundDecompositionRatio(site,
            //                            OtherData.MinCNenterSOM1, 
            //                            OtherData.MaxCNenterSOM1,
            //                            OtherData.MinContentN_SOM1);

            //        //Partition and schedule C and N flows 
            //        SiteVars.SOM2[site].TransferCarbon(SiteVars.SOM1soil[site], cFlowS2S1);
            //        SiteVars.SOM2[site].TransferNitrogen(SiteVars.SOM1soil[site], cFlowS2S1, som2c, ratioCNto1, site);
            //        //PlugIn.ModelCore.UI.WriteLine("AfterSOM2.  MineralN={0:0.00}.", SiteVars.MineralN[site]);
            //    }
                
            //}

            //---------------------------------------------------------------------
            // SOM3 decomposes to soil SOM1 with CO2 loss
           
            //double som3c = SiteVars.SOM3[site].Carbon; 
            
            //if (som3c > 0.0000001)
            //{
            //    //Determine C/N ratios for flows to SOM1.
            //    double ratioCNto1 = Layer.BelowgroundDecompositionRatio(site,
            //                            OtherData.MinCNenterSOM1, 
            //                            OtherData.MaxCNenterSOM1,
            //                            OtherData.MinContentN_SOM1);
                 
            //    double anerb = SiteVars.AnaerobicEffect[site];  

            //    //Compute total C flow out of SOM3C
            //    double totalCflow = som3c
            //                    * SiteVars.DecayFactor[site]
            //                    * PlugIn.Parameters.DecayRateSOM3 
            //                    * anerb 
            //                    * OtherData.MonthAdjust;


            //    //If decomposition can occur,
            //    if (SiteVars.SOM3[site].DecomposePossible(ratioCNto1, SiteVars.MineralN[site]))
            //    {
            //        //PlugIn.ModelCore.UI.WriteLine("BeforeSOM3 Decay.  C={0:0.00}.", SiteVars.SOM3[site].Carbon);
            //        // CO2 loss - Compute and schedule respiration flows.
            //        double co2loss = totalCflow * OtherData.FractionSOM3toCO2 * anerb;
            //        double netCFlow = totalCflow - co2loss;
            //        SiteVars.SOM3[site].Respiration(co2loss, site);

            //        // Decompose SOM3 to soil SOM1
            //        double cFlowS3S1 = netCFlow;

            //        // Partition and schedule C and N flows 
            //        SiteVars.SOM3[site].TransferCarbon(SiteVars.SOM1soil[site], cFlowS3S1);
            //        SiteVars.SOM3[site].TransferNitrogen(SiteVars.SOM1soil[site], cFlowS3S1, som3c, ratioCNto1, site);
            //        //PlugIn.ModelCore.UI.WriteLine("AfterSOM3 Decay.  C={0:0.00}.", SiteVars.SOM3[site].Carbon);
            //    }
            //}
        }
    }
}
