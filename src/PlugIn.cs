//  Authors: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Utilities;

using Landis.Library.InitialCommunities;
using Landis.Library.Succession;
using Landis.Library.LeafBiomassCohorts;
using Landis.Library.Climate;
using Landis.Library.Metadata;

using System;
using System.Collections.Generic;
using System.Linq;

namespace Landis.Extension.Succession.NECN
{
    public class PlugIn
        : Landis.Library.Succession.ExtensionBase
    {
        public static readonly string ExtensionName = "NECN Succession";
        private static ICore modelCore;
        public static IInputParameters Parameters;
        public static double[] ShadeLAI;
        public static double AnnualWaterBalance;

        private List<ISufficientLight> sufficientLight;
        public static string SoilCarbonMapNames = null;
        public static int SoilCarbonMapFrequency;
        public static string SoilNitrogenMapNames = null;
        public static int SoilNitrogenMapFrequency;
        public static string ANPPMapNames = null;
        public static int ANPPMapFrequency;
        public static string ANEEMapNames = null;
        public static int ANEEMapFrequency;
        public static string TotalCMapNames = null;
        public static int TotalCMapFrequency;
        public static int SuccessionTimeStep;
        public static double ProbEstablishAdjust;

        public static int FutureClimateBaseYear;
        public static int B_MAX;
        private ICommunity initialCommunity;

        //---------------------------------------------------------------------

        public PlugIn()
            : base(ExtensionName)
        {
        }

        //---------------------------------------------------------------------

        public override void LoadParameters(string dataFile,
                                            ICore mCore)
        {
            modelCore = mCore;
            //SiteVars.Initialize(); // this use functional type parameter
            InputParametersParser parser = new InputParametersParser();
            Parameters = Landis.Data.Load<IInputParameters>(dataFile, parser);

        }

        //---------------------------------------------------------------------

        public static ICore ModelCore
        {
            get
            {
                return modelCore;
            }
        }


        //---------------------------------------------------------------------

        public override void Initialize()
        {
            PlugIn.ModelCore.UI.WriteLine("Initializing {0} ...", ExtensionName);
            Timestep = Parameters.Timestep;
            SuccessionTimeStep = Timestep;
            sufficientLight = Parameters.LightClassProbabilities;
            ProbEstablishAdjust = Parameters.ProbEstablishAdjustment;
            MetadataHandler.InitializeMetadata(Timestep, modelCore, SoilCarbonMapNames, SoilNitrogenMapNames, ANPPMapNames, ANEEMapNames, TotalCMapNames); //,LAIMapNames, ShadeClassMapNames);

            FunctionalType.Initialize(Parameters);
            SpeciesData.Initialize(Parameters);
            SiteVars.Initialize(); // chihiro; this method use functional type data for initializing decay value
            Util.ReadSoilDepthMap(Parameters.SoilDepthMapName);
            Util.ReadSoilDrainMap(Parameters.SoilDrainMapName);
            Util.ReadSoilBaseFlowMap(Parameters.SoilBaseFlowMapName);
            Util.ReadSoilStormFlowMap(Parameters.SoilStormFlowMapName);
            Util.ReadFieldCapacityMap(Parameters.SoilFieldCapacityMapName);
            Util.ReadWiltingPointMap(Parameters.SoilWiltingPointMapName);
            Util.ReadPercentSandMap(Parameters.SoilPercentSandMapName);
            Util.ReadPercentClayMap(Parameters.SoilPercentClayMapName);
            Util.ReadSoilCNMaps(Parameters.InitialSOM1CSurfaceMapName,
                Parameters.InitialSOM1NSurfaceMapName,
                Parameters.InitialSOM1CSoilMapName,
                Parameters.InitialSOM1NSoilMapName,
                Parameters.InitialSOM2CMapName,
                Parameters.InitialSOM2NMapName,
                Parameters.InitialSOM3CMapName,
                Parameters.InitialSOM3NMapName);
            Util.ReadDeadWoodMaps(Parameters.InitialDeadSurfaceMapName, Parameters.InitialDeadSoilMapName);

            //Initialize climate.
            Climate.Initialize(Parameters.ClimateConfigFile, false, modelCore);
            FutureClimateBaseYear = Climate.Future_MonthlyData.Keys.Min();
            ClimateRegionData.Initialize(Parameters);

            ShadeLAI = Parameters.MaximumShadeLAI;
            OtherData.Initialize(Parameters);
            FireEffects.Initialize(Parameters);

            //  Cohorts must be created before the base class is initialized
            //  because the base class' reproduction module uses the core's
            //  SuccessionCohorts property in its Initialization method.
            Library.LeafBiomassCohorts.Cohorts.Initialize(Timestep, new CohortBiomass());

            // Initialize Reproduction routines:
            Reproduction.SufficientResources = SufficientLight;
            Reproduction.Establish = Establish;
            Reproduction.AddNewCohort = AddNewCohort;
            Reproduction.MaturePresent = MaturePresent;
            base.Initialize(modelCore, Parameters.SeedAlgorithm);

            // Delegate mortality routines:
            Landis.Library.LeafBiomassCohorts.Cohort.PartialDeathEvent += CohortPartialMortality;
            Landis.Library.LeafBiomassCohorts.Cohort.DeathEvent += CohortTotalMortality;

            InitializeSites(Parameters.InitialCommunities, Parameters.InitialCommunitiesMap, modelCore);

            if (Parameters.CalibrateMode)
                Outputs.CreateCalibrateLogFile();
            //Establishment.InitializeLogFile();

            B_MAX = 0;
            foreach (ISpecies species in ModelCore.Species)
            {
                if (SpeciesData.Max_Biomass[species] > B_MAX)
                    B_MAX = SpeciesData.Max_Biomass[species];
            }

            foreach (ActiveSite site in PlugIn.ModelCore.Landscape)
                Main.ComputeTotalCohortCN(site, SiteVars.Cohorts[site]);

            Outputs.WritePrimaryLogFile(0);
            Outputs.WriteShortPrimaryLogFile(0);


        }

        //---------------------------------------------------------------------

        public override void Run()
        {

            if (PlugIn.ModelCore.CurrentTime > 0)
                SiteVars.InitializeDisturbances();

            ClimateRegionData.AnnualNDeposition = new Landis.Library.Parameters.Ecoregions.AuxParm<double>(PlugIn.ModelCore.Ecoregions);

            //base.RunReproductionFirst();

            base.Run();

            if (Timestep > 0)
                ClimateRegionData.SetAllEcoregions_FutureAnnualClimate(ModelCore.CurrentTime);

            if (ModelCore.CurrentTime % Timestep == 0)
            {
                // Write monthly log file:
                // Output must reflect the order of operation:
                int[] months = new int[12] { 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5 };

                if (OtherData.CalibrateMode)
                    months = new int[12] { 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5 };

                for (int i = 0; i < 12; i++)
                {
                    int month = months[i];
                    Outputs.WriteMonthlyLogFile(month);
                }
                Outputs.WritePrimaryLogFile(PlugIn.ModelCore.CurrentTime);
                Outputs.WriteShortPrimaryLogFile(PlugIn.ModelCore.CurrentTime);
                // Outputs.WriteMaps(); Wataru, 2020.09.03
                if (!OtherData.ReduceOutputs) // Wataru, 2020.09.03 
                    Outputs.WriteMaps();      // Wataru, 2020.09.03
                Establishment.LogEstablishment();
            }

        }


        //---------------------------------------------------------------------

        public override byte ComputeShade(ActiveSite site)
        {
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];

            byte finalShade = 0;

            if (!ecoregion.Active)
                return 0;

            for (byte shade = 5; shade >= 1; shade--)
            {
                if (PlugIn.ShadeLAI[shade] <= 0)
                {
                    string mesg = string.Format("Maximum LAI has not been defined for shade class {0}", shade);
                    throw new System.ApplicationException(mesg);
                }
                if (SiteVars.LAI[site] >= PlugIn.ShadeLAI[shade])
                {
                    finalShade = shade;
                    break;
                }
            }

            //PlugIn.ModelCore.UI.WriteLine("Yr={0},      Shade Calculation:  B_MAX={1}, B_ACT={2}, Shade={3}.", PlugIn.ModelCore.CurrentTime, B_MAX, B_ACT, finalShade);

            return finalShade;
        }
        //---------------------------------------------------------------------


        /// <summary>
        /// Compute site shade above the grass species
        /// </summary>
        // Chihiro 2020.01.22
        //
        // Description:
        //     We just copied ComputeShade method and replaced SiteVars.LAI[site] to SiteVars.LAITree[site]
        //     to ignore grass species.
        //
        public static byte ComputeShadeTree(ActiveSite site)
        {
            IEcoregion ecoregion = PlugIn.ModelCore.Ecoregion[site];
            byte finalShade = 0;
            if (!ecoregion.Active)
                return 0;
            for (byte shade = 5; shade >= 1; shade--)
            {
                if (PlugIn.ShadeLAI[shade] <= 0)
                {
                    string mesg = string.Format("Maximum LAI has not been defined for shade class {0}", shade);
                    throw new System.ApplicationException(mesg);
                }
                if (SiteVars.LAITree[site] >= PlugIn.ShadeLAI[shade])
                {
                    finalShade = shade;
                    break;
                }
            }
            return finalShade;
        }
        //---------------------------------------------------------------------


        protected override void InitializeSite(ActiveSite site)//, ICommunity initialCommunity)
        {

            InitialBiomass initialBiomass = InitialBiomass.Compute(site, initialCommunity);
            SiteVars.MineralN[site] = Parameters.InitialMineralN;
        }


        //---------------------------------------------------------------------
        // This method does not trigger reproduction
        public void CohortPartialMortality(object sender, Landis.Library.BiomassCohorts.PartialDeathEventArgs eventArgs)
        {
            //PlugIn.ModelCore.UI.WriteLine("Cohort Partial Mortality:  {0}", eventArgs.Site);

            ExtensionType disturbanceType = eventArgs.DisturbanceType;
            ActiveSite site = eventArgs.Site;

            ICohort cohort = (Landis.Library.LeafBiomassCohorts.ICohort)eventArgs.Cohort;

            float fractionPartialMortality = (float)eventArgs.Reduction;
            float foliarInput = cohort.LeafBiomass * fractionPartialMortality;
            float woodInput = cohort.WoodBiomass * fractionPartialMortality;

            if (disturbanceType.IsMemberOf("disturbance:harvest"))
            {
                SiteVars.HarvestPrescriptionName = PlugIn.ModelCore.GetSiteVar<string>("Harvest.PrescriptionName");
                if (!Disturbed[site]) // this is the first cohort killed/damaged
                {
                    HarvestEffects.ReduceLayers(SiteVars.HarvestPrescriptionName[site], site);
                }
                woodInput -= woodInput * (float)HarvestEffects.GetCohortWoodRemoval(site);
                foliarInput -= foliarInput * (float)HarvestEffects.GetCohortLeafRemoval(site);
            }
            if (disturbanceType.IsMemberOf("disturbance:fire"))
            {

                SiteVars.FireSeverity = PlugIn.ModelCore.GetSiteVar<byte>("Fire.Severity");

                if (!Disturbed[site]) // this is the first cohort killed/damaged
                {
                    SiteVars.SmolderConsumption[site] = 0.0;
                    SiteVars.FlamingConsumption[site] = 0.0;
                    if (SiteVars.FireSeverity != null && SiteVars.FireSeverity[site] > 0)
                        FireEffects.ReduceLayers(SiteVars.FireSeverity[site], site);

                }

                double woodFireConsumption = woodInput * (float)FireEffects.ReductionsTable[(int)SiteVars.FireSeverity[site]].CoarseLitterReduction;
                double foliarFireConsumption = foliarInput * (float)FireEffects.ReductionsTable[(int)SiteVars.FireSeverity[site]].FineLitterReduction;

                SiteVars.SmolderConsumption[site] += woodFireConsumption;
                SiteVars.FlamingConsumption[site] += foliarFireConsumption;
                woodInput -= (float)woodFireConsumption;
                foliarInput -= (float)foliarFireConsumption;
            }

            ForestFloor.AddWoodLitter(woodInput, cohort.Species, site);
            ForestFloor.AddFoliageLitter(foliarInput, cohort.Species, site);

            Roots.AddCoarseRootLitter(Roots.CalculateCoarseRoot(cohort, cohort.WoodBiomass * fractionPartialMortality), cohort, cohort.Species, site);
            Roots.AddFineRootLitter(Roots.CalculateFineRoot(cohort, cohort.LeafBiomass * fractionPartialMortality), cohort, cohort.Species, site);

            //PlugIn.ModelCore.UI.WriteLine("EVENT: Cohort Partial Mortality: species={0}, age={1}, disturbance={2}.", cohort.Species.Name, cohort.Age, disturbanceType);
            //PlugIn.ModelCore.UI.WriteLine("       Cohort Reductions:  Foliar={0:0.00}.  Wood={1:0.00}.", HarvestEffects.GetCohortLeafRemoval(site), HarvestEffects.GetCohortLeafRemoval(site));
            //PlugIn.ModelCore.UI.WriteLine("       InputB/TotalB:  Foliar={0:0.00}/{1:0.00}, Wood={2:0.0}/{3:0.0}.", foliarInput, cohort.LeafBiomass, woodInput, cohort.WoodBiomass);
            Disturbed[site] = true;

            return;
        }
        //---------------------------------------------------------------------
        // Total mortality, including from disturbance or senescence.

        public void CohortTotalMortality(object sender, Landis.Library.BiomassCohorts.DeathEventArgs eventArgs)
        {

            //PlugIn.ModelCore.UI.WriteLine("Cohort Total Mortality: {0}", eventArgs.Site);

            ExtensionType disturbanceType = eventArgs.DisturbanceType;
            ActiveSite site = eventArgs.Site;

            ICohort cohort = (Landis.Library.LeafBiomassCohorts.ICohort)eventArgs.Cohort;
            double foliarInput = (double)cohort.LeafBiomass;
            double woodInput = (double)cohort.WoodBiomass;

            if (disturbanceType != null)
            {
                //PlugIn.ModelCore.UI.WriteLine("DISTURBANCE EVENT: Cohort Died: species={0}, age={1}, disturbance={2}.", cohort.Species.Name, cohort.Age, eventArgs.DisturbanceType);

                if (disturbanceType.IsMemberOf("disturbance:fire"))
                {
                    SiteVars.FireSeverity = PlugIn.ModelCore.GetSiteVar<byte>("Fire.Severity");
                    Landis.Library.Succession.Reproduction.CheckForPostFireRegen(eventArgs.Cohort, site);

                    if (!Disturbed[site])  // the first cohort killed/damaged
                    {
                        SiteVars.SmolderConsumption[site] = 0.0;
                        SiteVars.FlamingConsumption[site] = 0.0;
                        if (SiteVars.FireSeverity != null && SiteVars.FireSeverity[site] > 0)
                            FireEffects.ReduceLayers(SiteVars.FireSeverity[site], site);

                    }

                    double woodFireConsumption = woodInput * (float)FireEffects.ReductionsTable[(int)SiteVars.FireSeverity[site]].CoarseLitterReduction;
                    double foliarFireConsumption = foliarInput * (float)FireEffects.ReductionsTable[(int)SiteVars.FireSeverity[site]].FineLitterReduction;

                    SiteVars.SmolderConsumption[site] += woodFireConsumption;
                    SiteVars.FlamingConsumption[site] += foliarFireConsumption;
                    SiteVars.SourceSink[site].Carbon += woodFireConsumption * 0.47;
                    SiteVars.SourceSink[site].Carbon += foliarFireConsumption * 0.47;
                    woodInput -= woodFireConsumption;
                    foliarInput -= foliarFireConsumption;
                }
                else
                {
                    if (disturbanceType.IsMemberOf("disturbance:harvest"))
                    {
                        SiteVars.HarvestPrescriptionName = PlugIn.ModelCore.GetSiteVar<string>("Harvest.PrescriptionName");
                        if (!Disturbed[site])  // the first cohort killed/damaged
                        {
                            HarvestEffects.ReduceLayers(SiteVars.HarvestPrescriptionName[site], site);
                        }
                        double woodLoss = woodInput * (float)HarvestEffects.GetCohortWoodRemoval(site);
                        double foliarLoss = foliarInput * (float)HarvestEffects.GetCohortLeafRemoval(site);
                        SiteVars.SourceSink[site].Carbon += woodLoss * 0.47;
                        SiteVars.SourceSink[site].Carbon += foliarLoss * 0.47;
                        woodInput -= woodLoss;
                        foliarInput -= foliarLoss;
                    }

                    // If not fire, check for resprouting:
                    Landis.Library.Succession.Reproduction.CheckForResprouting(eventArgs.Cohort, site);
                }
            }


            //PlugIn.ModelCore.UI.WriteLine("Cohort Died: species={0}, age={1}, wood={2:0.00}, foliage={3:0.00}.", cohort.Species.Name, cohort.Age, wood, foliar);
            ForestFloor.AddWoodLitter(woodInput, cohort.Species, eventArgs.Site);
            ForestFloor.AddFoliageLitter(foliarInput, cohort.Species, eventArgs.Site);

            // Assume that ALL dead root biomass stays on site.
            Roots.AddCoarseRootLitter(Roots.CalculateCoarseRoot(cohort, cohort.WoodBiomass), cohort, cohort.Species, eventArgs.Site);
            Roots.AddFineRootLitter(Roots.CalculateFineRoot(cohort, cohort.LeafBiomass), cohort, cohort.Species, eventArgs.Site);

            if (disturbanceType != null)
                Disturbed[site] = true;

            return;
        }

        //---------------------------------------------------------------------
        //Grows the cohorts for future climate
        protected override void AgeCohorts(ActiveSite site,
                                           ushort years,
                                           int? successionTimestep)
        {
            Main.Run(site, years, successionTimestep.HasValue);

        }
        //---------------------------------------------------------------------
        /// <summary>
        /// Determines if there is sufficient light at a site for a species to
        /// germinate/resprout.
        /// This is a Delegate method to base succession.
        /// </summary>
        /// 
        // W.Hotta and Chihiro modified
        // 
        // Description:
        //     - Modify light probability based on the amount of nursery log on the site
        //     - Return both True/False and established location (surface, nlog or surface)
        //
        //
        // Psudo-code:
        // 
        //     double siteShade = site shade calculated by original NECN // This siteShade includes both tree and grass species
        //     double siteShadeTree = ComputeShadeTree(site)             // Only tree species, excluding grass species
        //     bool isSufficientlight = false
        //     
        //     double lightProbability = light probability of original NECN // This lightProbability includes both tree and grass species
        //     double lightProbabilityTree = light probability ignoring grass layer
        //     bool found = false
        //     
        //     int bestShadeClass = 0      // the most suitable shade class for the species
        //     string regenType = 'failed' // Identify where the cohort established
        //                                 // 'failed', 'nlog', or 'surface'
        //     
        //     // Compute the light probabilities on the site
        //     foreach lights in sufficientLight:
        //         identify lightProbability, lightProbabilityTree, and bestShadeClass for the species on the site
        //     
        //     // Compute the availability of nursery log on the site (nurseryLogAvailability)
        //     double nurseryLogAvailabilityModifier = 1.0;                    // user defined tuning parameter
        //     double nurseryLogAvailability = nurseryLogAvailabilityModifier 
        //                                     * ComputeNurseryLogAreaRatio(species, site);
        //     
        //     // Determine if the species can establish or not
        //     if The species is CWD-dependent: (Case 1)
        //         if there is sufficient light & enough nursery logs:
        //             // establish on the nursery log
        //             isSufficientlight = true
        //             regenType = 'nlog'
        //     
        //     else The species can establish on both forest floors & nursery logs (CWD independent): (Case 2)
        //         if there is sufficient light on forest floor considering both Tree and Grass species LAI:
        //             // establish on the surface soil
        //             isSufficientlight = true
        //             regenType = 'surface'
        //             
        //         else If (1) the site shade is darker than the best shade class for the species and 
        //                 (2) the light availability above grass species layer meets the species requirement:
        //             if Threre are sufficient amounts of downed logs:
        //                 // establish on the nursery log
        //                 isSufficientlight = true
        //                 regenType = 'nlog'
        //
        //     return isSufficientlight;
        //
        //
        public (bool, string) SufficientLight(ISpecies species, ActiveSite site) // 2020.10.31 Chhiro
        {

            //PlugIn.ModelCore.UI.WriteLine("  Calculating Sufficient Light from Succession.");
            byte siteShade = PlugIn.ModelCore.GetSiteVar<byte>("Shade")[site];
            byte siteShadeTree = ComputeShadeTree(site); // Shade by tree species excluding grass species; Chihiro
            bool isSufficientlight = false;

            double lightProbability = 0.0;
            double lightProbabilityTree = 0.0; // Light probability considering LAI of tree species only; Chihiro
            bool found = false;

            int bestShadeClass = 0; // the best shade class for the species; Chihiro
            string regenType = "failed"; // Identify where the cohort established; Chihiro

            foreach (ISufficientLight lights in sufficientLight)
            {

                //PlugIn.ModelCore.UI.WriteLine("Sufficient Light:  ShadeClass={0}, Prob0={1}.", lights.ShadeClass, lights.ProbabilityLight0);
                if (lights.ShadeClass == species.ShadeTolerance)
                {
                    if (siteShade == 0) lightProbability = lights.ProbabilityLight0;
                    if (siteShade == 1) lightProbability = lights.ProbabilityLight1;
                    if (siteShade == 2) lightProbability = lights.ProbabilityLight2;
                    if (siteShade == 3) lightProbability = lights.ProbabilityLight3;
                    if (siteShade == 4) lightProbability = lights.ProbabilityLight4;
                    if (siteShade == 5) lightProbability = lights.ProbabilityLight5;

                    // Light probability considering LAI of only tree species; Chihiro
                    if (siteShadeTree == 0) lightProbabilityTree = lights.ProbabilityLight0;
                    if (siteShadeTree == 1) lightProbabilityTree = lights.ProbabilityLight1;
                    if (siteShadeTree == 2) lightProbabilityTree = lights.ProbabilityLight2;
                    if (siteShadeTree == 3) lightProbabilityTree = lights.ProbabilityLight3;
                    if (siteShadeTree == 4) lightProbabilityTree = lights.ProbabilityLight4;
                    if (siteShadeTree == 5) lightProbabilityTree = lights.ProbabilityLight5;

                    // Identify the best shade class for the species; Chihiro
                    bestShadeClass = ComputeBestShadeClass(lights);
                    // if (PlugIn.ModelCore.CurrentTime == 1)
                    //     PlugIn.ModelCore.UI.WriteLine("MaxLight:{0},{1}",species.Name, bestShadeClass);

                    found = true;
                }
            }

            if (!found)
                PlugIn.ModelCore.UI.WriteLine("A Sufficient Light value was not found for {0}.", species.Name);

            // ------------------------------------------------------------------------
            // Modify light probability based on the amount of nursery log on the site
            // W.Hotta 2020.01.22
            //
            // Compute the availability of nursery log on the site
            //   Option1: function type is linear
            // double nurseryLogAvailabilityModifier = 1.0; // tuning parameter
            // double nurseryLogAvailability = nurseryLogAvailabilityModifier * ComputeNurseryLogAreaRatio(species, site);
            //   Option2: function type is power
            double nurseryLogAvailabilityModifier = 2.0; // tuning parameter (only even)
            double nurseryLogAvailability = 1 - Math.Pow(ComputeNurseryLogAreaRatio(species, site) - 1, nurseryLogAvailabilityModifier);
            // if (OtherData.CalibrateMode)
            if (OtherData.CalibrateMode && !OtherData.ReduceOutputs) //Wataru, 2020.09.03
            {
                PlugIn.ModelCore.UI.WriteLine("original_lightProbability:{0},{1},{2}", PlugIn.ModelCore.CurrentTime, species.Name, lightProbability);
                PlugIn.ModelCore.UI.WriteLine("siteShade:{0},{1}", siteShade, siteShadeTree);
                PlugIn.ModelCore.UI.WriteLine("siteLAI:{0}", SiteVars.LAI[site]);
            }

            // Case 1. CWD-dependent species (species which can only be established on nursery log)
            if (species.Name == "Picejezo" || species.Name == "Picegleh") // TODO: this should be a functional type parameter
            {
                lightProbabilityTree *= nurseryLogAvailability;
                isSufficientlight = modelCore.GenerateUniform() < lightProbabilityTree;
                if (isSufficientlight) regenType = "nlog";
            }
            // Case 2. CWD-independent species (species which can be established on both forest floor & nursery log)
            else
            {
                // 1. Can the cohort establish on forest floor? (lightProbability is considering both Tree and Grass species)
                if (modelCore.GenerateUniform() < lightProbability)
                {
                    isSufficientlight = true;
                    regenType = "surface";
                }
                else
                {
                    // 2. If (1) the site shade is darker than the best shade class for the species and 
                    //       (2) the light availability above grass species layer meets the species requirement,
                    if (siteShade > bestShadeClass && modelCore.GenerateUniform() < lightProbabilityTree)
                    {
                        // 3. check if threre are sufficient amounts of downed logs?
                        isSufficientlight = modelCore.GenerateUniform() < nurseryLogAvailability;
                        if (isSufficientlight) regenType = "nlog";
                    }
                }
            }

            //if (OtherData.CalibrateMode)
            if (OtherData.CalibrateMode && !OtherData.ReduceOutputs) //Wataru, 2020.09.03
            {
                PlugIn.ModelCore.UI.WriteLine("nurseryLogPenalty:{0},{1},{2}", PlugIn.ModelCore.CurrentTime, species.Name, nurseryLogAvailability);
                PlugIn.ModelCore.UI.WriteLine("modified_lightProbability:{0},{1},{2}", PlugIn.ModelCore.CurrentTime, species.Name, lightProbability);
                PlugIn.ModelCore.UI.WriteLine("regeneration_type:{0},{1},{2}", PlugIn.ModelCore.CurrentTime, species.Name, regenType);
            }
            // ---------------------------------------------------------------------

            return (isSufficientlight, regenType); // 2020.10.31 Chhiro
            //return isSufficientlight;
            //return modelCore.GenerateUniform() < lightProbability;
        }
        // Original SufficientLight method
        //public bool SufficientLight(ISpecies species, ActiveSite site)
        //{

        //    //PlugIn.ModelCore.UI.WriteLine("  Calculating Sufficient Light from Succession.");
        //    byte siteShade = PlugIn.ModelCore.GetSiteVar<byte>("Shade")[site];

        //    double lightProbability = 0.0;
        //    bool found = false;

        //    foreach (ISufficientLight lights in sufficientLight)
        //    {

        //        //PlugIn.ModelCore.UI.WriteLine("Sufficient Light:  ShadeClass={0}, Prob0={1}.", lights.ShadeClass, lights.ProbabilityLight0);
        //        if (lights.ShadeClass == species.ShadeTolerance)
        //        {
        //            if (siteShade == 0) lightProbability = lights.ProbabilityLight0;
        //            if (siteShade == 1) lightProbability = lights.ProbabilityLight1;
        //            if (siteShade == 2) lightProbability = lights.ProbabilityLight2;
        //            if (siteShade == 3) lightProbability = lights.ProbabilityLight3;
        //            if (siteShade == 4) lightProbability = lights.ProbabilityLight4;
        //            if (siteShade == 5) lightProbability = lights.ProbabilityLight5;
        //            found = true;
        //        }
        //    }

        //    if (!found)
        //        PlugIn.ModelCore.UI.WriteLine("A Sufficient Light value was not found for {0}.", species.Name);

        //    return modelCore.GenerateUniform() < lightProbability;

        //}


        //---------------------------------------------------------------------
        /// <summary>
        /// Compute the most suitable shade class for the species
        /// This function identifies the peak of the light establishment table.
        /// </summary>
        // Chihiro 2020.01.22
        //
        private static int ComputeBestShadeClass(ISufficientLight lights)
        {
            int bestShadeClass = 0;
            double maxProbabilityLight = 0.0;
            if (lights.ProbabilityLight0 > maxProbabilityLight) bestShadeClass = 0;
            if (lights.ProbabilityLight1 > maxProbabilityLight) bestShadeClass = 1;
            if (lights.ProbabilityLight2 > maxProbabilityLight) bestShadeClass = 2;
            if (lights.ProbabilityLight3 > maxProbabilityLight) bestShadeClass = 3;
            if (lights.ProbabilityLight4 > maxProbabilityLight) bestShadeClass = 4;
            if (lights.ProbabilityLight5 > maxProbabilityLight) bestShadeClass = 5;
            return bestShadeClass;
        }


        //---------------------------------------------------------------------
        /// <summary>
        /// Compute the ratio of projected area (= occupancy area) of nursery logs to the grid area.
        /// </summary>
        // W.Hotta & Chihiro;
        //
        // Description: 
        //     - Every SiteVars.CurrentDeadWoodC[site] is downed logs.
        //     - Only the downed logs (SiteVars.CurrentDeadWoodC[site]) which decay class is between 3 to 5 
        //       are suitable for establishment and treated as nursery logs.
        //     - The carbon stocks of the nursery logs are converted to volume 
        //       using a wood density of each decay class.
        //     - Then, the volume is converted to the projected area (occupation area) 
        //       using the mean height of downed logs derived from field data.
        //         - The shape of downed logs were assumed to be an elliptical cylinder
        //
        //
        // Psudo-code:
        //
        //     double hight_of_downed_logs      // static value from field observation  // Units: cm
        //     double density_of_decay_class_XX // Wood density (g cm^-3) of dead wood for each decay class.
        //                                      // XX ranges from 3-5
        //                                      // Decay class 3-5 is suitable for establishment.
        //                                      //   Reference: Unidentified spp category in Table 3 of Ugawa et al. (2012)
        //                                      //              https://www.ffpri.affrc.go.jp/pubs/bulletin/425/documents/425-2.pdf
        //     double decayClassXXAreaRatio     // the amount of downed logs which is suitable for tree establishment
        //                                      // computed by ComputeNurseryLogC function
        //                                      // XX ranges from 3-5
        //     
        //     // Compute the area of the downed logs which decay class is XX (XX ranges from 3 to 5)
        //     double Projected_area_of_the_downed_logs
        //       = 2 * biomass_of_the_downed_logs_which_decay_class_is_XX / (Math.PI * hight_of_downed_logs * density_of_decay_class_XX)
        //     
        //     // Compute the area occupancy of downed logs which decay class is XX
        //     double decayClassXXAreaRatio = Projected_area_of_the_downed_logs / Area_of_the_grid
        //       
        //     // Compute the sum of area occupancy of decay class 3-5
        //     return Math.Min(1.0, decayClass3AreaRatio + decayClass4AreaRatio + decayClass5AreaRatio) // decayClassAreaRatio should be 0-1
        //
        //
        private static double ComputeNurseryLogAreaRatio(ISpecies species, ActiveSite site)
        {
            // Hight of downed logs
            double hight = 28.64; // Units: cm

            // Wood density (g cm^-3) of dead wood for each decay class.
            // Decay class 3-5 is suitable for establishment.
            // Reference: Unidentified spp category in Table 3 of Ugawa et al. (2012)
            //            https://www.ffpri.affrc.go.jp/pubs/bulletin/425/documents/425-2.pdf
            double densityDecayClass0 = 0.421;
            double densityDecayClass3 = 0.255;
            double densityDecayClass4 = 0.178;
            double densityDecayClass5 = 0.112;

            // Compute the amount of nursery log carbon (gC m^-2)
            double[] nurseryLogC = ComputeNurseryLogC(site, densityDecayClass0, densityDecayClass3, densityDecayClass4, densityDecayClass5);

            // Compute the area ratio in the site of the nursery log occupies.
            // The shape of downed logs were assumed to be an elliptical cylinder
            // Variables:
            //   decayClassXAreaRatio (-)
            //   nurseryLogC[X] (gC m^-2)
            //   height (cm)
            //   densityDecayClass[X] (gC cm^-3)
            double decayClass3AreaRatio = 4 * 2 * nurseryLogC[0] / (Math.PI * hight * densityDecayClass3) * Math.Pow(10, -4); // Decay class 3
            double decayClass4AreaRatio = 4 * 2 * nurseryLogC[1] / (Math.PI * hight * densityDecayClass4) * Math.Pow(10, -4); // Decay class 4
            double decayClass5AreaRatio = 4 * 2 * nurseryLogC[2] / (Math.PI * hight * densityDecayClass5) * Math.Pow(10, -4); // Decay class 5
            //if (OtherData.CalibrateMode && species.Index == 0)
            if (OtherData.CalibrateMode && species.Index == 0 && !OtherData.ReduceOutputs)
            {
                PlugIn.ModelCore.UI.WriteLine("nurseryLogC:{0},{1},{2},{3}", PlugIn.ModelCore.CurrentTime, nurseryLogC[0], nurseryLogC[1], nurseryLogC[2]);
                PlugIn.ModelCore.UI.WriteLine("decayClassAreaRatios:{0},{1},{2},{3}", PlugIn.ModelCore.CurrentTime, decayClass3AreaRatio, decayClass4AreaRatio, decayClass5AreaRatio);
            }
            return Math.Min(1.0, decayClass3AreaRatio + decayClass4AreaRatio + decayClass5AreaRatio);
        }


        //---------------------------------------------------------------------
        /// <summary>
        /// Compute the amount of nursery log carbon based on its decay ratio
        /// </summary>
        // W.Hotta & Chihiro; 
        //
        // Description: 
        //     - In the process of decomposition of downed logs, 
        //       the volume remains the same, only the density changes.
        //
        //
        // Psudo-code:
        //
        //     // Define thresholds to identify decay class
        //     double retentionRatioThreshold3 = densityDecayClass3 / densityDecayClass0
        //     double retentionRatioThreshold4 = densityDecayClass4 / densityDecayClass0
        //     double retentionRatioThreshold5 = densityDecayClass5 / densityDecayClass0
        //     
        //     // Initialize nursery log carbon for each decay class
        //     double decayClass3 = 0.0;
        //     double decayClass4 = 0.0;
        //     double decayClass5 = 0.0;
        //     
        //     // Update the amount of carbon for each decayClass
        //     for (int i = 0; i < SiteVars.CurrentDeadWoodC[site].Length; i++):
        //         // Compute the ratio of the current dead wood C to the origindal dead wood carbon (i.e. the amount of carbon just after the focused dead wood was generated.)
        //         double retentionRatio = SiteVars.CurrentDeadWoodC[site][i] / SiteVars.OriginalDeadWoodC[site][i];
        //     
        //         // Identify the decay class of the current dead wood carbon 
        //         // and update the amount of C of each decay class
        //         if (retentionRatio >= retentionRatioThreshold4) and (retentionRatio < retentionRatioThreshold3):
        //             decayClass3 += SiteVars.CurrentDeadWoodC[site][i]
        //         elif (retentionRatio >= retentionRatioThreshold5) and (retentionRatio < retentionRatioThreshold4):
        //             decayClass4 += SiteVars.CurrentDeadWoodC[site][i]
        //         elif retentionRatio < retentionRatioThreshold5:
        //             decayClass5 += SiteVars.CurrentDeadWoodC[site][i]
        //     
        //     return new double[3] { decayClass3, decayClass4, decayClass5 }
        //
        //
        private static double[] ComputeNurseryLogC(ActiveSite site, double densityDecayClass0, double densityDecayClass3, double densityDecayClass4, double densityDecayClass5)
        {
            // Define thresholds to identify decay class
            double retentionRatioThreshold3 = densityDecayClass3 / densityDecayClass0;
            double retentionRatioThreshold4 = densityDecayClass4 / densityDecayClass0;
            double retentionRatioThreshold5 = densityDecayClass5 / densityDecayClass0;

            // Initialize nursery log carbon for each decay class
            double decayClass3 = 0.0;
            double decayClass4 = 0.0;
            double decayClass5 = 0.0;

            // Update the amount of carbon for each decayClass
            for (int i = 0; i < SiteVars.CurrentDeadWoodC[site].Length; i++)
            {
                // Compute the ratio of the current dead wood C to the origindal dead wood C
                double retentionRatio = SiteVars.CurrentDeadWoodC[site][i] / SiteVars.OriginalDeadWoodC[site][i];
                // PlugIn.ModelCore.UI.WriteLine("decayRatio:{0},{1}", PlugIn.ModelCore.CurrentTime, decayRatio);

                // Identify the decay class of the current dead wood carbon & update the amount of C of each decay class (i.e. the amount of carbon just after the focused dead wood was generated.)
                if (retentionRatio >= retentionRatioThreshold4 & retentionRatio < retentionRatioThreshold3)
                {
                    decayClass3 += SiteVars.CurrentDeadWoodC[site][i];
                }
                else if (retentionRatio >= retentionRatioThreshold5 & retentionRatio < retentionRatioThreshold4)
                {
                    decayClass4 += SiteVars.CurrentDeadWoodC[site][i];
                }
                else if (retentionRatio < retentionRatioThreshold5)
                {
                    decayClass5 += SiteVars.CurrentDeadWoodC[site][i];
                }
            }
            // PlugIn.ModelCore.UI.WriteLine("decayClasses:{0},{1},{2},{3}", PlugIn.ModelCore.CurrentTime, decayClass3, decayClass4, decayClass5);
            return new double[3] { decayClass3, decayClass4, decayClass5 };
        }


        //---------------------------------------------------------------------
        /// <summary>
        /// Add a new cohort to a site following reproduction or planting.  Does not include initial communities.
        /// This is a Delegate method to base succession.
        /// </summary>

        public void AddNewCohort(ISpecies species, ActiveSite site, string reproductionType, string estLoc)
        {
            float[] initialBiomass = CohortBiomass.InitialBiomass(species, SiteVars.Cohorts[site], site);
            SiteVars.Cohorts[site].AddNewCohort(species, 1, initialBiomass[0], initialBiomass[1], estLoc);
        }
        //---------------------------------------------------------------------

        /// <summary>
        /// Determines if a species can establish on a site.
        /// This is a Delegate method to base succession.
        /// </summary>
        public bool Establish(ISpecies species, ActiveSite site)
        {
            double establishProbability = Establishment.Calculate(species, site);

            return modelCore.GenerateUniform() < establishProbability;
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Determines if a species can establish on a site.
        /// This is a Delegate method to base succession.
        /// </summary>
        public bool PlantingEstablish(ISpecies species, ActiveSite site)
        {
            IEcoregion ecoregion = modelCore.Ecoregion[site];
            double establishProbability = Establishment.Calculate(species, site);

            return establishProbability > 0.0;
        }

        //---------------------------------------------------------------------

        /// <summary>
        /// Determines if there is a mature cohort at a site.
        /// This is a Delegate method to base succession.
        /// </summary>
        public bool MaturePresent(ISpecies species, ActiveSite site)
        {
            return SiteVars.Cohorts[site].IsMaturePresent(species);
        }

        //---------------------------------------------------------------------
        // Outmoded but required?

        //public static void SiteDisturbed(object sender,
        //                                 Landis.Library.BiomassCohorts.DisturbanceEventArgs eventArgs)
        //{
        //    PlugIn.ModelCore.UI.WriteLine("  Calculating Fire or Harvest Effects.");

        //    ExtensionType disturbanceType = eventArgs.DisturbanceType;
        //    ActiveSite site = eventArgs.Site;

        //    if (disturbanceType.IsMemberOf("disturbance:fire"))
        //    {
        //        SiteVars.FireSeverity = PlugIn.ModelCore.GetSiteVar<byte>("Fire.Severity");
        //        if (SiteVars.FireSeverity != null && SiteVars.FireSeverity[site] > 0)
        //            FireEffects.ReduceLayers(SiteVars.FireSeverity[site], site);
        //    }
        //    if (disturbanceType.IsMemberOf("disturbance:harvest"))
        //    {
        //        HarvestEffects.ReduceLayers(SiteVars.HarvestPrescriptionName[site], site);
        //    }
        //}

        //---------------------------------------------------------------------

        //public static float ReduceInput(float poolInput,
        //                                  Percentage reductionPercentage,
        //                                  ActiveSite site)
        //{
        //    float reduction = (poolInput * (float)reductionPercentage);

        //    SiteVars.SourceSink[site].Carbon += (double)reduction * 0.47;

        //    return (poolInput - reduction);
        //}

        public override void InitializeSites(string initialCommunitiesText, string initialCommunitiesMap, ICore modelCore)
        {
            ModelCore.UI.WriteLine("   Loading initial communities from file \"{0}\" ...", initialCommunitiesText);
            Landis.Library.InitialCommunities.DatasetParser parser = new Landis.Library.InitialCommunities.DatasetParser(Timestep, ModelCore.Species);
            Landis.Library.InitialCommunities.IDataset communities = Landis.Data.Load<Landis.Library.InitialCommunities.IDataset>(initialCommunitiesText, parser);

            ModelCore.UI.WriteLine("   Reading initial communities map \"{0}\" ...", initialCommunitiesMap);
            IInputRaster<uintPixel> map;
            map = ModelCore.OpenRaster<uintPixel>(initialCommunitiesMap);
            using (map)
            {
                uintPixel pixel = map.BufferPixel;
                foreach (Site site in ModelCore.Landscape.AllSites)
                {
                    map.ReadBufferPixel();
                    uint mapCode = pixel.MapCode.Value;
                    if (!site.IsActive)
                        continue;

                    //if (!modelCore.Ecoregion[site].Active)
                    //    continue;

                    //modelCore.Log.WriteLine("ecoregion = {0}.", modelCore.Ecoregion[site]);

                    ActiveSite activeSite = (ActiveSite)site;
                    initialCommunity = communities.Find(mapCode);
                    if (initialCommunity == null)
                    {
                        throw new ApplicationException(string.Format("Unknown map code for initial community: {0}", mapCode));
                    }

                    InitializeSite(activeSite); //, community);
                }
            }
        }
    }
    }

