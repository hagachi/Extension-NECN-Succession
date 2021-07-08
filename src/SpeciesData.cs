//  Author: Robert Scheller, Melissa Lucash

using Landis.Core;
using Landis.SpatialModeling;
using Landis.Utilities;
using Landis.Library.Succession;
using Landis.Library.Climate;
using Landis.Library.Parameters;

using System.Collections.Generic;
using System;
using System.IO;

namespace Landis.Extension.Succession.NECN
{


    public class SpeciesData 
    {
        public static Landis.Library.Parameters.Species.AuxParm<int> FuncType;
        public static Landis.Library.Parameters.Species.AuxParm<bool> NFixer;
        public static Landis.Library.Parameters.Species.AuxParm<int> GDDmin;
        public static Landis.Library.Parameters.Species.AuxParm<int> GDDmax;
        public static Landis.Library.Parameters.Species.AuxParm<int> MinJanTemp;
        public static Landis.Library.Parameters.Species.AuxParm<double> MaxDrought;
        public static Landis.Library.Parameters.Species.AuxParm<double> LeafLongevity;
        public static Landis.Library.Parameters.Species.AuxParm<bool> Epicormic;
        public static Landis.Library.Parameters.Species.AuxParm<double> LeafLignin;
        public static Landis.Library.Parameters.Species.AuxParm<double> WoodLignin;
        public static Landis.Library.Parameters.Species.AuxParm<double> CoarseRootLignin;
        public static Landis.Library.Parameters.Species.AuxParm<double> FineRootLignin;
        public static Landis.Library.Parameters.Species.AuxParm<double> LeafCN;
        public static Landis.Library.Parameters.Species.AuxParm<double> WoodCN;
        public static Landis.Library.Parameters.Species.AuxParm<double> CoarseRootCN;
        public static Landis.Library.Parameters.Species.AuxParm<double> LeafLitterCN;
        public static Landis.Library.Parameters.Species.AuxParm<double> FineRootCN;
        public static Landis.Library.Parameters.Species.AuxParm<int> Max_ANPP;
        public static Landis.Library.Parameters.Species.AuxParm<int> Max_Biomass;
        public static Landis.Library.Parameters.Species.AuxParm<double> KCompLimit;

        //---------------------------------------------------------------------
        public static void Initialize(IInputParameters parameters)
        {
            FuncType            = parameters.SppFunctionalType;
            NFixer              = parameters.NFixer;
            GDDmin              = parameters.GDDmin;
            GDDmax              = parameters.GDDmax;
            MinJanTemp          = parameters.MinJanTemp;
            MaxDrought          = parameters.MaxDrought;
            LeafLongevity       = parameters.LeafLongevity;
            Epicormic           = parameters.Epicormic;
            LeafLignin          = parameters.LeafLignin;
            WoodLignin          = parameters.WoodLignin ;
            CoarseRootLignin    = parameters.CoarseRootLignin ;
            FineRootLignin      = parameters.FineRootLignin ;
            LeafCN              = parameters.LeafCN;
            WoodCN              = parameters.WoodCN;
            CoarseRootCN        = parameters.CoarseRootCN;
            LeafLitterCN        = parameters.FoliageLitterCN;
            FineRootCN          = parameters.FineRootCN;
            Max_ANPP            = parameters.MaxANPP;
            Max_Biomass         = parameters.MaxBiomass;
            KCompLimit          = parameters.KCompLimit;

            foreach (ISpecies spp in PlugIn.ModelCore.Species)
            {
                try
                {
                    double maxLAI = FunctionalType.Table[SpeciesData.FuncType[spp]].MAXLAI;
                    //PlugIn.ModelCore.UI.WriteLine("Spp={0}, FT={1}", spp.Name, SpeciesData.FuncType[spp]);

                }
                catch (Exception)
                {
                    string mesg = string.Format("Species or Functional Type Missing: {0}", spp.Name);
                    throw new System.ApplicationException(mesg);
                }
            }

        }
    }
}
