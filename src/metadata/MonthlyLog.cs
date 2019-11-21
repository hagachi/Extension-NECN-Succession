using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Landis.Library.Metadata;

namespace Landis.Extension.Succession.NECN
{
    public class MonthlyLog
    {
        [DataFieldAttribute(Unit = FieldUnits.Year, Desc = "Simulation Year")]
        public int Time {set; get;}

        [DataFieldAttribute(Unit = FieldUnits.Month, Desc = "Simulation Month")]
        public int Month { set; get; }

        [DataFieldAttribute(Desc = "Climate Region Name")]
        public string ClimateRegionName { set; get; }

        [DataFieldAttribute(Desc = "Climate Region Index")]
        public int ClimateRegionIndex { set; get; }

        [DataFieldAttribute(Unit = FieldUnits.Count, Desc = "Number of Sites")]
        public int NumSites { set; get; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Precipitation", Format = "0.0")]
        public double ppt {get; set;}

        [DataFieldAttribute(Unit = FieldUnits.DegreeC, Desc = "Air Temperature", Format = "0.0")]
        public double airtemp { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Aboveground NPP C", Format = "0.00")]
        public double avgNPPtc { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Aboveground Heterotrophic Respiration", Format = "0.00")]
        public double avgResp { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_C_m2, Desc = "Net Ecosystem Exchange", Format = "0.00")]
        public double avgNEE { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "N Deposition and Fixation", Format = "0.00")]
        public double Ndep { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.g_N_m2, Desc = "N Leaching", Format = "0.00")]
        public double StreamN { get; set; }

        // chihiro; add additional outputs
        [DataFieldAttribute(Unit = "", Desc = "AmtNLeached", Format = "0.00")]
        public double AmtNLeached { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Water Movement", Format = "0.00")]
        public double WaterMovement { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Storm Flow", Format = "0.00")]
        public double StormFlow { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Base Flow", Format = "0.00")]
        public double BaseFlow { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "PET", Format = "0.00")]
        public double Pet { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Available Water", Format = "0.00")]
        public double AvailableWater { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Liquid Snow Pack", Format = "0.00")]
        public double LiquidSnowPack { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "Soil Water Content", Format = "0.00")]
        public double SoilWaterContent { get; set; }

        [DataFieldAttribute(Unit = "", Desc = "Decay Factor", Format = "0.00")]
        public double DecayFactor { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.DegreeC, Desc = "Soil Temperature", Format = "0.00")]
        public double SoilTemperature { get; set; }

        [DataFieldAttribute(Unit = FieldUnits.cm, Desc = "xh2o", Format = "0.00")]
        public double Xh2o { get; set; }

        [DataFieldAttribute(Unit = "", Desc = "Anaerobic Effect", Format = "0.00")]
        public double AnaerobicEffect { get; set; }
    }
}
