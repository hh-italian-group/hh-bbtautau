import ROOT
import argparse
import json
from numpy import *


class r_factor_calc:
    uncertainty_dictionary =  {"None": 0, "TauES": 1, "JetFull_Total": 2, "TopPt": 3, "JetFull_AbsoluteStat": 4, "JetFull_AbsoluteScale": 5,
        "JetFull_AbsoluteMPFBias": 6, "JetFull_AbsoluteFlavMap": 7, "JetFull_Fragmentation": 8, "JetFull_SinglePionECAL": 9,
        "JetFull_SinglePionHCAL": 10, "JetFull_FlavorQCD": 11, "JetFull_FlavorZJet": 12, "JetFull_FlavorPhotonJet": 13,
        "JetFull_FlavorPureGluon": 14, "JetFull_FlavorPureQuark": 15, "JetFull_FlavorPureCharm": 16,
        "JetFull_FlavorPureBottom": 17, "JetFull_TimePtEta": 18, "JetFull_RelativeJEREC1": 19, "JetFull_RelativeJEREC2": 20,
        "JetFull_RelativeJERHF": 21, "JetFull_RelativePtBB": 22, "JetFull_RelativePtEC1": 23, "JetFull_RelativePtEC2": 24,
        "JetFull_RelativePtHF": 25, "JetFull_RelativeBal": 26, "JetFull_RelativeFSR": 27, "JetFull_PileUpDataMC": 28,
        "JetFull_PileUpPtRef": 29, "JetFull_PileUpPtBB": 30, "JetFull_PileUpPtEC1": 31, "JetFull_PileUpPtEC2": 32,
        "JetFull_PileUpPtHF": 33, "JetFull_SubTotalPileUp": 34, "JetFull_SubTotalRelative": 35, "JetFull_SubTotalPt": 36,
        "JetFull_SubTotalScale": 37, "JetFull_SubTotalAbsolute": 38, "JetFull_SubTotalMC": 39, "JetFull_TotalNoFlavor": 40,
        "JetFull_TotalNoTime": 41, "JetFull_TotalNoFlavorNoTime": 42, "EleFakingTauES": 43, "JetReduced_Absolute": 44,
        "JetReduced_Absolute_year": 45, "JetReduced_BBEC1": 46, "JetReduced_BBEC1_year": 47, "JetReduced_EC2": 48,
        "JetReduced_EC2_year": 49, "JetReduced_FlavorQCD": 50, "JetReduced_HF": 51, "JetReduced_HF_year": 52,
        "JetReduced_RelativeBal": 53, "JetReduced_RelativeSample_year": 54, "JetReduced_Total": 55,
        "Lumi": 56, "QCDscale_W": 57, "QCDscale_WW": 58, "QCDscale_WZ": 59, "QCDscale_ZZ": 60, "QCDscale_EWK": 61,
        "QCDscale_ttbar": 62, "QCDscale_tW": 63, "QCDscale_ZH": 64, "QCDscale_ggHH": 65, "pdf_ggHH": 66, "BR_SM_H_bb": 67,
        "BR_SM_H_tautau": 68, "Eff_b": 69, "Eff_e": 70, "Eff_m": 71, "DY_0b_vLowPt": 72, "DY_0b_LowPt": 73, "DY_0b_Med1Pt": 74,
        "DY_0b_Med2Pt": 75, "DY_0b_HighPt": 76, "DY_0b_vHighPt": 77, "DY_1b_vLowPt": 78, "DY_1b_LowPt": 79, "DY_1b_Med1Pt": 80,
        "DY_1b_Med2Pt": 81, "DY_1b_HighPt": 82, "DY_1b_vHighPt": 83, "DY_2b_vLowPt": 84, "DY_2b_LowPt": 85, "DY_2b_Med1Pt": 86,
        "DY_2b_Med2Pt": 87, "DY_2b_HighPt": 88, "DY_2b_vHighPt": 89, "Qcd_norm": 90, "Qcd_sf_stat_unc": 91,
        "Qcd_sf_extrap_unc": 92, "TauTriggerUnc": 93, "EleTriggerUnc": 94, "MuonTriggerUnc": 95, "TauES_DM0": 96,
        "TauES_DM1": 97, "TauES_DM10": 98, "TauES_DM11": 99, "TauVSjetSF_DM0": 100, "TauVSjetSF_DM1": 101,
        "TauVSjetSF_3prong": 102, "TauVSjetSF_pt20to25": 103, "TauVSjetSF_pt25to30": 104, "TauVSjetSF_pt30to35": 105,
        "TauVSjetSF_pt35to40": 106, "TauVSjetSF_ptgt40": 107, "EleFakingTauES_DM0": 108, "EleFakingTauES_DM1": 109,
        "EleFakingTauES_3prong": 110, "TauVSeSF_barrel": 111, "TauVSeSF_endcap": 112, "MuFakingTauES_DM0": 113,
        "MuFakingTauES_DM1": 114, "MuFakingTauES_3prong": 115, "TauVSmuSF_etaLt0p4": 116, "TauVSmuSF_eta0p4to0p8": 117,
        "TauVSmuSF_eta0p8to1p2": 118, "TauVSmuSF_eta1p2to1p7": 119, "TauVSmuSF_etaGt1p7": 120, "MuFakingTauES": 121,
        "EleIdIsoUnc": 122, "MuonIdIsoUnc": 123, "TauTriggerUnc_DM0": 124, "TauTriggerUnc_DM1": 125, "TauTriggerUnc_DM10": 126,
        "TauTriggerUnc_DM11": 127, "L1_prefiring": 128, "PileUp": 129, "PileUpJetId_eff" : 130, "PileUpJetId_mistag" : 131,
        "TauCustomSF_DM0": 132, "TauCustomSF_DM1": 133, "TauCustomSF_DM10": 134, "TauCustomSF_DM11": 135,
        "VBFTriggerUnc_jets": 139, "VBFTauTriggerUnc_DM0": 140, "VBFTauTriggerUnc_DM1": 141, "VBFTauTriggerUnc_3prong": 142,
        "btag_lf": 143, "btag_hf": 144, "btag_hfstats1": 145, "btag_hfstats2": 146,
        "btag_lfstats1": 147, "btag_lfstats2": 148, "btag_cferr1": 149, "btag_cferr2": 150}
    variation_dictionary = { "Central" : 0, "Up" : 1, "Down" : -1 }
    unc_variations=["Up", "Down"]
    norm_unc_sources=["btag_lf","btag_hf","btag_hfstats1","btag_hfstats2","btag_lfstats1","btag_lfstats2","btag_cferr1","btag_cferr2"]
    jes_unc_sources=["JetReduced_Absolute","JetReduced_Absolute_year","JetReduced_BBEC1","JetReduced_BBEC1_year","JetReduced_EC2","JetReduced_EC2_year","JetReduced_FlavorQCD","JetReduced_HF","JetReduced_HF_year","JetReduced_RelativeBal","JetReduced_RelativeSample_year","JetReduced_Total"]
    les_unc_sources=["TauES", "TauES_DM0", "TauES_DM1", "TauES_DM10", "TauES_DM11", "EleFakingTauES_DM0", "EleFakingTauES_DM1", "MuFakingTauES"]
    unc_sources_dict={"Central": norm_unc_sources, "JES" : jes_unc_sources, "LES" : les_unc_sources}
    tune = 2
    weight_num = "weight"
    weight_den = "weight*weight_btag_IterativeFit"
    unc_sources_vector=[]
    unc_scales_vector=[]
    r_factors = []

    def __init__(self, input_dir, year, channel, unc_group, r_factors_d):
        self.input_dir = input_dir
        self.year = year
        self.channel = channel
        self.unc_group = unc_group
        self.r_factors_d = r_factors_d

    def open_dataframe(self):
        file = self.input_dir + "/"+self.year+"_"+self.channel+"_"+self.unc_group+".root"
        dataframe = ROOT.RDataFrame(self.channel, file)
        not_data = '! is_data'
        dataframe=dataframe.Filter(not_data)
        if self.tune == 1:
            dataframe=dataframe.Filter("is_TuneCP5==1")
        elif self.tune == 0:
            dataframe=dataframe.Filter("is_TuneCP5==0")
        return dataframe

    def add_dic(self):
        for (r_factor, unc_source, unc_scale) in zip (self.r_factors, self.unc_sources_vector, self.unc_scales_vector):
            print ("r factor related to %s and %s  is %.10f " % (unc_source, unc_scale,r_factor))
            if isnan(r_factor):
                print("warning, r-factor related to %s %s is nan" %(unc_source, unc_scale))
            if unc_source in self.r_factors_d:
                pass
            else:
                self.r_factors_d[unc_source]={}
            self.r_factors_d[unc_source][unc_scale]=r_factor
        return self.r_factors_d

    def evaluate_r_norm(self):
        d = self.open_dataframe()
        self.unc_sources_vector.append("None")
        self.unc_scales_vector.append("Central")
        w_before = d.Define("num", self.weight_num).Sum("num")
        w_after = d.Define("den", self.weight_den).Sum("den")
        #print(" r factor for None Central is %.10f " %(w_before.GetValue()/w_after.GetValue()))
        self.r_factors.append(w_before.GetValue()/w_after.GetValue())
        for unc_source in self.unc_sources_dict[self.unc_group]:
            for unc_scale in self.unc_variations:
                self.weight_den = self.weight_den+"*unc_"+unc_source+"_"+unc_scale
                self.unc_sources_vector.append(unc_source)
                self.unc_scales_vector.append(unc_scale)
                w_before = d.Define("num", self.weight_num).Sum("num")
                w_after = d.Define("den", self.weight_den).Sum("den")
                #print(" r factor for %s %s is %.10f " %(unc_source, unc_scale, w_before.GetValue()/w_after.GetValue()))
                self.r_factors.append(w_before.GetValue()/w_after.GetValue())
        self.add_dic()

    def evaluate_r_event(self):
        d = self.open_dataframe()
        for unc_source in self.unc_sources_dict[self.unc_group]:
            for unc_scale in self.unc_variations:
                filter_unc_source= "unc_source=="+str(self.uncertainty_dictionary[unc_source])
                filter_unc_scale= "unc_scale=="+str(self.variation_dictionary[unc_scale])
                w_before = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define("num", self.weight_num).Sum("num")
                w_after = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define("den", self.weight_den).Sum("den")
                self.r_factors.append(w_before.GetValue()/w_after.GetValue())
                self.unc_sources_vector.append(unc_source)
                self.unc_scales_vector.append(unc_scale)
                if(unc_source == "JetReduced_Total"):
                    self.weight_den = "weight*weight_btag_IterativeFit_withJES"
                    w_after = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define("den", self.weight_den).Sum("den")
                    #self.r_factors.append(w_before.GetValue()/w_after.GetValue())
                    #self.unc_sources_vector.append("JetReduced_Total_JES")
                    #self.unc_scales_vector.append(unc_scale)
                    print(" r factor with jes for %s %s is %.10f " %(unc_source, unc_scale, w_before.GetValue()/w_after.GetValue()))




    def evaluate_r(self):
        if(self.unc_group=="Central"):
            self.evaluate_r_norm()
        else:
            self.evaluate_r_event()
        self.add_dic()


parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('machine', type=str, default="pccms65", choices=['pccms65', 'local', 'gridui', 'other'])
parser.add_argument('--ch', required=False, type=str, default= "all", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "all", help= "year")
parser.add_argument('--unc_sources_group', required=False, type=str, default= "all", help="unc sources groups")
parser.add_argument('--input-dir', required=False, type=str, default="", help=" anatUples directory")
parser.add_argument('--tune', required=False, type=int, default=0, help=" 0 = use all, 1 = use isTune5 , 2= don't use tune5 ")
parser.add_argument('-n', required=False, type=bool, default=False, help=" don't write the file")


args = parser.parse_args()
channels=[]
if args.ch!="all":
    for ch in (args.ch).split(","):
      channels.append(ch)
else:
      channels.append("eTau")
      channels.append("muTau")
      channels.append("tauTau")
years=[]
if args.year!="all":
    for year in (args.year).split(","):
      years.append(year)
else:
      years.append("2016")
      years.append("2017")
      years.append("2018")

unc_sources=[]
if args.unc_sources_group!="all":
    for unc in (args.unc_sources_group).split(","):
      unc_sources.append(unc)
else:
      unc_sources.append("Central")
      unc_sources.append("LES")
      unc_sources.append("JES")


r_factors_dict={}
input_dir = ""
if args.machine == "pccms65":
    input_dir = pccms65_dir
    output = "/data/store/cms-it-hh-bbtautau//anaTuples/2020-12-01/"
elif args.machine == "local":
    #input_dir = "/mnt/data/Dottorato/anaTuples/2020-12-01/"
    input_dir = "/home/Valeria/Desktop/Dottorato/2020-12-22/"
elif args.machine == "gridui":
    input_dir = "/gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/anaTuples/2020-12-01/"
else:
    input_dir = args.input_dir

for year in years:
    for channel in channels:
        for unc_source in unc_sources:
            calculate_r_factors= r_factor_calc(input_dir, year, channel, unc_source, r_factors_dict)
            calculate_r_factors.tune=args.tune
            calculate_r_factors.evaluate_r()
        if(args.n==False):
            if(args.tune==1):
                with open("btag_correction_factors_"+channel+"_"+year+"_TuneCP5.json", "w") as write_file:
                    json.dump(r_factors_dict, write_file)
            elif(args.tune==0):
                with open("btag_correction_factors_"+channel+"_"+year+"_TuneCUETP8M1.json", "w") as write_file:
                    json.dump(r_factors_dict, write_file)
            else:
                with open("btag_correction_factors_"+channel+"_"+year+".json", "w") as write_file:
                    json.dump(r_factors_dict, write_file)
