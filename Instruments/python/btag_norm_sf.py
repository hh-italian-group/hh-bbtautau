import ROOT
import argparse
import json
from numpy import *
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


btag_unc_sources=["btag_lf","btag_hf","btag_hfstats1","btag_hfstats2","btag_lfstats1","btag_lfstats2","btag_cferr1","btag_cferr2"]
jes_unc_sources=["JetReduced_Absolute","JetReduced_Absolute_year","JetReduced_BBEC1","JetReduced_BBEC1_year","JetReduced_EC2","JetReduced_EC2_year","JetReduced_FlavorQCD","JetReduced_HF","JetReduced_HF_year","JetReduced_RelativeBal","JetReduced_RelativeSample_year","JetReduced_Total"]
les_unc_sources=["TauES", "TauES_DM0", "TauES_DM1", "TauES_DM10", "TauES_DM11", "EleFakingTauES_DM0", "EleFakingTauES_DM1", "MuFakingTauES"]
unc_variations=['Up', 'Down']

def add_dic(r_factor, unc_source,unc_scale, r_factors_dict):
    print ("r factor related to %s and %s  is %.10f " % (unc_source, unc_scale,r_factor))
    if isnan(r_factor):
        r_factor = 1
    r_factors_dict[unc_source][unc_scale]=r_factor

def evaluate_r_Central(file,channel, unc_sources, unc_scales, r_factors_dict, isTune5):
    #print file
    d = ROOT.RDataFrame(channel, file)
    not_data = '! is_data'
    d=d.Filter(not_data)
    if isTune5 == 1:
        d=d.Filter("is_TuneCP5==1")
    elif isTune5 == 2:
        d=d.Filter("is_TuneCP5==0")
    for unc_source in unc_sources:
        for unc_scale in unc_scales:
            num = "num_"+unc_source+"_"+unc_scale
            den = "den_"+unc_source+"_"+unc_scale
            weight_num = "weight"
            if unc_source=="None":
                weight_den = "weight*weight_btag_IterativeFit"
            else:
                weight_den = "weight*weight_btag_IterativeFit*unc_"+unc_source+"_"+unc_scale
            w_central_before = d.Define(num, weight_num).Sum(num)
            w_central_after = d.Define(den, weight_den).Sum(den)
            numerator= w_central_before.GetValue()
            denominator= w_central_after.GetValue()
            #print ("numerator is %.10f \n denominator is %.10f" %(numerator, denominator))
            r_factor=numerator/denominator
            add_dic(r_factor, unc_source, unc_scale ,r_factors_dict)


def evaluate_r_JESLES(file,channel, unc_sources, unc_scales, r_factors_dict, is_withJES, isTune5):
    #print file
    d = ROOT.RDataFrame(channel, file)
    not_data = '! is_data'
    d=d.Filter(not_data)
    if isTune5 == True:
        d=d.Filter("is_TuneCP5==1")
    for unc_source in unc_sources:
        for unc_scale in unc_scales:
            filter_unc_source= "unc_source=="+str(uncertainty_dictionary[unc_source])
            filter_unc_scale= "unc_scale=="+str(variation_dictionary[unc_scale])
            num = "num_"+unc_source+"_"+unc_scale
            den = "den_"+unc_source+"_"+unc_scale
            weight_num = "weight"
            weight_den = "weight*weight_btag_IterativeFit"
            w_central_before = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define(num, weight_num).Sum(num)
            w_central_after = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define(den, weight_den).Sum(den)
            r_factor=w_central_before.GetValue()/w_central_after.GetValue()
            add_dic(r_factor, unc_source, unc_scale ,r_factors_dict)
            if is_withJES and unc_source== "JetReduced_Total":
                num = "num_"+unc_source+"_"+unc_scale+"_withJES"
                den = "den_"+unc_source+"_"+unc_scale+"_withJES"
                weight_num = "weight"
                weight_den = "weight*weight_btag_IterativeFit_withJES"
                w_central_before = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define(num, weight_num).Sum(num)
                w_central_after = d.Filter(filter_unc_source).Filter(filter_unc_scale).Define(den, weight_den).Sum(den)
                r_factor=w_central_before.GetValue()/w_central_after.GetValue()
                add_dic(r_factor, unc_source+" withJES", unc_scale ,r_factors_dict)

parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('machine', type=str, default="pccms65", choices=['pccms65', 'local', 'gridui'])
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

central_unc_sources =['None']
central_variations = ['Central']


r_factors_d={}
for i in central_unc_sources:
    for j in central_variations:
        r_factors_d[i]={j : 1.}
for i in btag_unc_sources:
    for j in unc_variations:
        r_factors_d[i]={j : 1.}

for i in jes_unc_sources:
    for j in unc_variations:
        r_factors_d[i]={j : 1.}
for i in les_unc_sources:
    for j in unc_variations:
        r_factors_d[i]={j : 1.}

for j in unc_variations:
    r_factors_d["JetReduced_Total withJES"]={j : 1.}


input_dir = ""
output = ""
if args.machine == "pccms65":
    input_dir = pccms65_dir
    output = "/data/store/cms-it-hh-bbtautau//anaTuples/2020-12-01/"
elif args.machine == "local":
    input_dir = "/mnt/data/Dottorato/anaTuples/2020-12-01/"
elif args.machine == "gridui":
    input_dir = "/gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/anaTuples/2020-12-01/"
else:
    print "unknown input dir"

for year in years:
    for channel in channels:
        file = input_dir+year+"_"+channel+"_Central.root"
        evaluate_r_Central(file,channel,central_unc_sources, central_variations, r_factors_d, args.tune)
        evaluate_r_Central(file,channel,btag_unc_sources, unc_variations, r_factors_d, args.tune)
        file = input_dir+year+"_"+channel+"_JES.root"
        evaluate_r_JESLES(file,channel,jes_unc_sources, unc_variations, r_factors_d, True, args.tune)
        file = input_dir+year+"_"+channel+"_LES.root"
        evaluate_r_JESLES(file,channel,les_unc_sources, unc_variations, r_factors_d, False, args.tune)
        if(args.n==False):
            with open("btag_correction_factors_"+channel+"_"+year+".json", "w") as write_file:
                json.dump(r_factors_d, write_file)
