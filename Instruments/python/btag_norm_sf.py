import ROOT
import argparse
import json
from numpy import *
import getDataIdHash as g
uncertainty_dictionary =  {
    "None": 0, "TauES": 1, "JetFull_Total": 2, "TopPt": 3, "JetFull_AbsoluteStat": 4, "JetFull_AbsoluteScale": 5,
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
norm_unc_sources=["None","btag_lf","btag_hf","btag_hfstats1","btag_hfstats2","btag_lfstats1","btag_lfstats2","btag_cferr1","btag_cferr2"]
jes_unc_sources=["JetReduced_Absolute","JetReduced_Absolute_year","JetReduced_BBEC1","JetReduced_BBEC1_year","JetReduced_EC2","JetReduced_EC2_year","JetReduced_FlavorQCD","JetReduced_HF","JetReduced_HF_year","JetReduced_RelativeBal","JetReduced_RelativeSample_year","JetReduced_Total"]
les_unc_sources=["TauES", "TauES_DM0", "TauES_DM1", "TauES_DM10", "TauES_DM11", "EleFakingTauES_DM0", "EleFakingTauES_DM1", "MuFakingTauES"]
# for tests
# -----------
# norm_unc_sources=["btag_lf"]
# jes_unc_sources=["JetReduced_Absolute"]
# les_unc_sources=["TauES"]
# -----------
datasets_tuneCP5= ["TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic", "ST_tW_antitop", "ST_tW_top", "ST_t-channel_antitop", "ST_t-channel_top"]
unc_sources_dict={"Central": norm_unc_sources, "JES" : jes_unc_sources, "LES" : les_unc_sources}

def add_dic(r_factor_list, dictionary):
    for r_factor in r_factor_list:
        #print ("r factor related to %s and %s with tuneCP5 %s is %.10f " % (unc_source, unc_scale, tune,r_factor))
        r_factor.value = r_factor.GetRFactor()
        if isnan(r_factor.value):
            print("warning, r-factor related to %s %s is nan" %(unc_source, unc_scale))
        if not r_factor.unc_source in dictionary:
            dictionary[r_factor.unc_source]={}
        if not r_factor.unc_scale in dictionary[r_factor.unc_source]:
            dictionary[r_factor.unc_source][r_factor.unc_scale]={}
        dictionary[r_factor.unc_source][r_factor.unc_scale][r_factor.tune] = r_factor.value
    return dictionary

class r_factor_info:
    def __init__(self):
        pass

    def GetRFactor(self):
        if self.denominator.GetValue()==0.:
            print "warning, denominator for %s %s tune %s is 0" % (self.unc_source, self.unc_scale, self.tune)
            return 0.
        else:
            return self.numerator.GetValue()/self.denominator.GetValue()

class r_factor_calc:
    unc_sources_list=[]
    unc_scales_list=[]
    tunes_list = []
    dataIds = []

    def __init__(self, input_dir, year, channel, unc_group):
        self.input_dir = input_dir
        self.year = year
        self.channel = channel
        self.unc_group = unc_group
        self.file = self.input_dir + "/"+self.year+"_"+self.channel+"_"+self.unc_group+".root"

    def open_dataframe(self):
        dataframe = ROOT.RDataFrame(self.channel, self.file)
        not_data = '! is_data'
        dataframe=dataframe.Filter(not_data)
        return dataframe

    def condition_for_tune(self, d):
        dataIds = []
        id_collections = ['dataset' ]
        data_frames = g.LoadIdFrames(self.file, id_collections)
        for col in id_collections:
            df = data_frames[col]
            for i in datasets_tuneCP5:
                dataIds.append(g.GetDataIdHash(df, i, False)[i])
        cond_str = " double tune = ("
        cmp_str = ['dataset == {}'.format(i) for i in dataIds]
        cond_str = cond_str + ' || '.join(cmp_str)
        cond_str = cond_str + ") ? 1 : 0 ; return tune; "
        return cond_str

    def filter_tune(self, d, tune):
        if self.year=="2016" and self.channel != "muTau":
            d = d.Filter("is_TuneCP5=="+tune)
        else:
            cond_str =  self.condition_for_tune(d)
            d = d.Define("isTuneCP5", cond_str)
            d = d.Filter("isTuneCP5=="+tune)
        return d

    #must return the r_factor list and stop
    def calculate_r(self, dataframe, tune, unc_source, unc_scale):
        weight_num = "weight"
        weight_den = "weight*weight_btag_IterativeFit"
        r_factor = r_factor_info()
        tune_times=0
        r_factor.tune=tune
        r_factor.unc_source = unc_source
        r_factor.unc_scale = unc_scale
        if self.unc_group != "Central" :
            filter_unc_source= "unc_source=="+str(uncertainty_dictionary[unc_source])
            filter_unc_scale= "unc_scale=="+str(variation_dictionary[unc_scale])
            dataframe = dataframe.Filter(filter_unc_source).Filter(filter_unc_scale)
        else:
            if unc_source != "None":
                weight_den = weight_den+"*unc_"+unc_source+"_"+unc_scale
        r_factor.numerator = dataframe.Define("num_"+unc_source +"_"+unc_scale, weight_num).Sum("num_"+unc_source +"_"+unc_scale )
        r_factor.denominator = dataframe.Define("den_"+unc_source +"_"+unc_scale, weight_den).Sum("den_"+unc_source +"_"+unc_scale)
        if unc_source == "JetReduced_Total":
            pass
            #implement printing of r_factor with jes
        #r_factor.value = r_factor.GetRFactor()
        return r_factor

    def evaluate_r(self, r_factor_list):
        d = self.open_dataframe()
        d_tune0 = self.filter_tune(d, "0")
        #d_tune0 = d_tune0.Range(100) # for tests
        for unc_source in unc_sources_dict[self.unc_group]:
            if unc_source == "None":
                r_factor_list.append(self.calculate_r(d_tune0, "0", unc_source, "Central"))
            else:
                for unc_scale in unc_variations:
                    r_factor_list.append(self.calculate_r(d_tune0, "0", unc_source, unc_scale))
        if self.year=="2016" and self.channel != "muTau":
            d_tune1 = self.filter_tune(d, "1")
            #d_tune1 = d_tune1.Range(100) # for tests 
            for unc_source in unc_sources_dict[self.unc_group]:
                if unc_source == "None":
                    r_factor_list.append(self.calculate_r(d_tune1, "1", unc_source, "Central"))
                else:
                    for unc_scale in unc_variations:
                        r_factor_list.append(self.calculate_r(d_tune1, "1", unc_source, unc_scale))
        return r_factor_list


parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('machine', type=str, default="pccms65", choices=['pccms65', 'local', 'gridui', 'other'])
parser.add_argument('--ch', required=False, type=str, default= "all", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "all", help= "year")
parser.add_argument('--unc_sources_group', required=False, type=str, default= "all", help="unc sources groups")
parser.add_argument('-n', required=False, type=bool, default=False, help=" don't write the file")
args = parser.parse_args()

# ************ Fill the channels, years, and unc groups **************
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

# ************* Input and output directories *********************

input_dir = ""
out_dir = ""

dir_dict = {}
dir_dict["pccms65"]={}
dir_dict["pccms65"]["input"]="/data/store/cms-it-hh-bbtautau/anaTuples/2020-12-01/"
dir_dict["pccms65"]["output"]="hh-bbtautau/McCorrections/data/"

dir_dict["local"]={}
dir_dict["local"]["input"]="/Users/valeriadamante/Desktop/Dottorato/2020-12-22/"
dir_dict["local"]["output"]="hh-bbtautau/McCorrections/data/"

dir_dict["gridui"]={}
dir_dict["gridui"]["input"]="/gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/"
dir_dict["gridui"]["output"]="hh-bbtautau/McCorrections/data/"

input_dir = dir_dict[args.machine]["input"]
input_tree= input_dir+"/${YEAR}_${CH}_{UNC_GROUP}.root"
out_dir = dir_dict[args.machine]["output"]
r_factors_dict={}
r_factors_list= []

# ********** Fill the dict and the final file *************
for year in years:
    for channel in channels:
        for unc_source in unc_sources:
            #print year, channel, unc_source
            calculate_r_factors= r_factor_calc(input_dir, year, channel, unc_source)
            calculate_r_factors.evaluate_r(r_factors_list)
            r_factors_dict = add_dic(r_factors_list, r_factors_dict)
if(args.n==False):
    with open(out_dir+"btag_correction_factors_"+channel+"_"+year+".json", "w") as write_file:
        json.dump(r_factors_dict, write_file, indent=2)
