import ROOT
import argparse

parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('--ch', required=False, type=str, default= "eTau", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "2018", help= "year")
parser.add_argument('--input', required=False, type=str, default= " /eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-11-09/", help="Input file name")
#v_eos_directory = "/eos/home-v/vdamante/merged_trees/"
#v_eos_directory = os.getenv("merged_trees_dir")
#k_eos_directory= "/eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-11-07/"
#k_eos_directory = os.getenv("k_eos")
args = parser.parse_args()
fileName= args.input+args.year+"_"+args.ch+"_central.root"
treeName= args.ch
d = ROOT.RDataFrame(treeName, fileName)
#def GetRFactor(df):
not_data = '! is_data'
d=d.Filter(not_data)
w_before = d.Define( "num", "all_weights.at(0)").Sum('num')
w_after = d.Define( "den", "all_weights.at(0)*btag_weight_IterativeFit.at(0)").Sum('den')
r_factor=w_before.GetValue()/w_after.GetValue()
print r_factor
#return r_factor
