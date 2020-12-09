import ROOT
import argparse
import json

parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('--ch', required=False, type=str, default= "all", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "all", help= "year")
parser.add_argument('--unc_sources_group', required=False, type=str, default= "all", help="Input file name")

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

files = []
trees = []
for i in years:
    for j in channels:
        for k in unc_sources:
            filename = i+"_"+j+"_"+k+".root"
            files.append(filename)
            trees.append(j)

central_area = "/mnt/data/Dottorato/anaTuples/2020-11-30/"
r_factors={}
for file,tree in zip(files,trees):
    fileName= central_area+file
    treeName= tree
    print fileName,treeName
    d = ROOT.RDataFrame(treeName, fileName)
    not_data = '! is_data'
    d=d.Filter(not_data)
    w_before = d.Define( "num", "weight").Sum('num')
    w_after = d.Define( "den", "weight*weight_btag_IterativeFit").Sum('den')
    r_factor=w_before.GetValue()/w_after.GetValue()
    r_factors[file.strip('.root')]=r_factor
    print r_factor

with open("r_factors.json", "w") as write_file:
    json.dump(r_factors, write_file)
