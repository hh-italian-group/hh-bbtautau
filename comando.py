import numpy as np
from array import array
import sys
import os 
from datetime import date

# CH=eTau; for YEAR in 2018; do echo Processing $YEAR $CH; ./run.sh ProcessAnaTuple --sources hh-bbtautau/Analysis/config/$YEAR/sources.cfg --setup full -
#-channel $CH --period Run$YEAR --input /eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-11-07/${YEAR}_${CH}_v3.root --output output/shapes/11_11_2020/new_bweight/bjet_vars/${YEAR}_$CH --vars b1_pt,b2_pt,nbtag_Medium --draw 1 --n_threads 4 --input_friends /eos/home-v/vdamante/merged_trees/${YEAR}_${CH}_hh_btag.root --shapes 0; done

import argparse
parser = argparse.ArgumentParser(description='Create bash command')
parser.add_argument('--ch', required=False, type=str, default= "eTau", help= "channel")
parser.add_argument('--year', required=False, type=str, default= "2018", help= "year")
parser.add_argument('--input', required=False, type=str, default= "${YEAR}_${CH}_v3.root", help="Input file name")
parser.add_argument('--input_friends', required=False, type=str, default= "${YEAR}_${CH}_hh_btag.root", help="Output directory name")
parser.add_argument('--vars', required=False, type=str, default = "dnn_score+,", help="input variables separated by comma")
parser.add_argument('--btagdir', required=False, type=str, default = "/new_bweight/",)
#parser.add_argument('--output', required=True, type=str, default= "/new_bweight/dnn_vars/", help="Output directory name")
args = parser.parse_args()

#v_eos_directory = "/eos/home-v/vdamante/merged_trees/"
v_eos_directory = os.getenv("merged_trees_dir")
#k_eos_directory= "/eos/home-k/kandroso/cms-it-hh-bbtautau/anaTuples/2020-11-07/"
k_eos_directory = os.getenv("k_eos")
out_dir = "../output/shapes/"
today = date.today()
today_string = today.strftime("%d_%m_%Y")
out_directories = {
"b1_pt,b2_pt,nbtag_Medium":"bjet_vars",
"b1_DeepFlavour,b2_DeepFlavour,b1_HHbtag,b2_HHbtag":"btag_vars",
"dnn_score+,":"dnn_vars",
"SVfit_m":"svfit_vars",
"kinFit_convergence":"kinfit_vars",
"tau1_eta,tau1_pt,tau2_eta,tau2_pt,mt_1,mt_2":"tau_vars",
"dR_l1l2Pt_htautau":"tautau_vars",
}

output = out_dir+today_string+args.btagdir+out_directories[args.vars]+"/${YEAR}_$CH"
print output
#if(!os.path.exists(out_dir_final)):
#    os.mkdir(out_dir_final)

commandBase_1 = "; do echo Processing $YEAR $CH; ./run.sh ProcessAnaTuple --sources hh-bbtautau/Analysis/config/$YEAR/sources.cfg --setup full --channel $CH --period Run$YEAR"
commandBase_2 = " --draw 1 --n_threads 4 --input_friends "
commandBase_3 = " --shapes 0; done"

commandBase_1 = "CH="+args.ch +"; for YEAR in "+args.year+commandBase_1
input = k_eos_directory+args.input
input_friends = v_eos_directory+args.input_friends
command = commandBase_1 + " --input "+ input + " --output " + output + " --vars " + args.vars + commandBase_2 + input_friends + commandBase_3

print command
