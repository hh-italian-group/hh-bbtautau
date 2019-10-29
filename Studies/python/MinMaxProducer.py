# Calculates min and max values.
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import argparse
import numpy as np
import math
import json
import sys
sys.path.insert(0, "../include")
import InputsProducer
import ROOT


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", nargs='+')
parser.add_argument("-output", "--output_file")
args = parser.parse_args()

def ListToVector(files):
    v = ROOT.std.vector('string')()
    for file in files:
        v.push_back(file)
    return v
file_name = ListToVector(args.file)

def truncate(x):
    return int(x * 1000) / 1000.

def CalculateMinMax(file_name):
    data = InputsProducer.CreateRootDF(file_name, -1, False)
    X, Y,training_vars,var_pos,var_pos_z, var_name  = InputsProducer.CreateXY(data)

    val = var_pos['jet_{}_valid']
    n_valid = np.sum(X[:,:,val])
    min_max = {}
    # for var,pos in var_pos.items():
    min = np.amin(X[:, :, 10])
    max = np.amax(X[:, :, 10])
    min_max[1] = { 'min': truncate(float(min)), 'max': truncate(float(max)) }

    with open(args.output_file, 'w') as f:
        f.write(json.dumps(min_max, indent=4))

CalculateMinMax(file_name)
