#!/usr/bin/env python
# Calculates mean values and standard deviation.
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import numpy as np
import math
import json
import sys
sys.path.insert(0, "../include")
import InputsProducer

def truncate(x):
    return int(x * 1000) / 1000.

def CalculateMeanStd(file_name, output_file):
    data = InputsProducer.CreateRootDF(file_name, -1, False)
    X, Y,training_vars, var_pos, var_pos_z = InputsProducer.CreateXY(data)

    val = var_pos['jet_{}_valid']
    n_valid = np.sum(X[:,:,val])
    mean_std = {}
    for var,pos in var_pos.items():
        mean = np.sum(X[:, :, pos]) / n_valid
        std = math.sqrt(np.sum(X[:, :, val] * ((mean - X[:, :, pos]) ** 2)) / (n_valid - 1))
        mean_std[var] = { 'mean': truncate(float(mean)), 'std': truncate(float(std)) }

    with open(output_file, 'w') as f:
        f.write(json.dumps(mean_std, indent=4))
