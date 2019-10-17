import numpy as np
import math
import json
import sys
sys.path.insert(0, "../include")
import InputsProducer


def truncate(x):
    return int(x * 1000) / 1000.

def CalculateMinMax(file_name, output_file):
    data = InputsProducer.CreateRootDF(file_name, -1, False)
    X, Y, var_pos,training_vars = InputsProducer.CreateXY(data)

    val = var_pos['jet_{}_valid']
    n_valid = np.sum(X[:,:,val])
    min_max = {}
    for var,pos in var_pos.items():
        min = np.amin(X[:, :, pos])
        max = np.amax(X[:, :, pos])
        min_max[var] = { 'min': truncate(float(min)), 'max': truncate(float(max)) }

    f = open(output_file, 'w')
    f.write(json.dumps(min_max, indent=4))
    f.close()
