import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)

from tensorflow import keras
from tensorflow.keras.callbacks import CSVLogger

import argparse
import numba
import json
import numpy as np

import InputsProducer
import ParametrizedModel as pm
from CalculateWeigths import CreateSampleWeigts, CrossCheckWeights
# import Optimizer as op

parser = argparse.ArgumentParser()
parser.add_argument("-w", "--weights")
parser.add_argument("-f", "--file", nargs='+')
args = parser.parse_args()

# def CreateModel(weights):
file_name = pm.ListToVector(args.file)
with open('../config/params.json') as f:
    params = json.load(f)

data = InputsProducer.CreateRootDF(file_name, args.parity, False)
X, Y, Z, var_pos, var_pos_z, var_name = InputsProducer.CreateXY(data)

params = pm.TransformParams(params)

model = pm.HHModel(var_pos, '../config/mean_std_red.json', '../config/min_max_red.json', params)
opt = tf.keras.optimizers.Adam(learning_rate=params['learning_rate'])
model.compile(loss='binary_crossentropy',
              optimizer=opt,
              weighted_metrics=[pm.sel_acc_2, pm.sel_acc_3, pm.sel_acc_4])
model.build(X.shape)
model.load_weights(args.w, by_name=True)

@numba.njit
def Shuffle3D (X, var_idx):
    num_valid = np.count_nonzero(X[:, :, 0 ])
    X_vector  = np.zeros(num_valid)
    n = 0
    for evt in range(X.shape[0]):
        for jet in range (X.shape[1]):
            if X[evt, jet, 0 ] > 0.5:
                X_vector[n] = X[evt, jet, var_idx]
                n += 1

    X_s = np.copy(X)
    np.random.shuffle(X_vector)

    n = 0
    for evt in range(X.shape[0]):
        for jet in range (X.shape[1]):
            if X[evt, jet, 0 ] > 0.5:
                X_s[evt, jet, var_idx] = X_vector[n]
                n += 1
    return X_s

print("shape", X.shape)

ref_score = pm.sel_acc(Y, model.predict(X, batch_size=params['batch_size']), 2, 2, True, True)

def Compare(X, var_pos, ref_score, params):
    dscore = []
    dscore_var = {}
    for var_idx in range(1, len(var_pos)):
        X_shuffled = Shuffle3D(X, var_idx)
        var_score = pm.sel_acc(Y, model.predict(X_shuffled,  batch_size=params['batch_size'], steps=None), 2, 2,True, True)
        dscore_var[var_idx] = { var_name[var_idx] : ref_score - var_score}
        dscore.append(ref_score - var_score)

    # dscore.sort(reverse = True)
    return dscore, dscore_var

dscore, dscore_var =  Compare(X, var_pos, ref_score, params)
print(dscore_var)
