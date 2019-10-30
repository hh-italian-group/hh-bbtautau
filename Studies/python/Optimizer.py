# Definition of the model optimization using Bayesian optimization as a method
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)

from tensorflow import keras
from tensorflow.keras.callbacks import CSVLogger

from bayes_opt import BayesianOptimization
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events

import sys
import argparse
import pandas as pd
import numpy as np
import json
import ParametrizedModel as pm
import InputsProducer
from CalculateWeigths import CreateSampleWeigts

import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("-results", "--results")
parser.add_argument("-max", "--max_params")
parser.add_argument("-n_iter", "--n_iter", type=int)
parser.add_argument("-init_points", "--init_points", type=int)
parser.add_argument("-f", "--file", nargs='+')
args = parser.parse_args()

def ListToVector(files):
    v = ROOT.std.vector('string')()
    for file in files:
        v.push_back(file)
    return v
file_name = ListToVector(args.file)

def CreateGetLoss(file_name, cfg_mean_std, cfg_min_max, n_epoch):
    data = InputsProducer.CreateRootDF(file_name, 0, True)
    X, Y, Z, var_pos, var_pos_z, var_name = InputsProducer.CreateXY(data)
    w = CreateSampleWeigts(X, Z)
    Y = Y.reshape(Y.shape[0:2])
    def GetLoss(**params):
        params = pm.TransformParams(params)
        tf.random.set_seed(12345)
        model = pm.HHModel(var_pos, cfg_mean_std, cfg_min_max, params)
        model.call(X[0:1,:,:])
        if params['optimizers'] == 'Nadam':
            opt = tf.keras.optimizers.Nadam(learning_rate=10 ** params['learning_rate_exp'])
        elif params['optimizers'] == 'Adam':
            opt = tf.keras.optimizers.Adam(learning_rate=10 ** params['learning_rate_exp'])
        elif params['optimizers'] == 'SGD':
            opt = tf.keras.optimizers.SGD(learning_rate=10 ** params['learning_rate_exp'])
        elif params['optimizers'] == 'RMSprop':
            opt = tf.keras.optimizers.RMSprop(learning_rate=10 ** params['learning_rate_exp'])

        model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  weighted_metrics=[pm.sel_acc_2, pm.sel_acc_3, pm.sel_acc_4])
        model.build(X.shape)

        history = model.fit(X, Y,sample_weight=w,validation_split=0.5, epochs=n_epoch, batch_size=params['batch_size'], verbose=0)
        #model.save_weights()
        tf.keras.backend.clear_session()
        return np.amax(history.history['val_sel_acc_2'])
    return GetLoss

get_loss = CreateGetLoss(file_name, '../config/mean_std_red.json','../config/min_max_red.json', 20)

param_ranges = {'num_den_layers_pre': (0, 3), 'num_neurons_den_layers_pre': (10,50), 'dropout_rate_den_layers_pre': (0, 0.5),
                'num_den_layers_post': (0,3), 'num_neurons_den_layers_post': (10,50), 'dropout_rate_den_layers_post':(0, 0.5),
                'num_lstm_layers': (3,7), 'num_neurons_lstm_layer':(50, 200), 'activation_dense_pre': (0,2),
                'activation_dense_post': (0,2), 'dropout_rate_lstm': (0, 0.5), 'batch_size': (100,200),
                'learning_rate_exp': (-5, -1), 'lstm_activation': (0,2), 'lstm_recurrent_activation': (0,2),
                'optimizers': (0, 3)

}

optimizer = BayesianOptimization(
    f=get_loss,
    pbounds=param_ranges,
    random_state=1, verbose=1
)

last_params={"num_den_layers_pre" : 0 , "num_neurons_den_layers_pre": 81, "dropout_rate_den_layers_pre":  0.27,
            "num_den_layers_post": 0, "num_neurons_den_layers_post": 44, "dropout_rate_den_layers_post": 0.056,
            "num_lstm_layers": 5, "num_neurons_lstm_layer":74, "activation_dense_pre": 0, "activation_dense_post": 1,
            "dropout_rate_lstm": 0.15, "batch_size":114, "learning_rate_exp":-3, "lstm_activation": 1,
            "lstm_recurrent_activation": 0, "optimizers": 2}

optimizer.probe(params=last_params, lazy=True)

# optimizer.probe(params=[ 0 , 81, 0.27, 0, 44, 0.056, 5, 74, 0, 1, 0.15, 114, -3, 1, 0, 2 ], lazy=True)


logger = JSONLogger(path=args.results)
optimizer.subscribe(Events.OPTMIZATION_STEP, logger)

optimizer.maximize(
    init_points=int(args.init_points),
    n_iter=int(args.n_iter)
)
max = optimizer.max

with open(args.max_params, 'w') as f:
    f.write(json.dumps(max, indent=4))
