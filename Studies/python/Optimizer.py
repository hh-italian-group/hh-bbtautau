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

def TransformParams(params):
    activations = { 0: 'sigmoid', 2: 'relu', 1: 'tanh' }
    new_params = {}
    new_params['num_den_layers_pre'] = int(round(params['num_den_layers_pre']))
    new_params['num_neurons_den_layers_pre'] = int(round(params['num_neurons_den_layers_pre']))
    new_params['dropout_rate_den_layers_pre'] = params['dropout_rate_den_layers_pre']
    new_params['num_den_layers_post'] = int(round(params['num_den_layers_post']))
    new_params['num_neurons_den_layers_post'] = int(round(params['num_neurons_den_layers_post']))
    new_params['dropout_rate_den_layers_post'] = params['dropout_rate_den_layers_post']
    new_params['num_lstm_layers'] = int(round(params['num_lstm_layers']))
    new_params['num_neurons_lstm_layer'] = int(round(params['num_neurons_lstm_layer']))
    new_params['activation_dense_pre'] = activations[int(round(params['activation_dense_pre']))]
    new_params['activation_dense_post'] = activations[int(round(params['activation_dense_post']))]
    new_params['dropout_rate_lstm'] = params['dropout_rate_lstm']
    new_params['learning_rate'] = params['learning_rate']
    new_params['batch_size'] = int(round(params['batch_size']))

    return new_params

def CreateGetLoss(file_name, cfg_mean_std, cfg_min_max, n_epoch):
    data = InputsProducer.CreateRootDF(file_name, 0, True)
    X, Y,Z,var_pos, var_pos_z = InputsProducer.CreateXY(data)
    w = CreateSampleWeigts(X, Z)
    Y = Y.reshape(Y.shape[0:2])
    def GetLoss(**params):
        params = TransformParams(params)
        tf.random.set_seed(12345)
        model = pm.HHModel(var_pos, cfg_mean_std, cfg_min_max, params)
        model.call(X[0:1,:,:])
        opt = tf.keras.optimizers.Adam(learning_rate=params['learning_rate'])
        model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  weighted_metrics=[pm.sel_acc_2, pm.sel_acc_3, pm.sel_acc_4])
        model.build(X.shape)

        history = model.fit(X, Y,sample_weight=w,validation_split=0.5, epochs=n_epoch, batch_size=params['batch_size'], verbose=0)
        #model.save_weights()
        tf.keras.backend.clear_session()
        return np.amax(history.history['val_sel_acc_2'])
    return GetLoss

get_loss = CreateGetLoss(file_name, '../config/mean_std_red.json','../config/min_max_red.json', 10)

param_ranges = {   'num_den_layers_pre': (0, 5), 'num_neurons_den_layers_pre': (10,100), 'dropout_rate_den_layers_pre': (0, 0.5),
                   'num_den_layers_post': (0,5), 'num_neurons_den_layers_post': (10,100), 'dropout_rate_den_layers_post':(0, 0.5),
                   'num_lstm_layers': (1,10), 'num_neurons_lstm_layer':(10, 200), 'activation_dense_pre': (0,2),
                   'activation_dense_post': (0,2), 'dropout_rate_lstm': (0, 0.5), 'batch_size': (100,500), 'learning_rate': (0.1, 0.0001)

}

optimizer = BayesianOptimization(
    f=get_loss,
    pbounds=param_ranges,
    random_state=1, verbose=1
)

logger = JSONLogger(path=args.results)
optimizer.subscribe(Events.OPTMIZATION_STEP, logger)

optimizer.maximize(
    init_points=int(args.init_points),
    n_iter=int(args.n_iter)
)
max = optimizer.max

with open(args.max_params, 'w') as f:
    f.write(json.dumps(max, indent=4))
