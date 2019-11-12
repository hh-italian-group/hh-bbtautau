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

import BayesianOptimizationCustom as bo

import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("-results", "--results")
parser.add_argument("-training_variables", "--training_variables")
parser.add_argument("-init_points_to_probe", "--init_points_to_probe")
parser.add_argument("-params", "--params")
parser.add_argument("-n_iter", "--n_iter", type=int)
parser.add_argument("-n_epochs", "--n_epochs", type=int)
parser.add_argument("-load_points", "--load_points")
parser.add_argument("-f", "--file", nargs='+')
parser.add_argument('-val_split', '--val_split', nargs='?', default=0.25, type=float)
parser.add_argument('-seed', '--seed', nargs='?', default=12345, type=int)
parser.add_argument("-random_state", "--random_state",nargs='?', type=int, default=1)
parser.add_argument("-prev_point", "--prev_point", nargs='?')

args = parser.parse_args()

args = parser.parse_args()

def ListToVector(files):
    v = ROOT.std.vector('string')()
    for file in files:
        v.push_back(file)
    return v
file_name = ListToVector(args.file)

def CreateGetLoss(file_name, cfg_mean_std, cfg_min_max, n_epoch):
    np.random.seed(args.seed)
    data = InputsProducer.CreateRootDF(file_name, 0, True, True)
    X, Y, Z, var_pos, var_pos_z, var_name = InputsProducer.CreateXY(data, args.training_variables)
    w = CreateSampleWeigts(X, Z)
    Y = Y.reshape(Y.shape[0:2])
    def GetLoss(**params):
        tf.random.set_seed(args.seed)
        model = pm.HHModel(var_pos, cfg_mean_std, cfg_min_max, params)
        model.call(X[0:1,:,:])
        opt = getattr(tf.keras.optimizers, params['optimizers'])(learning_rate=10 ** params['learning_rate_exp'])

        model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  weighted_metrics=[pm.sel_acc_2])
        model.build(X.shape)

        history = model.fit(X, Y,sample_weight=w, validation_split=args.val_split, epochs=args.n_epochs,
                            batch_size=params['batch_size'], verbose=0)
        tf.keras.backend.clear_session()
        return np.amax(history.history['val_sel_acc_2'])
    return GetLoss

get_loss = CreateGetLoss(file_name, '../config/mean_std_red.json','../config/min_max_red.json', args.n_epochs)

optimizer = bo.BayesianOptimizationCustom(args.params, args.init_points_to_probe, get_loss,
                                          '{}_target.json'.format(args.results), '{}_opt.json'.format(args.results),
                                          args.n_iter, args.random_state)
#Include point from previus optimization
if args.load_points == True:
    bo.LoadPoints('target_{}'.format(args.prev_point), 'opt_{}'.format(args.prev_point))


params, result = optimizer.maximize(args.n_iter)
