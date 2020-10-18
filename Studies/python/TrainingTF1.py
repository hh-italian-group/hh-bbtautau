# Definition of the training for the HH-btag NN
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
from keras import Model
config = tf.ConfigProto(allow_soft_placement=False, device_count = {'CPU' : 1, 'GPU' : 1})
config.gpu_options.allow_growth = True
session = tf.Session(config=config)

from keras.callbacks import CSVLogger
from keras.callbacks import Callback
from keras.optimizers import adam
from keras import backend as K
K.set_session(session)

import argparse
import json
import numpy as np

import InputsProducer
import ParametrizedModelTF1 as pm
from CalculateWeigths import CreateSampleWeigts, CrossCheckWeights
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("-params", "--params")
parser.add_argument("-output", "--output")
parser.add_argument("-training_variables", "--training_variables")
parser.add_argument("-n_epoch", "--n_epoch", type=int)
parser.add_argument("-patience", "--patience", type=int)
parser.add_argument("-validation_split", "--validation_split", type=float)
parser.add_argument("-parity", "--parity", type=int)
parser.add_argument("-f", "--file", nargs='+')
parser.add_argument('-seed', '--seed', nargs='?', default=12345, type=int)

args = parser.parse_args()

file_name = pm.ListToVector(args.file)

with open(args.params) as f:
    params = json.load(f)

class WeightsSaver(Callback):
  def __init__(self, N):
    self.N = N
    self.epoch = 0

  def on_epoch_end(self, epoch, logs={}):
    if self.epoch % self.N == 0:
        name = '{}_par{}_weight_epoch_{}.h5'.format(args.output, args.parity, self.epoch )
        # name = ('weights_epoch_%04d.h5') % self.epoch
        self.model.save_weights(name)
    self.epoch += 1

def PerformTraining(file_name, n_epoch, params):
    np.random.seed(args.seed)
    data = InputsProducer.CreateRootDF(file_name, 0, True, True)
    X, Y, Z, var_pos, var_pos_z, var_name = InputsProducer.CreateXY(data, args.training_variables)
    print(var_pos)
    w = CreateSampleWeigts(X, Z)
    Y = Y.reshape(Y.shape[0:2])
    tf.set_random_seed(args.seed)

    model = pm.HHModel(var_pos, 10,'../config/mean_std_red.json', '../config/min_max_red.json', params)
    opt = adam(lr=10 ** params['learning_rate_exp'])
    model.compile(loss=K.binary_crossentropy,  optimizer=opt)
    model.summary()

    early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', patience=args.patience)
    csv_logger = CSVLogger('{}_par{}_training_history.csv'.format(args.output, args.parity), append=False, separator=',')
    save_best_only =  tf.keras.callbacks.ModelCheckpoint(filepath='{}_par{}_best_weights.h5'.format(args.output, args.parity),
                                                         monitor='val_loss',  mode='min', save_weights_only=True,
                                                         save_best_only=True, verbose=1)

    model.fit(X, Y, sample_weight=w, validation_split=args.validation_split, epochs=args.n_epoch, batch_size=100,
                  callbacks=[csv_logger, save_best_only,early_stop,  WeightsSaver(1)],verbose=2)

    model.save('{}_par_{}_all_model'.format(args.output,args.parity ))

    model.save_weights('{}_par{}_final_weights.h5'.format(args.output, args.parity))

with open('{}_par{}_params.json'.format(args.output, args.parity), 'w') as f:
    f.write(json.dumps(params, indent=4))

PerformTraining(file_name, args.n_epoch, params)
