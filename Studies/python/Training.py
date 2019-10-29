import tensorflow as tf
gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)

from tensorflow import keras
from tensorflow.keras.callbacks import CSVLogger

import argparse
import json

import InputsProducer
import ParametrizedModel as pm
from CalculateWeigths import CreateSampleWeigts, CrossCheckWeights
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("-params", "--params")
parser.add_argument("-output", "--output")
parser.add_argument("-n_epoch", "--n_epoch", type=int)
parser.add_argument("-parity", "--parity", type=int)
parser.add_argument("-f", "--file", nargs='+')

args = parser.parse_args()

file_name = pm.ListToVector(args.file)

with open(args.params) as f:
    params = json.load(f)

def PerformTraining(file_name, n_epoch, params):
    print(params)
    data = InputsProducer.CreateRootDF(file_name, 0, True)
    X, Y, Z, var_pos, var_pos_z, var_name = InputsProducer.CreateXY(data)
    w = CreateSampleWeigts(X, Z)
    Y = Y.reshape(Y.shape[0:2])
    tf.random.set_seed(12345)

    params = pm.TransformParams(params)

    model = pm.HHModel(var_pos, '../config/mean_std_red.json', '../config/min_max_red.json', params)
    model.call(X[0:1,:,:])
    opt = tf.keras.optimizers.Adam(learning_rate=params['learning_rate'])
    model.compile(loss='binary_crossentropy',
                  optimizer=opt,
                  weighted_metrics=[pm.sel_acc_2, pm.sel_acc_3, pm.sel_acc_4])
    model.build(X.shape)

    model.summary()

    early_stop = tf.keras.callbacks.EarlyStopping(monitor='val_sel_acc_2', mode='max', patience=10)
    csv_logger = CSVLogger('{}_par{}_training_history.csv'.format(args.output, args.parity), append=False, separator=',')
    save_best_only =  tf.keras.callbacks.ModelCheckpoint(filepath='{}_par{}_best_weights.h5'.format(args.output, args.parity), monitor='val_sel_acc_2',
                                                         mode='max', save_best_only=True, verbose=1)

    model.fit(X, Y, sample_weight=w, validation_split=0.25, epochs=n_epoch, batch_size=params['batch_size'],
              callbacks=[early_stop, csv_logger, save_best_only],verbose=2)

    model.save_weights('{}_par{}_final_weights.h5'.format(args.output, args.parity))

with open('{}_par{}_params.json'.format(args.output, args.parity), 'w') as f:
    f.write(json.dumps(params, indent=4))

PerformTraining(file_name, args.n_epoch, params)
