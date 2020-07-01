# Definition of the model of the signal selection NN using TF1
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
import numpy as np
import json
import ROOT

from keras.layers import Input, Layer, Lambda, Dense, Dropout, LSTM, GRU, SimpleRNN, TimeDistributed, Concatenate, BatchNormalization, Embedding
from keras import Model

class StdLayer(Layer):
    def __init__(self, file_name, var_pos, n_sigmas, **kwargs):
        with open(file_name) as json_file:
            data_json = json.load(json_file)
        n_vars = len(var_pos)
        self.vars_std = [1] * n_vars
        self.vars_mean = [1] * n_vars
        self.vars_apply = [False] * n_vars
        for var, ms in data_json.items():
            pos = var_pos[var]
            self.vars_mean[pos] = ms['mean']
            self.vars_std[pos] = ms['std']
            self.vars_apply[pos] = True
        self.vars_mean = tf.constant(self.vars_mean, dtype=tf.float32)
        self.vars_std = tf.constant(self.vars_std, dtype=tf.float32)
        self.vars_apply = tf.constant(self.vars_apply, dtype=tf.bool)
        self.n_sigmas = tf.constant(n_sigmas, dtype=tf.float32)
        super(StdLayer, self).__init__(**kwargs)

    def call(self, X):
        Y = tf.clip_by_value(( X - self.vars_mean ) / self.vars_std, -self.n_sigmas, self.n_sigmas)
        vars_apply = tf.logical_and(tf.ones_like(X, dtype=tf.bool), self.vars_apply)
        return tf.where(vars_apply, Y, X)

class ScaleLayer(Layer):
    def __init__(self, file_name, var_pos, interval_to_scale, **kwargs):
        with open(file_name) as json_file:
            data_json = json.load(json_file)
        self.a = interval_to_scale[0]
        self.b = interval_to_scale[1]
        n_vars = len(var_pos)
        self.vars_max = [1] * n_vars
        self.vars_min = [1] * n_vars
        self.vars_apply = [False] * n_vars
        for var, mm in data_json.items():
            pos = var_pos[var]
            self.vars_min[pos] = mm['min']
            self.vars_max[pos] = mm['max']
            self.vars_apply[pos] = True
        self.vars_min = tf.constant(self.vars_min, dtype=tf.float32)
        self.vars_max = tf.constant(self.vars_max, dtype=tf.float32)
        self.vars_apply = tf.constant(self.vars_apply, dtype=tf.bool)
        self.y = (self.b - self.a) / (self.vars_max - self.vars_min)

        super(ScaleLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        super(ScaleLayer, self).build(input_shape)

    def call(self, X):
        Y = tf.clip_by_value( (self.y * ( X - self.vars_min))  + self.a , self.a, self.b)
        vars_apply = tf.logical_and(tf.ones_like(X, dtype=tf.bool), self.vars_apply)
        return tf.where(vars_apply, Y, X)

class NormToTwo(Layer):
    def __init__(self, **kwargs):
        self.supports_masking = True
        super(NormToTwo, self).__init__(**kwargs)

    def compute_mask(self, inputs, mask):
        return mask

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1])

    def call(self, x, mask=None):
        if mask is None:
            raise RuntimeError("Mask is none")
        input_shape = tf.shape(x)
        x = tf.reshape(x, shape=(input_shape[0], input_shape[1]))
        x = x * tf.cast(mask, dtype=tf.float32)
        s = tf.reshape(tf.reduce_sum(x, axis = 1), shape=(input_shape[0], 1))
        x =  x = (2 * x ) / s
        return x

class Slice(Layer):
    def __init__(self, **kwargs):
        self.supports_masking = True
        super(Slice, self).__init__(**kwargs)

    def compute_mask(self, inputs, mask):
        return mask

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], input_shape[2] - 1)

    def call(self, x, mask=None):
        if mask is None:
            raise RuntimeError("Slice Mask is none")
        x = x[:,:,1:]
        return x


class Masking(Layer):
    def __init__(self, **kwargs):
        super(Masking, self).__init__(**kwargs)
        self.supports_masking = True

    def compute_mask(self, inputs, mask=None):
        return inputs[:, :, 0] > 0.5

    def call(self, inputs):
        boolean_mask = inputs[:, :, 0] > 0.5
        input_shape = tf.shape(inputs)
        boolean_mask = tf.reshape(boolean_mask, [input_shape[0], input_shape[1], 1])
        return inputs * tf.cast(boolean_mask, tf.float32)

    def compute_output_shape(self, input_shape):
        return input_shape

def HHModel(var_pos, n_jets, mean_std_json, min_max_json, params):
    n_vars = len(var_pos)
    input = Input(shape=(n_jets, n_vars), name="input")
    normalize = StdLayer(mean_std_json, var_pos, 5, name='std_layer')(input)
    scale = ScaleLayer(min_max_json, var_pos, [-1,1], name='scale_layer')(normalize)

    masked_scale = Masking(name='masking')(scale)
    x = Slice(name='slice')(masked_scale)

    # Dense Pre Block
    # for i in range(params['num_den_layers_pre']):
    #     x = masked_scale if i == 0 else dense
    #     dense = TimeDistributed(Dense(params['num_units_den_layers_pre'], activation=params['activation_dense_pre']), name='dense_pre_{}'.format(i)) (x)
    #     BatchNormalization(name='batch_normalization_pre_{}'.format(i))(x)
    last_pre = x
    # #  RNN Block
    for i in range(params['num_rnn_layers']):
        x = LSTM(params['num_units_rnn_layer'], implementation=2, recurrent_activation='sigmoid',return_sequences=True,name='rnn_{}'.format(i)) (x)
        x = BatchNormalization(name='batch_normalization_rnn_{}'.format(i))(x)
        if params['dropout_rate_rnn'] > 0.0:
            x = Dropout(params['dropout_rate_rnn'], name='dropout_rnn_pre_{}'.format(i))(x)
        if i < params['num_rnn_layers'] - 1 :
            x = Concatenate(name='concatenate_{}'.format(i))([last_pre, x])

    # Dense Post Block
    for i in range(params['num_den_layers_post']):
        x = TimeDistributed(Dense(params['num_units_den_layers_post'], activation=params['activation_dense_post']),
                            name='dense_pos_{}'.format(i)) (x)
        x = BatchNormalization(name='batch_normalization_post_{}'.format(i))(x)
        if params['dropout_rate_rnn'] > 0.0:
            x = Dropout(params['dropout_rate_rnn'], name='dropout_rnn_pre_{}'.format(i))(x)

    output = TimeDistributed(Dense(1, activation="sigmoid"), name='output') (x)
    norm = NormToTwo(name='NormToTwo')(output)

    return Model([input], [norm], "HHModel")

def ListToVector(files):
    v = ROOT.std.vector('string')()
    for file in files:
        v.push_back(file)
    return v

def sel_acc(y_true, y_pred, n_positions, n_exp, do_ratio, return_num=False):
    n = 10
    pred_sorted = np.argsort(y_pred, axis=1)
    pred_sorted = np.fliplr(pred_sorted)
    n_evt = np.shape(y_true)[0]
    evt_id = np.arange(n_evt)
    matches_vec = []
    for n in range(n_positions):
        index = pred_sorted[:, n]
        matches_vec.append(y_true[evt_id, index])
    matches_sum = np.sum(matches_vec, axis=0)
    valid = np.equal(matches_sum, n_exp)
    if do_ratio:
        n_valid = np.sum(valid)
        ratio = n_valid / n_evt
        return ratio
    return valid

def sel_acc_2(y_true, y_pred):
    return sel_acc(y_true, y_pred, 2, 2, False)
def sel_acc_3(y_true, y_pred):
    return sel_acc(y_true, y_pred, 3, 2, False)
def sel_acc_4(y_true, y_pred):
    return sel_acc(y_true, y_pred, 4, 2, False)
