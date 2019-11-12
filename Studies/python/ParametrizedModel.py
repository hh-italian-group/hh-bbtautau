
# Definition of the model of the signal selection NN
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
import numpy as np
import json
import ROOT

from tensorflow.keras.layers import Layer, Dense, Dropout, LSTM, GRU, SimpleRNN, TimeDistributed, Concatenate, BatchNormalization
from tensorflow.keras.models import Model

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
        self.n_sigmas = n_sigmas

        super(StdLayer, self).__init__(**kwargs)

    def build(self, input_shape):
        super(StdLayer, self).build(input_shape)

    def call(self, X):
        Y = tf.clip_by_value(( X - self.vars_mean ) / self.vars_std, -self.n_sigmas, self.n_sigmas)
        return tf.where(self.vars_apply, Y, X)

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
        return tf.where(self.vars_apply, Y, X)


class HHModel(Model):
    def __init__(self, var_pos, mean_std_json, min_max_json, params):
        super(HHModel, self).__init__()

        self.normalize = StdLayer(mean_std_json, var_pos, 5, name='std_layer')
        self.scale = ScaleLayer(min_max_json, var_pos, [-1,1], name='scale_layer')

        # Pre Dense Block
        self.dense_pre = []
        self.dropout_dense_pre = []
        self.batch_norm_dense_pre = []
        for i in range(params['num_den_layers_pre']):
            self.dense_pre.append(TimeDistributed(Dense(params['num_units_den_layers_pre'],
                                                        activation=params['activation_dense_pre']),
                                                        name='dense_pre_{}'.format(i)))
            self.batch_norm_dense_pre.append(tf.keras.layers.BatchNormalization(name='batch_normalization_pre_{}'.format(i)))
            if params['dropout_rate_den_layers_pre'] > 0:
                self.dropout_dense_pre.append(Dropout(params['dropout_rate_den_layers_pre'],
                                                      name='dropout_dense_pre_{}'.format(i)))

        # RNN Dense Block
        self.rnn = []
        self.dropout_rnn = []
        self.concatenate = []
        self.batch_norm_rnn = []
        for i in range(params['num_rnn_layers']):
            self.rnn.append(getattr(tf.keras.layers, params['rnn_type'])(params['num_units_rnn_layer'],
                             return_sequences=True ,name='rnn_{}'.format(i)))
            self.batch_norm_rnn.append(tf.keras.layers.BatchNormalization(name='batch_normalization_rnn_{}'.format(i)))
            if params['dropout_rate_rnn'] > 0:
                self.dropout_rnn.append(Dropout(params['dropout_rate_rnn'], name='dropout_rnn_pre_{}'.format(i)))
            if i < params['num_rnn_layers'] - 1 :
                self.concatenate.append(Concatenate(name='concatenate_{}'.format(i)))

        # Post Dense Block
        self.dense_post = []
        self.dropout_dense_post = []
        self.batch_norm_dense_post = []
        for i in range(params['num_den_layers_post']):
            self.dense_post.append(TimeDistributed(Dense(params['num_units_den_layers_post'],
                                                         activation=params['activation_dense_post']),
                                                          name='dense_pos_{}'.format(i)))
            self.batch_norm_dense_post.append(tf.keras.layers.BatchNormalization(name='batch_normalization_post_{}'.format(i)))
            if params['dropout_rate_den_layers_post'] > 0:
                self.dropout_dense_post.append(Dropout(params['dropout_rate_den_layers_post'],
                                                       name='dropout_dense_pre_{}'.format(i)))

        self.final_dense = TimeDistributed(Dense(1, activation="sigmoid"), name='output')

    def call(self, inputs):
        x = self.normalize(inputs)
        x = self.scale(x)
        mask = tf.convert_to_tensor(inputs[:, :, 0] > 0.5, dtype=tf.bool)

        x = x[:, :, 1:]
        for i in range(len(self.dense_pre)):
            x = self.dense_pre[i](x, mask=mask)
            if len(self.batch_norm_dense_pre) > i:
                x = self.batch_norm_dense_pre[i](x)
            if len(self.dropout_dense_pre) > i:
                x = self.dropout_dense_pre[i](x)

        last_pre = x

        for i in range(len(self.rnn)):
            x = self.rnn[i](x, mask=mask)
            if len(self.batch_norm_rnn) > i:
                x = self.batch_norm_rnn[i](x)
            if len(self.dropout_rnn) > i:
                x = self.dropout_rnn[i](x)
            if i < len(self.rnn) - 1 :
                x = self.concatenate[i]([last_pre, x])

        for i in range(len(self.dense_post)):
            x = self.dense_post[i](x, mask=mask)
            if len(self.batch_norm_dense_post) > i:
                x = self.batch_norm_dense_post[i](x)
            if len(self.dropout_dense_post) > i:
                x = self.dropout_dense_post[i](x)
        x = self.final_dense(x, mask=mask)

        input_shape = tf.shape(inputs)
        x = tf.reshape(x, shape=(input_shape[0], input_shape[1]))
        x = x * tf.cast(mask, dtype=tf.float32)
        s = tf.reshape(tf.reduce_sum(x, axis = 1), shape=(input_shape[0], 1))
        x = (2 * x ) / s
        return x

def ListToVector(files):
    v = ROOT.std.vector('string')()
    for file in files:
        v.push_back(file)
    return v

def sel_acc(y_true, y_pred, n_positions, n_exp, do_ratio, return_num=False):
    pred_sorted = tf.argsort(y_pred, axis=1, direction='DESCENDING')
    n_evt = tf.shape(y_true)[0]
    evt_id = tf.range(n_evt)
    matches_vec = []
    for n in range(n_positions):
        index = tf.transpose(tf.stack([evt_id, tf.reshape(pred_sorted[:, n], shape=(n_evt,))]))
        matches_vec.append(tf.gather_nd(y_true, index))
    matches_sum = tf.add_n(matches_vec)
    valid = tf.cast(tf.equal(matches_sum, n_exp), tf.float32)
    if do_ratio:
        n_valid = tf.reduce_sum(valid)
        ratio = n_valid / tf.cast(n_evt, tf.float32)
        if return_num:
            ratio = ratio.numpy()
        return ratio
    return valid

def sel_acc_2(y_true, y_pred):
    return sel_acc(y_true, y_pred, 2, 2, False)
def sel_acc_3(y_true, y_pred):
    return sel_acc(y_true, y_pred, 3, 2, False)
def sel_acc_4(y_true, y_pred):
    return sel_acc(y_true, y_pred, 4, 2, False)
