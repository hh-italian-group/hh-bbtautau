# Definition of the model of the signal selection NN
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
import numpy as np
import json

from tensorflow.keras.layers import Layer, Dense, Dropout, LSTM, TimeDistributed, Concatenate
from tensorflow.keras.models import Model

class StdLayer(Layer):
    def __init__(self, file_name, var_pos, n_sigmas):
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

        super(StdLayer, self).__init__()

    def build(self, input_shape):
        super(StdLayer, self).build(input_shape)

    def call(self, X):
        Y = tf.clip_by_value(( X - self.vars_mean ) / self.vars_std, -self.n_sigmas, self.n_sigmas)
        return tf.where(self.vars_apply, Y, X)

class ScaleLayer(Layer):
    def __init__(self, file_name, var_pos, interval_to_scale):
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

        super(ScaleLayer, self).__init__()

    def build(self, input_shape):
        super(ScaleLayer, self).build(input_shape)

    def call(self, X):
        Y = tf.clip_by_value( (self.y * ( X - self.vars_min))  + self.a , self.a, self.b)
        return tf.where(self.vars_apply, Y, X)

class HHModel(Model):
    def __init__(self, var_pos, mean_std_json, min_max_json):
        super(HHModel, self).__init__()
        self.normalize = StdLayer(mean_std_json, var_pos, 5)
        self.scale = ScaleLayer(min_max_json, var_pos, [-1,1])
        self.lstm_1 = LSTM(16, return_sequences=True)
        self.dropout_1 = Dropout(0.2)
        self.concat_1 = Concatenate()
        self.lstm_2 = LSTM(16, return_sequences=True)
        self.dropout_2 = Dropout(0.2)
        self.dense_1 = TimeDistributed(Dense(10, activation="sigmoid"))
        self.dropout_3 = Dropout(0.2)
        self.dense_2 = TimeDistributed(Dense(1, activation="sigmoid"))

    def call(self, inputs):
        trans_x = self.normalize(inputs)
        trans_x = self.scale(trans_x)
        mask = tf.convert_to_tensor(inputs[:, :, 0] > 0.5, dtype=tf.bool)
        x = self.lstm_1(trans_x[:, :, 1:], mask=mask)
        x = self.dropout_1(x)
        x = self.concat_1([trans_x, x])
        x = self.lstm_2(x, mask=mask)
        x = self.dropout_2(x)
        x = self.dense_1(x, mask=mask)
        x = self.dropout_3(x)
        x = self.dense_2(x, mask=mask)
        input_shape = tf.shape(inputs)
        return tf.reshape(x, shape=(input_shape[0], input_shape[1]))

def sel_acc(y_true, y_pred, n_positions, n_exp, do_ratio):
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
        return n_valid / tf.cast(n_evt, tf.float32)
    return valid

def sel_acc_2(y_true, y_pred):
    return sel_acc(y_true, y_pred, 2, 2, False)
def sel_acc_3(y_true, y_pred):
    return sel_acc(y_true, y_pred, 3, 2, False)
def sel_acc_4(y_true, y_pred):
    return sel_acc(y_true, y_pred, 4, 2, False)
