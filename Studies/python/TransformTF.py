# Inclusion of the HH-btag method to the framework
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-params", "--params")
parser.add_argument("-output", "--output")
parser.add_argument("-weights", "--weights")

args = parser.parse_args()

import ParametrizedModelTF1 as pm
import json


import tensorflow as tf
from tensorflow.python.framework.graph_io import write_graph
from tensorflow.python.framework.graph_util import convert_variables_to_constants
from tensorflow.tools.graph_transforms import TransformGraph
import keras
from keras.models import load_model
from keras import backend as K

config = tf.ConfigProto(intra_op_parallelism_threads=2,
                        inter_op_parallelism_threads=2,
                        allow_soft_placement=True,
                        device_count = {'CPU' : 1, 'GPU' : 0})

session = tf.Session(config=config)
K.set_session(session)
K.set_learning_phase(0)

var_pos = {'jet_{}_valid': 0, 'jet_{}_pt': 1, 'jet_{}_eta': 2, 'rel_jet_{}_M_pt': 3, 'rel_jet_{}_E_pt': 4,
           'jet_{}_htt_deta': 5, 'jet_{}_deepFlavour': 6, 'jet_{}_htt_dphi': 7, 'sample_year': 8, 'channelId': 9,
           'htt_pt': 10, 'htt_eta': 11, 'htt_met_dphi': 12, 'rel_met_pt_htt_pt': 13, 'htt_scalar_pt': 14
           }

with open(args.params) as f:
    params = json.load(f)

n_jets = 10
model = pm.HHModel(var_pos, n_jets, '../config/mean_std_red.json', '../config/min_max_red.json', params)

print(model.summary())

model.load_weights(args.weights, by_name=True)

model.save_weights('model_weights_tf1.h5')
def node_names(nodes):
    return [ node.name.split(':')[0] for node in nodes ]

input_nodes = node_names(model.inputs)
output_nodes = node_names(model.outputs)
print(model.outputs)

with K.get_session() as sess:
    const_graph = convert_variables_to_constants(sess, sess.graph.as_graph_def(), output_nodes)
    final_graph = const_graph

out_file= '{}.pb'.format(args.output)
out_dir = '.'
write_graph(final_graph, out_dir, out_file, as_text=False)
