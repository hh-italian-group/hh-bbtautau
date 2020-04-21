# Script to calculate the predictions using the model in TF1
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import tensorflow as tf
import numpy as np
import argparse
import time
import InputsProducer

parser = argparse.ArgumentParser()
parser.add_argument("-model", "--model")
parser.add_argument("-input", "--x")
parser.add_argument("-output", "--output")

args = parser.parse_args()

model_file= args.model

X = np.load(args.x)

with tf.gfile.GFile(model_file, 'rb') as f:
    graph_def = tf.GraphDef()
    graph_def.ParseFromString(f.read())


with tf.Graph().as_default() as graph:
    tf.import_graph_def(graph_def, name="HHModel")
    config = tf.ConfigProto(intra_op_parallelism_threads=2,
                            inter_op_parallelism_threads=2,
                            allow_soft_placement=True,
                            device_count = {'CPU' : 1, 'GPU' : 0})
    sess = tf.Session(graph=graph, config=config)

    #Debug part
    # for n in graph.as_graph_def().node:
    #    print(n.name, n)
    # layer_names = [n.name for n in graph.as_graph_def().node]
    # print(layer_names)
    # for op in graph.get_operations():
    #     if op.name == 'HHModel/batch_normalization_rnn_0/batchnorm_1/mul_1':
    #         print(op)
    #         print("n_inputs={}".format(len(op.inputs)))
    #         for x in op.inputs:
    #             print(x)
    # raise RuntimeError("stop")
    x_graph = graph.get_tensor_by_name('HHModel/input:0')
    y_graph = graph.get_tensor_by_name('HHModel/NormToTwo/div:0')

    sess.run(tf.global_variables_initializer())

    start = time.time()
    pred = sess.run(y_graph, feed_dict={ x_graph: X})
    end = time.time()
    print('Passed time:',end - start)

    np.save(args.output, pred)
