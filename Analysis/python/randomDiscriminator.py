#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description='Create output with random DNN scores.')
parser.add_argument('--input', required=True, type=str, help="Input ROOT file with AnaTuple")
parser.add_argument('--output', required=True, type=str, help="Output ROOT file with DNN scores")
parser.add_argument('--tree', required=True, type=str, help="AnaTuple tree")
args = parser.parse_args()

import os
import numpy as np
import pandas
import uproot
import root_pandas

if os.path.isfile(args.output):
    os.remove(args.output)

print("Processing tree '{}'...".format(args.tree))
with uproot.open(args.input) as file:
    tree_obj = file[args.tree]
    n_total = tree_obj.numentries

df = pandas.DataFrame(data = { 'dnn_score': np.random.rand(n_total) }, dtype=np.float32)
df.to_root(args.output, mode='a', key=args.tree)

print("All entries has been processed.")
