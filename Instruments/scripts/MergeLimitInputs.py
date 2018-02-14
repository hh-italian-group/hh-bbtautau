#!/usr/bin/env python
# Merge files with shapes for limits.
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import os
import re
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Merge files with shapes for limits.')
parser.add_argument('--input', required=True, dest='input', type=str, metavar='PATH', help="input path")
parser.add_argument('--output', required=True, dest='output', type=str, metavar='PATH', help="output path")
parser.add_argument('--channels', required=False, dest='channels', type=str, default="eTau,muTau,tauTau",
                    metavar='LIST', help="channel list")
parser.add_argument('--var', required=True, dest='var', type=str, metavar='VAR',
                    help="variable used for signal extraction")
args = parser.parse_args()

channels = re.split(',', args.channels)
file_patterns = set()

for ch in channels:
    pattern = '^{}_{}_.*\.root$'.format(ch, args.var)
    ch_files = [f for f in os.listdir(args.input) if re.match(pattern, f)]
    ch_patterns = [ f[len(ch):] for f in ch_files ]
    file_patterns.update(ch_patterns)

if not os.path.exists(args.output):
    os.makedirs(args.output)

for pattern in file_patterns:
    cmd_line = 'hadd -f9 {}/all{}'.format(args.output, pattern)

    for ch in channels:
        f = '{}/{}{}'.format(args.input, ch, pattern)
        if os.path.isfile(f):
            cmd_line += ' {}'.format(f)

    print('% {}'.format(cmd_line))
    subprocess.call([cmd_line], shell=True)
