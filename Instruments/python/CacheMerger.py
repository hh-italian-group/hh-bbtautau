import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-inputs", "--inputs", nargs='+')
parser.add_argument("-output", "--output")

args = parser.parse_args()

files = ''
for file in args.inputs:
    files  += file + ' '

cmd_line = 'hadd -f {} {}'.format(args.output, files)
print(cmd_line)
os.system(cmd_line)
