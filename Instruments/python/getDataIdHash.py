#!/usr/bin/env python
# Get dataId hash
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import uproot
import pandas
import re

def LoadDataIdFrame(file):
    """Load pandas DataFrame from root file."""
    file = uproot.open(file)
    aux = file['aux']
    df = aux.arrays('*', outputtype=pandas.DataFrame)

    names = [ name.decode("utf-8") for name in df.dataId_names[0] ]
    new_df = pandas.DataFrame(data = {
                              'id' : df.dataIds[0],
                              'name': names
                              })
    return new_df

def GetDataIdHash(df, dataId, use_regex=False):
    """Retrieve hash for the given dataId."""
    id_df = df[df.name.str.contains(dataId, regex=use_regex)]
    result = {}
    for n in range(id_df.shape[0]):
        result[id_df.name.values[n]] = id_df.id.values[n]
    return result

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get dataId hash.')
    parser.add_argument('input', type=str, help="Input file with anaTuple")
    parser.add_argument('dataId', type=str, help="dataId")
    args = parser.parse_args()
    df = LoadDataIdFrame(args.input)
    ids = GetDataIdHash(df, args.dataId, False)
    if len(ids) == 0:
        ids = GetDataIdHash(df, args.dataId, True)
    if len(ids) != 0:
        for name in ids:
            print('{:<20}\t{}'.format(ids[name], name))
    else:
        print("ERROR: Data Id not found")
