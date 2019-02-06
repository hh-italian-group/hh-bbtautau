#!/usr/bin/env python
# Calculate dataId
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

import uproot
import pandas
import sys

def LoadDataIdFrame(file):
    """Load pandas DataFrame from root file."""
    file = uproot.open(file)

    aux = file['aux']

    df = aux.arrays('*', outputtype=pandas.DataFrame)\

    names = [ name.decode("utf-8") for name in df.dataId_names[0] ]
    new_df = pandas.DataFrame(data = {
                              'id' : df.dataIds[0],
                              'name': names
                              })
    return new_df

def GetDataIdHash(df, dataId):
    """Retrieve hash for the given dataId."""
    id_df = df[df.name == dataId]
    if id_df.shape[0] != 1:
        raise RuntimeError("ERROR: Data Id not defined")
    return id_df.id.values[0]

if __name__ == '__main__':
    if len(sys.argv) == 3 :
        file = sys.argv[1]
        dataId = sys.argv[2]
        df = LoadDataIdFrame(file)
        id = GetDataIdHash(df, dataId)
        print(id)
    else:
        print("Two arguments are needed")
