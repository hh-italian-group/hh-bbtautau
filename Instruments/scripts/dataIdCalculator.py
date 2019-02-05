#!/usr/bin/env python
# Calculate dataId
# This file is part of https://github.com/hh-italian-group/h-tautau.

import uproot
import pandas
import sys

file = uproot.open('tauTau_tuple.root')

aux = file['aux']

df = aux.arrays('*', outputtype=pandas.DataFrame)\

new_df = pandas.DataFrame(data = {
                          'id' : df.dataIds[0],
                          'name': df.dataId_names[0]
                          })
new_df.head()

syncDataIds = "2j/NoCuts/OS_Isolated/Central/DY_nlo"

try:
	id = new_df[new_df.name == syncDataIds].id.values[0]
	print id

except Exception as e:
	print (sys.stderr, "ERROR: Data Id not defined")
