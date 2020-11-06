#!/usr/bin/env python

import argparse
import re
import io
import ROOT
from uncertainties import ufloat

parser = argparse.ArgumentParser(description='Plot NLL as a function of nuissance parameter.')
parser.add_argument('--input', required=True, type=str, metavar='PATH', help="input shapes")
parser.add_argument('--output', required=True, type=str, metavar='PATH', help="output csv file with the yield table")
parser.add_argument('--threshold', required=False, default=None, type=float, help="threshold")
parser.add_argument('--bin', required=False, default=None, type=int, help="consider only a single bin")
args = parser.parse_args()

years = [2016, 2017, 2018]
unc_scales = [ 'Up', 'Down' ]
unc_sources = [ 'scale_t', 'scale_j' ]

energy_scales = []
for year in years:
    for unc_source in unc_sources:
        for unc_scale in unc_scales:
            es = 'CMS_{}_13TeV{}{}'.format(unc_source, year, unc_scale)
            energy_scales.append(es)

ignore_processes = [ '^ggGraviton_hh_ttbb_M.*', '^ggRadion_hh_ttbb_M.*', '^ggh_hh_ttbb_kl([^1]|1.+)']
ignore_categories = [ '.*[OS]S_(antiiso|iso|LooseIsolated)$', '.*res2lb$' ]

def CollectItems(dir, class_name=''):
    items = {}
    dir_iter = ROOT.TIter(dir.GetListOfKeys())
    key = dir_iter.Next()
    while key != None:
        root_class = ROOT.gROOT.GetClass(key.GetClassName())
        if len(class_name) == 0 or (root_class != None and root_class.InheritsFrom(class_name)):
            items[key.GetName()] = dir.Get(key.GetName())
        key = dir_iter.Next()
    return items

def IsUncertainty(name):
    for es in energy_scales:
        if len(name) > len(es) and name[-len(es):] == es:
            return True
    return False

def MatchesIgnore(name, ignore_list):
    for pattern in ignore_list:
        if re.match(pattern, name) != None:
            return True
    return False

def JoinByProcess(items):
    processes = {}
    for process,hist in items.items():
        if IsUncertainty(process) or MatchesIgnore(process, ignore_processes): continue
        p_dict = {}
        p_dict['central'] = hist
        for es in energy_scales:
            name = '{}_{}'.format(process, es)
            p_dict[es] = items.get(name)
        processes[process] = p_dict
    return processes


shape_file = ROOT.TFile(args.input)
category_dirs = CollectItems(shape_file, 'TDirectory')

with io.open(args.output, 'w', encoding='utf-8') as output:
    output.write(u'\ufeff')

    columns = [ 'central' ]
    columns.extend(energy_scales)
    header_str = u'Process,'
    for c in columns:
        header_str += u'{},{} err,'.format(c, c)
    output.write(header_str + '\n')

    for cat,dir in sorted(category_dirs.items()):
        if MatchesIgnore(cat, ignore_categories): continue
        output.write(u'\n' + cat + u'\n')
        items = CollectItems(dir, 'TH1')
        processes = JoinByProcess(items)
        for process, p_dict in sorted(processes.items()):
            p_str = u'{},'.format(process)
            for col in columns:
                hist = p_dict[col]
                if hist == None:
                    p_str += u'-,-,'
                    continue

                if args.bin is not None:
                    if args.bin < 0:
                        max_bin = hist.GetNbinsX() + 1 + args.bin
                    else:
                        max_bin = args.bin
                    min_bin = max_bin
                else:
                    min_bin, max_bin = 1, hist.GetNbinsX()

                err = ROOT.Double(0)
                integral = hist.IntegralAndError(min_bin, max_bin, err)
                if args.threshold != None and integral < float(args.threshold):
                    p_str += u'0,-,'
                else:
                    #process_yield = ufloat(integral, err)
                    #p_str += u'{:.2uP},'.format(process_yield)
                    p_str += u'{:.4f},{:.4f},'.format(integral, err)
            output.write(p_str + u'\n')
