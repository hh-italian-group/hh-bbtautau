#!/usr/bin/env python

import argparse
import re
import io
import ROOT
from uncertainties import ufloat

parser = argparse.ArgumentParser(description='Plot NLL as a function of nuissance parameter.')
parser.add_argument('--input', required=True, type=str, metavar='PATH', help="input shapes")
parser.add_argument('--output', required=True, type=str, metavar='PATH', help="output csv file with the yield table")
parser.add_argument('--threshold', required=False, default=None, metavar='PATH', help="output table")
args = parser.parse_args()

energy_scales = [ 'CMS_scale_t_13TeVUp', 'CMS_scale_t_13TeVDown', 'CMS_scale_j_13TeVUp', 'CMS_scale_j_13TeVDown',
                  'CMS_topPt_13TeVUp', 'CMS_topPt_13TeVDown' ]

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
    output.write('\ufeff')

    columns = [ 'central' ]
    columns.extend(energy_scales)
    header_str = u'Process,'
    for c in columns:
        header_str += u'{},'.format(c)
    output.write(header_str + '\n')

    for cat,dir in sorted(category_dirs.items()):
        if MatchesIgnore(cat, ignore_categories): continue
        output.write('\n' + cat + '\n')
        items = CollectItems(dir, 'TH1')
        processes = JoinByProcess(items)
        for process, p_dict in sorted(processes.items()):
            p_str = u'{},'.format(process)
            for col in columns:
                hist = p_dict[col]
                if hist == None:
                    p_str += u'-,'
                    continue
                err = ROOT.Double(0)
                integral = hist.IntegralAndError(1, hist.GetNbinsX(), err)
                if args.threshold != None and integral < float(args.threshold):
                    p_str += u'0,'
                else:
                    process_yield = ufloat(integral, err)
                    p_str += u'{:.2uP},'.format(process_yield)
            output.write(p_str + '\n')
