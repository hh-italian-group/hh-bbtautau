import argparse
import sys, os
import subprocess
import pandas as pd

pathCsv = 'sample_list.csv'

# convert csv to dataframe
dfIn = pd.read_csv(pathCsv, names=['a'], skiprows=0)

# take column that contains sample names
samples = dfIn['a']

defList = ["e","mu","tau"]
# defYear = ["2017","2018"]
defYear = ["2016","2017","2018"]

parser = argparse.ArgumentParser()
# parser.add_argument("-s", "--sample", nargs='+')
parser.add_argument("-l", "--lep", nargs='+', default = defList)
parser.add_argument("-y", "--year", nargs='+', default = defYear)

args = parser.parse_args()

# pwd = os.path.dirname(os.path.realpath(__file__))
pwd = "/gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau"
spin = -1
mass = -1
node = -1
# for sample in args.sample:
for i, sample in enumerate(samples):
    print(sample)
    for lepton in args.lep:
        if lepton not in defList:
            print(" \"%s\" is not a valid lepton" % lepton)
            sys.exit()
        # print(lepton)
        for year in args.year:
            if year not in defYear:
                print(" \"%s\" is not a valid year" % year)
                sys.exit()
            path = "%s/Tuples%s_v5/Full/" %(pwd, year)
            sample_name_and_path = "%s%s" %(path, sample)

            #ggF part
            if sample[:6] == "GluGlu" :
                sample_v2 = sample[6:]
                process = sample_v2[:4]
                #Non resonant
                if process == 'ToHH' :
                    process_name  = 'ggHH_NonRes'
                    node_pos = sample_v2.find('node_')
                    node = sample_v2[node_pos+5:]
                    if node == 'SM' :
                        node = 0
                    elif node == 'box' :
                        node = 1
                #Resonant
                elif process == 'ToBu' or process == 'ToRa' :
                    process_name  = 'ggHH_Res'
                    #Graviton
                    if process == 'ToBu':
                        spin = 2
                        m_pos = sample_v2.find('-')
                        mass = sample_v2[m_pos+1:-7]
                    #Radion
                    elif process == 'ToRa' :
                        spin = 0
                        m_pos = sample_v2.find('-')
                        mass = sample_v2[m_pos+1:-7]
            #VBF part
            elif sample[:3] == "VBF" :
                sample_v2 = sample[3:]
                process = sample_v2[:4]
                #Non resonant
                if process == 'HHTo' :
                    process_name  = 'VBFHH_NonRes'
                    CV_pos = sample_v2.find('CV_')
                    C2V_pos = sample_v2.find('C2V_')
                    C3_pos = sample_v2.find('C3_')

                    CV_node = sample_v2[CV_pos+3:C2V_pos-1]
                    if len(CV_node) == 3 :
                        CV_node_1 = sample_v2[CV_pos+3:C2V_pos-3]
                        CV_node_2 = sample_v2[CV_pos+5:C2V_pos-1]
                    elif len(CV_node) == 1 :
                        CV_node_1 = sample_v2[CV_pos+3:C2V_pos-1]
                        CV_node_2 = 0

                    C2V_node = sample_v2[C2V_pos+4:C3_pos-1]
                    C3_node = sample_v2[C3_pos+3:]

                    if CV_node_1 == '0' and CV_node_2 == '5' and C2V_node == '1' and C3_node == '1' :
                        node = 1
                    elif CV_node_1 == '1' and CV_node_2 == '5' and C2V_node == '1' and C3_node == '1' :
                        node = 2
                    elif CV_node_1 == '1' and C2V_node == '1' and C3_node == '0' :
                        node = 3
                    elif CV_node_1 == '1' and C2V_node == '1' and C3_node == '1' :
                        node = 4
                    elif CV_node_1 == '1' and C2V_node == '1' and C3_node == '2' :
                        node = 5
                    elif CV_node_1 == '1' and C2V_node == '2' and C3_node == '1' :
                        node = 6
                    else :
                        print(" \"%s\" is not a valid sample" )
                        sys.exit()
                #Resonant
                elif process == 'ToBu' or process == 'ToRa' :
                    process_name  = 'VBFHH_Res'
                    #Graviton
                    if process == 'ToBu':
                        spin = 2
                        m_pos = sample_v2.find('-')
                        mass = sample_v2[m_pos+1:-7]
                    #Radion
                    elif process == 'ToRa' :
                        spin = 0
                        m_pos = sample_v2.find('-')
                        mass = sample_v2[m_pos+1:-7]

            #  args =  sample_path, channel, output_file, pdg_info, process_name, spin, mass, node, year,
            cmd = "./run.sh GenStudy --inputPath=%s.root --channel=%sTau --outputFile=test_%sTau.root --particleNameTypeFile=pdg_name_type_charge.txt --sample_type=%s --spin=%s --mass_point=%s --node=%s --year=%s --new_output_file=NN_samples/"  \
                  % (sample_name_and_path, lepton, lepton, process_name, spin, mass, node, year)
            print(cmd)
            os.system(cmd)
