# How to run
### Run merge analyzer for DYjets
./run.sh NJets_HT_BinFileMerger summary /Users/Tita/Desktop hh-bbtautau/Instruments/config/dyjets_splitting.cfg hh-bbtautau/Instruments/config/DYJets_file.cfg hh-bbtautau/Analysis/config/dyjets_weights.cfg

### Run merge analyzer for Wjets
./run.sh NJets_HT_BinFileMerger summary /Users/Tita/Desktop hh-bbtautau/Instruments/config/wjets_splitting.cfg hh-bbtautau/Instruments/config/WJets_file.cfg hh-bbtautau/Analysis/config/wjets_weights.cfg

### Run merge analyzer for TTbar
./run.sh TTFileMerger muTau /eos/user/k/kandroso/cms-it-hh-bbtautau/Tuples2016_v2 hh-bbtautau/Instruments/config/ttbar_splitting.cfg hh-bbtautau/Instruments/config/TTbar_file.cfg hh-bbtautau/Analysis/config/ttbar_weights_muTau.cfg

### Run merge analyzer for SM samples
./run.sh SMFileMerger all_events /Users/Tita/Desktop hh-bbtautau/Instruments/config/SM_node_file.cfg weight_SM.root

### Run analyzer for SM samples to check the reweighting
./run.sh SMWeight_t weight_SM.root /Users/Tita/Desktop hh-bbtautau/Instruments/config/SM_node_file.cfg eTau output_SM_reweight_eTau.root > output_SM_reweight_eTau.txt

### Run analyzer for DYjets samples to check the reweighting
./run.sh DYWeight_t summary /Users/Tita/Desktop hh-bbtautau/Analysis/config/dyjets_weights.cfg hh-bbtautau/Instruments/config/DYJets_file.cfg

### Run analyzer to test new configuration file for Run2
./run.sh AnalyzerConfig_t hh-bbtautau/Analysis/config/new_sources_Run2.cfg

### Run TupleSkimmer on gridui
./run.sh TupleSkimmer hh-bbtautau/Instruments/config/Skimmer.cfg /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2 skimmer_output #samples #setup

###Run timinig analyzer on gridui
find /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full*  -name "*.root" -not -name "*_part?.root" -exec ./run.sh TimeAnalyzer_t '{}' \; | grep -v '\[' | grep -v 'Scanning dependencies' > time_test.csv


