# How to run
### Run merge analyzer for DYjets
./run.sh NJets_HT_BinFileMerger --tree_name summary --input_path /Users/Tita/Desktop --cfg_name hh-bbtautau/Instruments/config/dyjets_splitting.cfg --file_cfg_name hh-bbtautau/Instruments/config/DYJets_file.cfg --output_file hh-bbtautau/Analysis/config/dyjets_weights.cfg

### Run merge analyzer for Wjets
./run.sh NJets_HT_BinFileMerger --tree_name summary --input_path /Users/Tita/Desktop --cfg_name hh-bbtautau/Instruments/config/wjets_splitting.cfg --file_cfg_name hh-bbtautau/Instruments/config/WJets_file.cfg --output_file hh-bbtautau/Analysis/config/wjets_weights.cfg

### Run merge analyzer for TTbar
./run.sh TTFileMerger --tree_name muTau --input_path /eos/user/k/kandroso/cms-it-hh-bbtautau/Tuples2016_v2 --cfg_name hh-bbtautau/Instruments/config/ttbar_splitting.cfg --file_cfg_name hh-bbtautau/Instruments/config/TTbar_file.cfg --output_file hh-bbtautau/Analysis/config/ttbar_weights_muTau.cfg

### Run merge analyzer for SM samples
./run.sh SMFileMerger --tree_name all_events --input_path /Users/Tita/Desktop --file_cfg_name hh-bbtautau/Instruments/config/SM_node_file.cfg --output_file weight_SM.root

### Run analyzer for SM samples to check the reweighting
./run.sh SMWeight_t --weight_file weight_SM.root --input_path /Users/Tita/Desktop --file_cfg_name hh-bbtautau/Instruments/config/SM_node_file.cfg --tree_name eTau --output_file output_SM_reweight_eTau.root > output_SM_reweight_eTau.txt

### Run analyzer for DYjets samples to check the reweighting
./run.sh DYWeight_t --tree_name summary --input_path /Users/Tita/Desktop --cfg_name hh-bbtautau/Analysis/config/dyjets_weights.cfg --file_cfg_name hh-bbtautau/Instruments/config/DYJets_file.cfg

### Run analyzer to test new configuration file for Run2
./run.sh AnalyzerConfig_t --file_cfg_name hh-bbtautau/Analysis/config/new_sources_Run2.cfg

### Run TupleSkimmer on gridui
./run.sh TupleSkimmer --cfg hh-bbtautau/Instruments/config/Skimmer.cfg --inputPath /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2 --outputPath skimmer_output --jobs DYJets,EWK (or all) --setup_name mh_setup

###Run timinig analyzer on gridui
find /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full*  -name "*.root" -not -name "*_part?.root" -exec ./run.sh ExeTimeAnalyzer '{}' \; | grep -v '\[' | grep -v 'Scanning dependencies' > time_test.csv


