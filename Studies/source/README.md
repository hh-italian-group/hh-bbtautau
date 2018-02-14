# How to run

### VariableDistribution: 1D and 2D distributions for a given channel(--tree_name) and a given value of spin(--spin)
```shell
./run.sh VariableDistribution --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file VariableDistributionSkimmed_tauTau_spin0.root --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --tree_name tauTau --spin 0 --skimmed true
```

### FindOptimalBandwidth: for single resonant samples
```
./run.sh FindOptimalBandwidth --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM false
```
### FindOptimalBandwidth: for resonant samples joined together in a given range(--which_range)
```
./run.sh FindOptimalBandwidth --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg  --number_threads 4 --range true --which_range 250 --number_events 15000 --is_SM false
```
### FindOptimalBandwidth: for non resonant sample(--is_SM)
```
./run.sh FindOptimalBandwidth --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM true
```

### FindMutual: for single resonant samples
```
./run.sh FindMutual --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM false
```
### FindMutual: for resonant samples joined together in a given range(--which_range)
```
./run.sh FindMutual --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range true --which_range 250 --number_events 15000 --is_SM false
```
### FindMutual: for non resonant sample(--is_SM)
```
./run.sh FindMutual --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM true
```

### FindJSD: for single resonant samples for sgn vs bkg(--bkg_vs_sgn)
```
./run.sh FindJSD --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM false --bkg_vs_sgn true
```
### FindJSD: for resonant samples joined together in a given range(--which_range) for sgn vs bkg(--bkg_vs_sgn)
```
./run.sh FindJSD --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range true --which_range 250 --number_events 15000 --is_SM false --bkg_vs_sgn true
```
### FindJSD: for non resonant sample(--is_SM) for sgn vs bkg(--bkg_vs_sgn)
```
./run.sh FindJSD --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range false --which_range 0 --number_events 15000 --is_SM true --bkg_vs_sgn true
```
### FindJSD: for resonant samples joined together in a given range(--which_range) for sgn vs sgn(--bkg_vs_sgn)
```
./run.sh FindJSD --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --number_threads 4 --range true --which_range 250 --number_events 15000 --is_SM false --bkg_vs_sgn false
```

### MvaRangesCompatibility
```shell
 ./run.sh MvaRangesCompatibility --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file RangesCompatibility.root  --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand  --jsd_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/JensShan --mutual_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/Mutual --number_threads 4 --number_variables 20 --number_events 300000 --number_sets 1 --set 0
```

### MvaVariableSelection: for resonant samples joined together in ranges
```shell
./run.sh MvaVariableSelection --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light  --output_file VariableSelection.root  --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand --jsd_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/JensShan --mutual_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/Mutual --number_threads 4 --number_variables 20 --number_events 300000 --number_sets 1 --set 0 --is_SM false
```
### MvaVariableSelection: for non resonant sample(--is_SM)
```shell
./run.sh MvaVariableSelection --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light  --output_file VariableSelection.root  --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --optband_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/OptBand --jsd_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/JensShan --mutual_folder /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/Mutual --number_threads 4 --number_variables 20 --number_events 300000 --number_sets 1 --set 0 --is_SM true
```

### MvaTraining: for resonant samples joined together in a given range (--range), with half part of the events blind(--blind) and the other one divided in a given number of parts(--subdivisions) and one of these is for testing(--which_test)
```shell
./run.sh MvaTraining --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Training --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg  --number_events 300000 --range Low --number_variables 20 --which_test 0 --seed 12345678 --all_data false --blind 1 --subdivisions 2 --is_SM false
```
### MvaTraining: for non resonant sample(--is_SM), with half part of the events blind(--blind) and the other one divided in a given number of parts(--subdivisions) and one of these is for testing(--which_test)
```shell
./run.sh MvaTraining --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Training --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg  --number_events 300000 --range SM --number_variables 20 --which_test 0 --seed 12345678 --all_data false --blind 1 --subdivisions 2 --is_SM true
```

### MvaAnalyzer
```shell
./run.sh MvaAnalyzer --output_file Analyzer_Low15_muTau --which_test chi --cut 0.05 --cross_validation true --input_file Training_Low_*
```

### MvaEvaluation: for resonant samples joined together in a given range (--min, --max) with a defined set of variables(--range) written in the config file (--cfg_file), with half part of the events blind(--blind). The evaluation is done for a given channel (--channel) and a given value of spin (--spin). .xml file(--file_xml) and method name (--method_name) are needed.
```shell
./run.sh MvaEvaluation --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Eval --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --file_xml  /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Analysis/config/mva/BDT-LowMass_train0-29_Grad_1.weights.xml --method_name BDT::Grad_1 --number_events 200000 --number_sets 2 --seed 12345678 --min 250 --max 320 --channel muTau --spin 0 --isLegacy false --isLow false --range Low20 --is_SM false --all_data 0 --blind 1 --subdivision 2 --which_test 1
```
### MvaEvaluation: for non resonant sample with a defined set of variables(--range) written in the config file (--cfg_file), with half part of the events blind(--blind). The evaluation is done for a given channel (--channel) and a given value of spin (--spin). .xml file(--file_xml) and method name (--method_name) are needed.
```shell
./run.sh MvaEvaluation --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Eval --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --file_xml  /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Analysis/config/mva/BDT-SM_train0-44_Grad_3.weights.xml --method_name BDT::Grad_3 --number_events 200000 --number_sets 2 --seed 12345678 --min 0 --max 10 --channel muTau --spin 0 --isLegacy false --isLow false --range SM20 --is_SM true --all_data 0 --blind 1 --subdivision 4 --which_test 3
```
### MvaEvaluation: for legacy resonant samples joined together in a given range (--min, --max) with a defined set of variables(--range) written in the config file (--cfg_file), with half part of the events blind(--blind). Legacy has to be specified(--isLegacy) as well as if it is low range(--isLow == true). The evaluation is done for a given channel (--channel) and a given value of spin (--spin). .xml file(--file_xml) and method name (--method_name) are needed.
```shell
./run.sh MvaEvaluation --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Eval --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --file_xml  /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Analysis/config/mva/HIG-17-002-BDT-LowMass.xml --method_name BDT::Grad_1 --number_events 200000 --number_sets 2 --seed 12345678 --min 250 --max 350 --channel muTau --spin 0 --isLegacy true --isLow true --range LowAN --is_SM false --all_data 0 --blind 1
```
### MvaEvaluation: for legacy non resonant sample  with a defined set of variables(--range) written in the config file (--cfg_file), with half part of the events blind(--blind). Legacy has to be specified(--isLegacy) as well as if it is low range or SM(--isLow == true). The evaluation is done for a given channel (--channel) and a given value of spin (--spin). .xml file(--file_xml) and method name (--method_name) are needed.
```shell
./run.sh MvaEvaluation --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v3/Skimmed_light --output_file Eval --cfg_file /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Studies/config/mva_config.cfg --file_xml  /gpfs/ddn/cms/user/giraldi/workspace/CMSSW_9_0_0/src/hh-bbtautau/Analysis/config/mva/HIG-17-002-BDT-LowMass.xml --method_name BDT::Grad_3 --number_events 200000 --number_sets 2 --seed 12345678 --min 0 --max 10 --channel muTau --spin 0 --isLegacy true --isLow true --range SMAN --is_SM true --all_data 0 --blind 1
```
