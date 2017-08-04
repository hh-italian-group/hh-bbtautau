# How to run
### VariableDistribution
```shell
./run.sh VariableDistribution --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full --output_file VariableDistribution_muTau.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau
```
### MvaRangesCompatibility
```shell
./run.sh MvaRangesCompatibility --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full --output_file MassVariables__muTau_2.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau  --number_threads 4 --number_variables 20 --number_events 10000 --number_sets 3 --set 2
```
### MvaVariableSelection
```shell
./run.sh MvaVariableSelection --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full --output_file MassVariables__muTau_2.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name muTau  --number_threads 4 --number_variables 20 --number_events 10000 --number_sets 3 --set 2
```
### MvaTraining
```shell
./run.sh MvaTraining --input_path /gpfs/ddn/cms/user/androsov/store/cms-it-hh-bbtautau/Tuples2016_v2/Full --output_file LowMass --cfg_file hh-bbtautau/Studies/config/mva_config.cfg  --number_events 20000 --range Low --seed 12345678
```
### MvaAnalyzer
```shell
./run.sh MvaAnalyzer  --output_file analyzer_2345678.root --input_file BDT_lm_2345678.root
```
### MvaEvaluation
```shell
./run.sh MvaEvaluation ~/Desktop/tuples Eval_LowFra hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_LOW.weights.xml BDT::500t_PU_mass_newvars_LOW 200000 250 350 2 12345678 1234567 true true LowF
./run.sh MvaEvaluation ~/Desktop/tuples Eval_HighFra hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_HIGH_oldvars.weights.xml BDT::500t_PU_mass_newvars_HIGH_oldvars 20000 400 900 2 12345678 1234567 true false HighF
./run.sh MvaEvaluation ~/Desktop/tuples prova_low25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryLowMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 250 320 2 12345678 1234567 false false Low25
./run.sh MvaEvaluation ~/Desktop/tuples prova_low13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryLowMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 250 320 2 12345678 1234567 false false Low13
./run.sh MvaEvaluation ~/Desktop/tuples prova_medium25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryMediumMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 340 400 2 12345678 1234567 false false Medium25
./run.sh MvaEvaluation ~/Desktop/tuples prova_medium13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryMediumMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 340 400 2 12345678 1234567 false false Medium13
./run.sh MvaEvaluation ~/Desktop/tuples prova_high25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryHighMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 450 900 2 12345678 1234567 false false High25
./run.sh MvaEvaluation ~/Desktop/tuples prova_high13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryHighMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 450 900 2 12345678 1234567 false false High13
```
