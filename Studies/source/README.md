# How to run
### VariableDistribution
```shell
./run.sh VariableDistribution --input_path ~/Downloads/SkimmedLight --output_file VariableDistributionSkimmed_tauTau_spin0.root --cfg_file hh-bbtautau/Studies/config/mva_config.cfg --tree_name tauTau --spin 0 --skimmed true
```
### FindOptimalBandwidth
```
./run.sh FindOptimalBandwidth ~/Downloads/SkimmedLight/ hh-bbtautau/Studies/config/mva_config.cfg 4 20000
```
### FindJSD
```
./run.sh FindJSD ~/Downloads/SkimmedLight hh-bbtautau/Studies/config/mva_config.cfg OptimalBandwidth20000 4 100
```
### FindMutual
```
./run.sh FindMutual ~/Downloads/SkimmedLight hh-bbtautau/Studies/config/mva_config.cfg OptimalBandwidth20000 4 20000
```
### MvaRangesCompatibility
```shell
 ./run.sh MvaRangesCompatibility ~/Downloads/SkimmedLight  Ranges.root  hh-bbtautau/Studies/config/mva_config.cfg OptimalBandwidthTOT JensenShannonDivergence200 MutualInformationDistance100 4 20 300000 2 1
```
### MvaVariableSelection
```shell
./run.sh MvaVariableSelection ~/Downloads/SkimmedLight  variables.root  hh-bbtautau/Studies/config/mva_config.cfg OptimalBandwidthTOT JensenShannonDivergence200 MutualInformationDistance100 4 20 300 1 0
```
### MvaTraining
```shell
./run.sh MvaTraining --input_path ~/Downloads/SkimmedLight --output_file Low21 --cfg_file hh-bbtautau/Studies/config/mva_config.cfg  --number_events 200000 --range Low --number_variables 21 --which_test 3 --seed 12345678
```
### MvaAnalyzer
```shell
./run.sh MvaAnalyzer Analyzer_Low15_muTau chi 0.05 false blind-unblind\(28Agosto\)/evaluation/muTau/EvalLow15*
```
### MvaEvaluation
```shell
./run.sh MvaEvaluation ~/Downloads/SkimmedLight EvalSM16_muTau hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactorySM16_blind_2_Grad_2.weights.xml BDT::Grad_2 200000   0 10 2 12345678 false false SM16 muTau 0
./run.sh MvaEvaluation ~/Downloads/skimmed Eval_LowFra hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_LOW.weights.xml BDT::500t_PU_mass_newvars_LOW 200000 250 350 2 12345678 1234567 true true LowF
./run.sh MvaEvaluation ~/Downloads/skimmed Eval_HighFra hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_HIGH_oldvars.weights.xml BDT::500t_PU_mass_newvars_HIGH_oldvars 20000 400 900 2 12345678 true false HighF
./run.sh MvaEvaluation ~/Downloads/skimmed Eval_SMFra hh-bbtautau/Studies/config/mva_config.cfg ~/Desktop/xml/n/TMVAClassification_500t_PU_mass_newvars_LOW.weights.xml BDT::500t_PU_mass_newvars_LOW 200000 0 10 2 12345678 true true SMF
./run.sh MvaEvaluation ~/Desktop/tuples prova_low25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryLowMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 250 320 2 12345678 1234567 false false Low25
./run.sh MvaEvaluation ~/Desktop/tuples prova_low13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryLowMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 250 320 2 12345678 1234567 false false Low13
./run.sh MvaEvaluation ~/Desktop/tuples prova_medium25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryMediumMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 340 400 2 12345678 1234567 false false Medium25
./run.sh MvaEvaluation ~/Desktop/tuples prova_medium13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryMediumMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 340 400 2 12345678 1234567 false false Medium13
./run.sh MvaEvaluation ~/Desktop/tuples prova_high25_500 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryHighMassprova25_500_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 450 900 2 12345678 1234567 false false High25
./run.sh MvaEvaluation ~/Desktop/tuples prova_high13 hh-bbtautau/Studies/config/mva_config.cfg ~/workspace/hh-analysis/mydataloader/weights/myFactoryHighMassprova13_Grad_shrinkage_0.1_12345678.weights.xml BDT::Grad_shrinkage_0.1_12345678 20000 450 900 2 12345678 1234567 false false High13

