# Auto generated configuration file
# using: 
# Revision: 1.381.2.2 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: STEP2 --step RAW2DIGI,L1Reco,RECO,VALIDATION:validation_prod,DQM:DQMOfflinePOGMC --conditions START53_V7A::All --pileup 2012_Summer_50ns_PoissonOOTPU --datamix NODATAMIXER --eventcontent AODSIM,DQM --datatier AODSIM,DQM -n 10 --no_exec --mc --filein tc_RAWSIM.root --fileout tc_AODSIM.root --python_filename ggH250_RECO_cfg.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        'file:../../RAW_SIM/HHTauTauBB/HH_TauTauBB_RAWSIM.root'
        #'root://osg-se.sprace.org.br//store/user/caber//AbelianZprime_M-1500_hZ_tautauqq_LHCC-Madgraph5_22072013/AbelianZprime_M-1500_hZ_tautauqq_LHCC-Madgraph5_22072013_RAWSIM/06e6bebb6525c1c46ccfcc56d82513c0/HZ_TAUTAUJJ_RAWSIM_14_1_0mr.root'
        #'root://osg-se.sprace.org.br//store/user/caber/AbelianZprime_M-2500_hZ_tautauqq_LHCC_NEW-Madgraph5_31072013/AbelianZprime_M-2500_hZ_tautauqq_LHCC-Madgraph5_02082013_RAWSIM/06e6bebb6525c1c46ccfcc56d82513c0/HZ_TAUTAUJJ_RAWSIM_28_1_a0v.root'
    )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.2 $'),
    annotation = cms.untracked.string('STEP2 nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('HH_TauTauBB_AODSIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AODSIM')
    )
)
 
# Additional output definition

# Other statements
process.mix.playback = True
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
process.GlobalTag.globaltag = 'START53_V7A::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.dqmoffline_step = cms.Path(process.DQMOfflinePOGMC)
process.validation_step = cms.EndPath(process.validation_prod)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)
 

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.validation_step,process.dqmoffline_step,process.endjob_step,process.AODSIMoutput_step)

