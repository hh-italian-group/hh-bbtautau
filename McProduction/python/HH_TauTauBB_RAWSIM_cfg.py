# Auto generated configuration file
# using: 
# Revision: 1.381.2.2 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: REDIGI --step DIGI,L1,DIGI2RAW,HLT:7E33v2 --conditions START53_V7A::All --pileup 2012_Summer_50ns_PoissonOOTPU --datamix NODATAMIXER --eventcontent RAWSIM --datatier GEN-SIM-RAW --no_exec --mc --filein tc_GENSIM.root --fileout tc_RAWSIM.root --python_filename ggH250_RAWSIM_cfg.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_7E33v2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(
        'file:/gpfs/ddn/srm/cms/store/user/mgrippo/Graviton_500_GENSIM/HH_TauTauBB_GENSIM.root' 
    )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.2 $'),
    annotation = cms.untracked.string('REDIGI nevts:3'),
    name = cms.untracked.string('PyReleaseValidation')
)

#pile-up samples
process.mix.input.fileNames = cms.untracked.vstring(
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/2C39E694-ED5D-E111-88F8-003048F0E826.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/F63E6236-505E-E111-B9F1-00266CF330B8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/BC3A84DE-2D5E-E111-86C0-003048C68A9A.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/F6D158F1-CF5D-E111-A9CA-0030487D5D8D.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/98131E75-2F5E-E111-A90C-003048C693FE.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/E09D4118-0A5E-E111-BDDB-0030487D5DB7.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/22C99365-525E-E111-A636-00266CF330B8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/D856B590-EA5D-E111-95B9-002481E94C56.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/D6D0A79F-705E-E111-8FED-002481E101DA.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/30FFE650-205E-E111-9FAF-003048C6617E.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/98DCBF67-7C5E-E111-B94E-003048F02CB8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/707CA27F-DD5D-E111-92A6-00266CF33318.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/DE55E439-1C5E-E111-B2B0-00266CF32CD0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/54B42ECA-355E-E111-98BD-0030487F1BD1.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/C8E28663-1A5E-E111-86C8-003048C68A92.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/8046216E-595E-E111-8137-0030487D8581.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/02446A08-515E-E111-82C7-00266CF330B8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/569AD9E5-5E5E-E111-B738-0030487D5D7B.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/C06CE153-6B5E-E111-957D-0030487F929D.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/1C8292E8-0E5E-E111-9CA2-0030487D811F.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/024BBC60-CD5D-E111-A41A-0025901D4C32.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/20DB039F-E45D-E111-A8BE-0030487FA629.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0000/BEAC9726-3B5E-E111-A197-0030487F1A57.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/EE4904D1-4D62-E111-8B1A-0030487D5DBF.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/2E39C88E-0964-E111-8E30-002481E0DBE0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/9A0E6D37-5962-E111-8E64-002481E0D6EE.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/0A2E8C12-8863-E111-9DB5-003048C68FEC.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/C61DB9DA-8763-E111-8BA0-003048D4610C.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/32962EBD-8763-E111-BEBC-003048D439A0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/C4BCB4D4-CC63-E111-BFA0-0030487D5EAF.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/78F58D12-8863-E111-9BAE-003048C68FEC.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/848EF462-0864-E111-9883-0025901D4AF0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/B2E79055-0A65-E111-907D-00266CF32684.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/F88F3BAF-8763-E111-8438-003048D437C4.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/943BAE89-0A65-E111-8E13-003048F0E184.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/1A7AAFEB-BA64-E111-9BFF-003048C69288.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/AAB1F78B-AB64-E111-88C7-003048C69412.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/BC62210F-4664-E111-A1DB-003048C6931E.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/6E7D7C62-F863-E111-A278-00266CFFA750.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/44213513-E561-E111-9578-0030487F1665.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/D2379210-7C62-E111-B501-003048D4DCD8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/786F9F94-0A65-E111-BF59-00266CF32E78.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/1205F429-6262-E111-9F30-0030487D83B9.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/94FEFBB6-3B64-E111-A964-003048D436C6.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/CA23348B-1064-E111-A2F3-0030487D5DC7.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/188B3433-7D62-E111-A4C6-0030487F938F.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0001/4647374F-9F64-E111-845C-003048D4364C.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/50581B54-9867-E111-80EE-00266CFFA658.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/20174A7C-9967-E111-9EF5-00266CFFA6DC.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/943FE2CF-D066-E111-89DB-003048D462FE.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/5A55BC92-9867-E111-91C1-0025901D4AF0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/00024E74-D96A-E111-9CE8-0030487D5DC3.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/32A86CBF-5C69-E111-B259-0025901D4936.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/5A5ECDAB-E267-E111-BDEB-003048C69406.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/D64D3E92-9867-E111-BEDB-00266CFFA418.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/98E497AC-E967-E111-93FB-003048CF6338.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/B0FC282C-E26A-E111-8240-0030487F4B8B.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/5267F11C-9867-E111-939E-00266CF32E70.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/4EFBF2B4-4D6A-E111-9EF3-00266CFFA344.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/6230C266-196B-E111-B29A-00266CF2ABA8.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/7654E151-3968-E111-A2D3-003048C69046.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/3A62F816-0C6A-E111-81C4-0030487F16FB.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/A4CCF79A-3368-E111-9772-002481E0E440.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/044C8BDF-516A-E111-AFFF-003048D436D2.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/7ECE561F-9867-E111-AA31-00266CF32E70.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/6681B621-2768-E111-AE10-003048C662D4.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/AA8EF528-8E69-E111-8941-003048C692BA.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/70B77665-1F68-E111-B9C3-003048D43958.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/FA16C8E7-9967-E111-9098-00266CF327C0.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/10C4B861-9967-E111-A6B0-0025901D493A.root",
"/store/mc/Summer12/MinBias_TuneZ2star_8TeV-pythia6/GEN-SIM/START50_V13-v3/0002/CAFF85A4-8C69-E111-A71F-002481E14D72.root",
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('HH_TauTauBB_RAWSIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)

# Additional output definition

# Other statements
# customise the HLT menu for running on MC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC
process = customizeHLTforMC(process)

process.GlobalTag.globaltag = 'START53_V7A::All'

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])

