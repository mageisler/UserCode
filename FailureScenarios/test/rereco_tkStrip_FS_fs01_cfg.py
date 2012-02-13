# Auto generated configuration file
# using: 
# Revision: 1.342 
# Source: /cvs/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: tkStrip_fs08 -s SIM,DIGI,L1,DIGI2RAW,HLT,RAW2DIGI,RECO,ENDJOB --no_exec -n -1 --filein file:/tmp/mgeisler/QCD_Pt_15to3000_TuneZ2_Flat_8TeV_pythia6_GEN-SIM_FS08.root --fileout file:/tmp/mgeisler/fs_08_test.root --process RERECO --conditions auto:startup
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_50ns_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
        'file:/user/geisler/FS2011/RAW/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_GEN-SIM-RECO.root'
] ); 

secFiles.extend( [
        'file:/user/geisler/FS2011/RAW/RelValMinBias_CMSSW_4_4_2-START44_V7-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'
] ); 

# Input source
process.source = cms.Source("PoolSource",
    fileNames = readFiles,
    secondaryFileNames = cms.untracked.vstring(),   
    #eventsToProcess = cms.untracked.VEventRange('1:3810-1:3813'),
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.2 $'),
    annotation = cms.untracked.string('tkStrip_fs08 nevts:-1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,   
    fileName = cms.untracked.string('QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_RERECO_FS01.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)
		
process.AddOutput = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep TrackingParticles_mergedtruth_*_*',
	'keep TrackingVertexs_mergedtruth_*_*',
	'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*',
	'keep StripDigiSimLinkedmDetSetVector_*_*_*',
	'drop *_*_*_RECO',
	'drop PileupSummaryInfos_addPileupInfo__HLT'),
)
		
process.RECOSIMoutput.outputCommands.extend(process.AddOutput.outputCommands)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START44_V10::All'
		
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("SiStripBadModuleRcd"),
           tag = cms.string("SiStripBadComponents_FailureScenario_PS1"),
           connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_STRIP")
          )
) 
	
#process.mix.input.fileNames = secFiles
		
from SimGeneral.TrackingAnalysis.trackingParticles_cfi import *

# Path and EndPath definitions
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.validation_step = cms.Path(process.trackingParticles)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.RECOSIMoutput_step])

