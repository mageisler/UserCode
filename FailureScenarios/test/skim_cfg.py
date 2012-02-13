# configuration file to skim the datafiles to a usable size  

import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process('SKIM')

# import of standard configurations
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('MGeisler.FailureScenarios.SkimOutputContent_cff')
process.load('MGeisler.FailureScenarios.Draft_FileList_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring(),   
)


# Additional input definition	

# Output definition
process.SKIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(10000),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring('drop *_*_*_*'),   
    fileName = cms.untracked.string('test.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)
		
#Additional settings
input='/user/geisler/FS2011/FileLists/InputFileList.txt'
abbr="Draft"
		
if len(sys.argv) > 2:
    OutDir = 'file:/user/geisler/FS2011/RECO/Skim/'
    abbr = sys.argv[2]
    if len(sys.argv) > 3:
        process.maxEvents.input= int(sys.argv[3])
        if len(sys.argv) > 4:
            if sys.argv[4] == "True":
	        OutDir = 'file:/user/geisler/FS2011/RECO/Tmp/'
    print "input set to "+abbr
    if abbr == "Draft":
        OutName= OutDir +'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_Draft.root'
    if abbr == "FS01":
        OutName= OutDir +'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_FS01.root'
    if abbr == "FS08":
        OutName= OutDir +'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_FS08.root'
    if abbr == "FS09":
        OutName= OutDir +'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_FS09.root'
    if abbr == "FS10":
        OutName= OutDir +'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_FS10.root'
    process.SKIMoutput.fileName = OutName
    print "outputfile set to "+OutName+" \n"

from MGeisler.FailureScenarios.sourceFiles_cfi import *
process.source.fileNames.extend(GetFileNames(abbr, input))
	
process.SKIMoutput.outputCommands.extend(process.SkimOutputContent.outputCommands)

# Other statements
process.GlobalTag.globaltag = 'START44_V10::All'

# Path and EndPath definitions
process.output_step = cms.EndPath(process.SKIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.output_step)