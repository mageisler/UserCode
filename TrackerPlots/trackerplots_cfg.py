import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackerPlots")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/FS2011/RECO/QCD_Pt-15to3000_TuneZ2_Flat_PU_S6_START42_WO.root')
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.test1 = cms.EDAnalyzer("TrackerPlots",
    FileName = cms.untracked.string("Test")	
)


process.p = cms.Path(process.test1)
