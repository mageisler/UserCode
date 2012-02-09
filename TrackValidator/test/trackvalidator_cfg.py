import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root')
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('Demo.root'),
    closeFileFast = cms.untracked.bool(True),
)

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
process.load("MGeisler.PF_PU_AssoMap.cuts_cff")
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

########### track selection configuration ########
process.cutsRecoTracks.minRapidity = cms.double(-2.5)
process.cutsRecoTracks.maxRapidity = cms.double(2.5)
process.cutsRecoTracks.quality = cms.vstring('highPurity')
process.cutsRecoTracks.tip = cms.double(3.)
process.cutsRecoTracks.lip = cms.double(30.)
process.cutsRecoTracks.ptMin = cms.double(1.)
				
process.demo = cms.EDAnalyzer('TrackValidator',
	tcLabel = cms.VInputTag(cms.InputTag("cutsRecoTracks")),
    	tcRefLabel = cms.InputTag("generalTracks"),
    	PULabel = cms.InputTag("addPileupInfo"),
    	TPLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
	ignoremissingtrackcollection=cms.bool(False),
	UseLogPt=cms.bool(True),
	generalTpSelector = TrackingParticleSelectionGeneral,
)


process.p = cms.Path(
      process.cutsRecoTracks
    * process.demo
)
