import FWCore.ParameterSet.Config as cms

from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

trackvalidatorDemo = cms.EDAnalyzer('TrackValidator',
	tcLabel = cms.VInputTag(cms.InputTag("generalTracks")),
    	tcRefLabel = cms.InputTag("generalTracks"),
	pfLabel = cms.VInputTag(),
    	pfRefLabel = cms.InputTag("particleFlow"),
    	PULabel = cms.InputTag("addPileupInfo"),
    	TPLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
	ignoremissingtrackcollection=cms.bool(False),
	photonPtMin=cms.double(1.0),
	photonEtaMin=cms.double(-2.4),
	photonEtaMax=cms.double(2.4),
	photonLip=cms.double(30.),
	photonTip=cms.double(3.),
	UseLogPt=cms.bool(False),
	generalTpSelector = TrackingParticleSelectionGeneral,
	useJetWeighting = cms.untracked.bool(True),
	genJetCollLabel = cms.untracked.string("ak5GenJets"),
	jetCollLabel = cms.untracked.string("ak5PFJets"),
)
