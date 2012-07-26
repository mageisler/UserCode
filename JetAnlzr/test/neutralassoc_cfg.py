import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/user/geisler/QCD_Pt-600to800_TuneZ2star_8TeV_PU_S7_START52_V9-v1_GEN-SIM-RECODEBUG.root'
    )
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string("NeutralAssoc.root"),
    closeFileFast = cms.untracked.bool(True),
)	
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(False),
)
		
process.demo = cms.EDAnalyzer('NeutralAssoc',
    	PileUp = cms.InputTag("addPileupInfo"),
    	PFCandidates = cms.InputTag("particleFlow"),
    	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
    	VertexCollection = cms.InputTag("selectedPrimaryVertexQuality"),
    	BeamSpot = cms.InputTag("offlineBeamSpot"),
)


process.p = cms.Path(process.selectedPrimaryVertexQuality + process.demo)
