import FWCore.ParameterSet.Config as cms

process = cms.Process("REASSOCTEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/user/geisler/QCD_Pt-600to800_TuneZ2star_8TeV_PU_S7_START52_V9-v1_GEN-SIM-RECODEBUG.root'
    )
)
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

process.TFileService = cms.Service('TFileService',
    fileName = cms.string("ReassociationTest.root"),
    closeFileFast = cms.untracked.bool(True),
)
				
process.Uncleaned = cms.EDAnalyzer('ReassocTest',
	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
	TpSelector = TrackingParticleSelectionGeneral,
        TrackCollection = cms.InputTag('generalTracks'),
        VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	ConversionsCollection = cms.InputTag('allConversions'),
	V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	UseBeamSpotCompatibility = cms.untracked.bool(False),
	BeamSpot = cms.InputTag('offlineBeamSpot'),
        VertexAssOneDim = cms.untracked.bool(True),
        VertexAssClosest = cms.untracked.bool(True),
        VertexAssUseAbsDistance = cms.untracked.bool(False),
        ignoreMissingCollection = cms.bool(True),
	#Set the pt cut for the tracks
        TrackPtCut = cms.double(1000.),
        BeamSpotCompatibilityCut = cms.double(5.),
        nTrackWeight = cms.double(0.01),
    	PULabel = cms.InputTag("addPileupInfo"),
        GetCleanedCollections = cms.untracked.bool(False),
)
				
process.Cleaned = cms.EDAnalyzer('ReassocTest',
	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
	TpSelector = TrackingParticleSelectionGeneral,
        TrackCollection = cms.InputTag('generalTracks'),
        VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	ConversionsCollection = cms.InputTag('allConversions'),
	V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	UseBeamSpotCompatibility = cms.untracked.bool(False),
	BeamSpot = cms.InputTag('offlineBeamSpot'),
        VertexAssOneDim = cms.untracked.bool(True),
        VertexAssClosest = cms.untracked.bool(True),
        VertexAssUseAbsDistance = cms.untracked.bool(False),
        ignoreMissingCollection = cms.bool(True),
	#Set the pt cut for the tracks
        TrackPtCut = cms.double(1000.),
        BeamSpotCompatibilityCut = cms.double(5.),
        nTrackWeight = cms.double(0.01),
    	PULabel = cms.InputTag("addPileupInfo"),
        GetCleanedCollections = cms.untracked.bool(True),
)


process.p = cms.Path(process.Uncleaned + process.Cleaned)
