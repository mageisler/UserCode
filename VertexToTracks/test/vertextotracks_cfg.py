import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-RECO.root'),
    secondaryFileNames =  cms.untracked.vstring('file:/user/geisler/RelValZmumuJets_Pt_20_300_PU_START53_V6-v1_GEN-SIM-DIGI-RAW-HLTDEBUG.root'),
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V7F::All'
		
### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.myProducerLabel = cms.EDProducer('VertexToTracks',	
          TrackCollection = cms.InputTag('generalTracks'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  BeamSpot = cms.InputTag('offlineBeamSpot'),
          GetCleanedCollections = cms.bool(False),
          ConversionsCollection = cms.InputTag('allConversions'),
          V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
          V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
          NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),	
          ignoreMissingCollection = cms.bool(True),
          MaxNumberOfAssociations = cms.int32(3),
)

  
process.p = cms.Path(process.myProducerLabel)
