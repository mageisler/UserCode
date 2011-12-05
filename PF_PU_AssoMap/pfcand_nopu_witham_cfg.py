import FWCore.ParameterSet.Config as cms

process = cms.Process("PFCAND")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root')
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All'

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### AssociationMap-specific includes
process.load("MGeisler.PF_PU_AssoMap.PF_PU_AssoMap_cff")

process.PFCand = cms.EDProducer('PFCand_NoPU_WithAM',
          PFCandidateCollection = cms.InputTag('particleFlow'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
          ConversionsCollection = cms.InputTag('allConversions'),
          V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
          V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
          NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex')
)

  
process.p = cms.Path(  process.Tracks2Vertex
		     * process.PFCand
)
		
#process.myOutput = cms.OutputModule("PoolOutputModule",
     	#fileName = cms.untracked.string('myOutput.root')
#)
  
#process.e = cms.EndPath( process.myOutput )
