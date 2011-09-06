import FWCore.ParameterSet.Config as cms

process = cms.Process("PFPUAssoMap")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/m/mgeisler/ValidationFiles/QCD_TuneZ2_7TeV_Spring11_0.root')
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
process.load("MGeisler.PF_PU_AssoMap.cuts_cff")
process.load("MGeisler.PF_PU_AssoMap.PF_PU_AssoMap_cff")
		
#process.Tracks2Vertex.GsfElectronCollection = cms.untracked.string('default')

process.FirstVertexTrackCollection = cms.EDProducer('FirstVertexTracks',
          TrackCollection = cms.untracked.string('cutsRecoTracks'),
          GsfElectronCollection = cms.untracked.string('gsfElectrons'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
)

### TrackAssociation-specific includes	
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')	
from MGeisler.PF_PU_AssoMap.TrackingParticleSelection_cfi import *

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

process.p = cms.Path(process.Tracks2Vertex*process.cutsRecoTracks*process.FirstVertexTrackCollection)
		
process.outpath = cms.EndPath(process.out)
