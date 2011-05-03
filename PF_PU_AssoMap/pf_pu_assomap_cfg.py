import FWCore.ParameterSet.Config as cms

process = cms.Process("PFPUAssoMap")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/m/mgeisler/ValidationFiles/QCD_TuneZ2_7TeV_Spring11_0.root')
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)

process.Tracks2Vertex = cms.EDProducer('PF_PU_AssoMap',
          TrackCollection = cms.InputTag('generalTracks'),
          GsfElectronCollection = cms.InputTag('gsfElectrons'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
          VertexQuality = cms.untracked.bool(True),
          VertexMinNdof = cms.untracked.double(4.),
          ClosestVertex = cms.untracked.bool(True),
)

process.FirstVertexTrackCollection = cms.EDProducer('FirstVertexTracks',
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.Tracks2Vertex*process.FirstVertexTrackCollection)
		
process.outpath = cms.EndPath(process.out)
