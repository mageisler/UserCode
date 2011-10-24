import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root'
    )
)

process.myProducerLabel = cms.EDProducer('PU_counter',
          OutName = cms.untracked.string('Lumi_GluGluToHiggsToGG.root')
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('Test.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath()
