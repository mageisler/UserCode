import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root'
"file:/user/kuessel/CMSSW/CMSSW_4_2_5/src/UserCode/RWTH3b/GridTools/mergedTTJetsSummer11_1.root"
    )
)

process.myProducerLabel = cms.EDProducer('PU_counter',
                                         OutName = cms.untracked.string('blablabla.root'),
                                         NumberOfBins=cms.untracked.int32(36)
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('Test.root')
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath()
