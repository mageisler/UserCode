import FWCore.ParameterSet.Config as cms

TrackingParticleSelection = cms.PSet(
    lip = cms.double(30.0),
    chargedOnly = cms.bool(True),
    stableOnly = cms.bool(False),
    pdgId = cms.vint32(),
    signalOnly = cms.bool(True),
    minRapidity = cms.double(-2.4),
    minHit = cms.int32(0),
    ptMin = cms.double(0.005),
    maxRapidity = cms.double(2.4),
    tip = cms.double(60)
)
