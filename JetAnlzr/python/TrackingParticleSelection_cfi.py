import FWCore.ParameterSet.Config as cms


TrackingParticleSelection = cms.PSet(
    lip = cms.double(30.),
    chargedOnly = cms.bool(True),
    stableOnly = cms.bool(False),
    pdgId = cms.vint32(),
    signalOnly = cms.bool(True),
    minRapidity = cms.double(-2.4),
    minHit = cms.int32(6),
    ptMin = cms.double(1.0),
    maxRapidity = cms.double(2.4),
    tip = cms.double(3.)
)
