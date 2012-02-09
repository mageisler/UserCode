import FWCore.ParameterSet.Config as cms

TrackingParticleSelectionGeneral = cms.PSet(
    lip = cms.double(30.0),
    chargedOnly = cms.bool(True),
    pdgId = cms.vint32(),
    signalOnly = cms.bool(True),
    stableOnly = cms.bool(False),
    minRapidity = cms.double(-2.4),
    minHit = cms.int32(6),
    ptMin = cms.double(1.),
    maxRapidity = cms.double(2.4),
    tip = cms.double(3.0)
)