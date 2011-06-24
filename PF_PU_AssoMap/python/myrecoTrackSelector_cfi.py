import FWCore.ParameterSet.Config as cms
 
myrecoTrackSelector = cms.EDFilter("RecoTrackSelector",
    src = cms.InputTag("generalTracks"),
    maxChi2 = cms.double(10000.0),
    tip = cms.double(3.5),
    minRapidity = cms.double(-2.5),
    lip = cms.double(30.0),
    ptMin = cms.double(0.9),
    maxRapidity = cms.double(2.5),
    quality = cms.vstring('highPurity'),
    algorithm = cms.vstring(),
    minHit = cms.int32(3),
    min3DHit = cms.int32(0),
    beamSpot = cms.InputTag("offlineBeamSpot")
)