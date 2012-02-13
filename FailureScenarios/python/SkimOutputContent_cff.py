import FWCore.ParameterSet.Config as cms
		
SkimOutputContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
	'keep PileupSummaryInfos_addPileupInfo_*_*',
	'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_RERECO',
	'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_RERECO',
        'keep TrackingParticles_mergedtruth_*_*',
	'keep TrackingVertexs_mergedtruth_*_*',
	'keep TrackingRecHitsOwned_ckf*_*_*',
	'keep TrackingRecHitsOwned_electronGsfTracks_*_*',
	'keep TrackingRecHitsOwned_generalTracks_*_*',
	'keep TrackingRecHitsOwned_pixelTracks_*_*',
	'keep recoBeamSpot_offlineBeamSpot_*_*',
	'keep recoConversions_allConversions_*_*',
	'keep recoElectronSeeds_electronMergedSeeds_*_*',
	'keep recoGsfElectrons_gsfElectrons_*_*',
	'keep recoGsfElectrons_pfElectronTranslator_*_*',
	'keep recoGsfElectronCores_gsfElectronCores_*_*',
	'keep recoGsfElectronCores_pfElectronTranslator_*_*',
	'keep recoGsfTracks_electronGsfTracks_*_*', 
	'keep recoGsfTrackExtras_electronGsfTracks_*_*',
	'keep recoMuons_muons_*_*',
	'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*',
	'keep recoPFMETs_pfMet_*_*',
	'keep recoTracks_ckf*_*_*',
	'keep recoTracks_generalTracks_*_*',
	'keep recoTracks_pixelTracks_*_*',
	'keep recoTracks_*Muons*_*_*',
	'keep recoTrackExtras_ckf*_*_*',
	'keep recoTrackExtras_electronGsfTracks_*_*',
	'keep recoTrackExtras_generalTracks_*_*',
	'keep recoTrackExtras_pixelTracks_*_*',
	'keep recoTrackExtras_*Muons*_*_*',
	'keep recoVertexs_offlinePrimaryVertices_*_*',
	'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*',
	'keep recoVertexs_pixelVertices_*_*',
	'keep *_ak5JetTracksAssociator*_*_*',
	'keep *_ak5PFJets_*_*',
	'keep *_ak5TrackJets_*_*',
	'keep *_dedx*_*_*',
	'keep *_csc2DRecHits_*_*',
	'keep *_cscSegments_*_*',
	'keep *_dt1D*RecHits_*_*',
	'keep *_dt4D*Segments_*_*',
	'keep *_ecalPreshowerRecHit_*_*',
	'keep *_ecalRecHit_*_*',
	'keep *_g4SimHits_*_RERECO',
	'keep *_generalV0Candidates_*_*',
	'keep *_generator_*_*',
	'keep *_hbhereco__*',
	'keep *_hfreco__*',
	'keep *_horeco__*',
        'keep *_regionalCosmicTracks_*_*',
        'keep *_rpcRecHits_*_*',
        'keep *_siPixelClusters_*_*',
        'keep *_siPixelDigis_*_*',
	'keep *_siStripClusters_*_*',
	'keep *_siStripDigis_*_*', 
        'keep *_towerMaker_*_*',
	'keep *_particleFlow_*_*'),
)
		
SkimOutputContent2 = cms.PSet(
    outputCommands = cms.untracked.vstring(
        'keep recoCaloClusters_*_*_*',
	'keep recoPhotonCores_photonCore_*_*',
	'keep recoPhotons_photons_*_*',
        'keep recoSuperClusters_*_*_*',
        'keep *_trackerDrivenElectronSeeds_*_*',
	'keep *_particleFlow_*_*',
	'keep *_pfElectronTranslator_pf_*',),
) 
 