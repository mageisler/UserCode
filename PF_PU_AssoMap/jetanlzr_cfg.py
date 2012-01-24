import FWCore.ParameterSet.Config as cms

process = cms.Process("JetAnlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root')
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)
		
OutFile = cms.string('Test.root')
		
### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_44_V12::All'

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

### AssociationMap-specific includes
process.load("MGeisler.PF_PU_AssoMap.PF_PU_AssoMap_cff")

########### track selection configuration ########
process.cutsRecoTracks.minRapidity = cms.double(-2.5)
process.cutsRecoTracks.maxRapidity = cms.double(2.5)
process.cutsRecoTracks.quality = cms.vstring('highPurity','tight','loose')
process.cutsRecoTracks.tip = cms.double(120.)
process.cutsRecoTracks.lip = cms.double(280.)
process.cutsRecoTracks.ptMin = cms.double(0.1)
				  
process.FirstVertexTrackCollection = cms.EDProducer('FirstVertexTracks',
          TrackCollection = cms.InputTag('cutsRecoTracks'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
)

########### configuration MultiTrackValidator 1 ########
process.multiTrackValidator.outputFile = OutFile
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.skipHistoFit=cms.untracked.bool(False)
process.multiTrackValidator.label = ['cutsRecoTracks','FirstVertexTrackCollection']
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.minpT = cms.double(0.1)
process.multiTrackValidator.maxpT = cms.double(3000.0)
process.multiTrackValidator.nintpT = cms.int32(40)
process.multiTrackValidator.UseAssociators = cms.bool(True)
process.multiTrackValidator.runStandalone = cms.bool(True)

process.PFCandidatesAM = cms.EDProducer('PFCand_NoPU_WithAM',
          PFCandidateCollection = cms.InputTag('particleFlow'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
          ConversionsCollection = cms.InputTag('allConversions'),
          V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
          V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
          NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex')
)
		
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
		
process.kt6PFJetsAM = kt4PFJets.clone(
      	src = cms.InputTag("PFCandidatesAM"),
    	rParam = cms.double(0.6),
    	doAreaFastjet = cms.bool(True),
    	doRhoFastjet = cms.bool(True),
    	Ghost_EtaMax = cms.double(6.5)
)

process.jetanalyzer = cms.EDAnalyzer('JetAnlzr',
	genJets = cms.string("kt6GenJets"),
	recoJets = cms.vstring("kt6PFJets"),
	PileUpInfo = cms.string("addPileupInfo")
)

# paths	
		
### produce the association map
process.T2V = cms.Sequence(
      process.Tracks2Vertex
)

### produce a collection of tracks associated to the first vertex 		
process.FVTC = cms.Sequence(
      process.FirstVertexTrackCollection
)	

### do the efficiency and fake rate analyses for the charged particles 				 
process.MTV = cms.Sequence(
      process.multiTrackValidator
)	

### produce a collection of PFCandidates associated to the first vertex		
process.PFCand = cms.Sequence(
      process.PFCandidatesAM
)	

### produce a jet collection from the PFCandidates		
process.PFJ = cms.Sequence(
      process.kt6PFJetsAM
)	

### do the jet analysis				  
process.JA = cms.Sequence(
      process.jetanalyzer
)

process.p = cms.Path(
    #* process.cutsRecoTracks
    #* process.T2V
    #* process.FVTC
    #* process.MTV
    #* process.PFCand
    #* process.PFJ
     process.JA
)					  
	       
process.schedule = cms.Schedule(
      process.p
)
