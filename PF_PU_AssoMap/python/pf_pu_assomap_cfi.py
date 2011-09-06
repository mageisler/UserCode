import FWCore.ParameterSet.Config as cms

from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *
from SimTracker.TrackAssociation.CosmicParametersDefinerForTP_cfi import *
from Validation.RecoTrack.MTVHistoProducerAlgoForTrackerBlock_cfi import *

Tracks2Vertex = cms.EDProducer('PF_PU_AssoMap',
	
	 
	  #Set the Input Collections
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
          TrackCollection = cms.untracked.string('generalTracks'),
          GsfElectronCollection = cms.untracked.string('gsfElectrons'),
		  
	    
	  #Configuration for the Vertices
          VertexQuality = cms.untracked.bool(True),
          VertexMinNdof = cms.untracked.double(4.),
		  
	  #Configuration for the 2nd association
          VertexAssOneDim = cms.untracked.bool(True),
          VertexAssClosest = cms.untracked.bool(True),
          VertexAssUseAbsDistance = cms.untracked.bool(False),
		  
	  #Configuration for the GsfElectrons
          UseGsfElectronVertex = cms.untracked.bool(True),
	  UseCtfAssVertexForGsf = cms.untracked.bool(False),
	   
	  #Configuration for the reassociation of gamma conversion particles
	  ConversionsCollection = cms.InputTag('allConversions'),
	   
	  #Configuration for the reassociation of particles from V0 decays
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	   
	  #Configuration for the reassociation of particles from nuclear interactions
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
)
