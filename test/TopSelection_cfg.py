
############################
#
#   ADAPT TEMPLATE
#
############################

# import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *

## remove MC matching from the default sequence to make it run on real data
removeMCMatching(process, ['All'])

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = cms.ESSource("PoolDBESSource",
                      DBParameters = cms.PSet(
                        authenticationPath = cms.untracked.string('.'),
                        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
                        idleConnectionCleanupPeriod = cms.untracked.int32(10),
                        messageLevel = cms.untracked.int32(0),
                        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
                        enableConnectionSharing = cms.untracked.bool(True),
                        connectionRetrialTimeOut = cms.untracked.int32(60),
                        connectionTimeOut = cms.untracked.int32(60),
                        connectionRetrialPeriod = cms.untracked.int32(10)
                      ),
                      BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
                      connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'),
                      globaltag = cms.string('GR_R_38X_V13::All') 
		      #globaltag = cms.string('GR_R_37X_V1::All')
                    )
process.options.wantSummary = False
		
process.load( "SimGeneral.HepPDTESSource.pythiapdt_cfi" )

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
          #fileNames = cms.untracked.vstring ("/store/data/Run2010B/MinimumBias/RECO/TRKFailureScenario_1_v1/0092/*.root"),
          fileNames = cms.untracked.vstring ("/store/data/Run2010B/MinimumBias/RECO/TRKFailureScenario_11_v1/0092/FAAB0A4A-73D3-DF11-B2D1-002618B27F8A.root"),
          secondaryFileNames = cms.untracked.vstring()
)

process.maxEvents.input = 1000
		
		
# Tell the process what filename to use to save the output
process.add_(cms.Service("TFileService",
              fileName = cms.string('TopSelection_FS11.root' ),
	      closeFileFast = cms.untracked.bool(False)  
	    ) 
)
		

#process.MessageLogger.categories +=  (['MC'])
#process.MessageLogger.categories +=  (['topchargeAnalyzerFullFW'])

############################
#
#   CONFIGURE PAT OBJECTS
#
############################

process.patMuons.usePV = False
process.patElectrons.usePV = False
process.patMuons.embedTrack =True

process.patJetFlavourAssociation.physicsDefinition= True

## Set the right sample for JEC factors here
process.patJetCorrFactors.corrSample      = 'Spring10'

from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
switchJetCollection( process,
                     jetCollection=cms.InputTag('ak5CaloJets'),
                     jetCorrLabel=('AK5', 'Calo'))

process.patJetGenJetMatch.matched="sisCone5GenJets"

############################
#
#   SELECTION V4
#
############################

process.load( "Top.Selection.topSelectionV4_cfi" )


############################
#
#   PATHS
#
############################

process.p = cms.Path(
      	process.patDefaultSequence * process.ttSemimuEventSelection 
)



############################
#
#   OUTPUT
#
############################

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('TopSelection_FS_11.root'),
                               # save only events passing the full path                               
		               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('keep *')                                
                               )
		
				
process.outpath = cms.EndPath(process.out)


############################
#
#   CONFIGURE PAT TRIGGER 
#
############################

process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
#process.patTriggerEvent.processName = 'REDIGI'
#process.patTrigger.processName = 'REDIGI'

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerAll( process, True )
