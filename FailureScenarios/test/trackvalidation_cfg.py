import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TRACKVALIDATION")

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    default = cms.untracked.PSet(
	limit = cms.untracked.int32(100)
    ),
    makeTriggerResults = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(False),
    Rethrow = cms.untracked.vstring(
	'Unknown', 
        'DictionaryNotFound', 
        'InsertFailure', 
        'Configuration', 
        'LogicError', 
        'UnimplementedFeature', 
        'InvalidReference', 
        'NullPointerError', 
        'NoProductSpecified', 
        'EventTimeout', 
        'EventCorruption', 
        'ModuleFailure', 
        'ScheduleExecutionFailure', 
        'EventProcessorFailure', 
        'FileInPathError',
        'FileReadError', 
        'FatalRootError', 
        'ReadBasketBuffers', 
    ),
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(),
	secondaryFileNames = cms.untracked.vstring(),
)

process.maxEvents = cms.untracked.PSet( 
	input = cms.untracked.int32(-1) 
)

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START44_V9::All'

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
process.load("MGeisler.PF_PU_AssoMap.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.load("MGeisler.PF_PU_AssoMap.PF_PU_AssoMap_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

########### configuration cutsRecoTrack ########
process.cutsRecoTracks.quality = cms.vstring('highPurity')
process.cutsRecoTracks.minRapidity = cms.double(-2.5)
process.cutsRecoTracks.maxRapidity = cms.double(2.5)
process.cutsRecoTracks.tip = cms.double(3.)
process.cutsRecoTracks.lip = cms.double(30.)
process.cutsRecoTracks.ptMin = cms.double(1.0)

########### configuration Tracks2Vertex ######## 
process.Tracks2Vertex.VertexAssOneDim = cms.untracked.bool(False)
process.Tracks2Vertex.VertexAssUseAbsDistance = cms.untracked.bool(False)
process.Tracks2Vertex.VertexCollection = cms.InputTag('offlinePrimaryVertices')

process.FirstVertexTrackCollection = cms.EDProducer('FirstVertexTracks',
          TrackCollection = cms.InputTag('cutsRecoTracks'),
          VertexCollection = cms.InputTag('offlinePrimaryVertices'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex'),
)	

process.Vtxtracks = cms.EDProducer('VertexTracks',
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  TrackCollection = cms.untracked.string('generalTracks'),
) 

########### configuration MultiTrackValidator ########
process.multiTrackValidator.outputFile = cms.string('generalTracks')
process.multiTrackValidator.associators = ['TrackAssociatorByHits']
process.multiTrackValidator.skipHistoFit=cms.untracked.bool(False)
process.multiTrackValidator.label = ['cutsRecoTracks','FirstVertexTrackCollection','Vtxtracks']
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.minpT = cms.double(0.1)
process.multiTrackValidator.maxpT = cms.double(3000.0)
process.multiTrackValidator.nintpT = cms.int32(40)
process.multiTrackValidator.UseAssociators = cms.bool(True)
process.multiTrackValidator.runStandalone = cms.bool(True)
	
#Additional settings
InFile=''
OutName=''
		
if len(sys.argv) > 2:
    InDir = 'file:/user/geisler/FS2011/RECO/Skim/'
    abbr= sys.argv[2]
    if abbr == "Draft":
        OutName="RECO-Tracks_Draft"
    if abbr == "FS01":
        OutName="RECO-Tracks_FS01"
    if abbr == "FS08":
        OutName="RECO-Tracks_FS08"
    if abbr == "FS09":
        OutName="RECO-Tracks_FS09"
    if abbr == "FS10":
        OutName="RECO-Tracks_FS10" 
    process.multiTrackValidator.outputFile = "../files/FS2012_" + OutName + "_TV.root" 
    print " output is ../files/FS2012_" + OutName + "_TV.root  \n"
    if len(sys.argv) > 3:
        process.maxEvents.input= int(sys.argv[3])
        if len(sys.argv) > 4:
            if sys.argv[4] == "True":
	        InDir = '/store/user/mgeisler/FailureScenarios/RECO/Skim/'

InFile = InDir+'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_'+abbr+'.root'

process.source.fileNames = cms.untracked.vstring(InFile)
print "Input file is " + InFile

# sequences, paths && schedule
process.T2V = cms.Sequence(
      process.Tracks2Vertex
)

process.FVTC = cms.Sequence(
      process.FirstVertexTrackCollection
)	

process.VT = cms.Sequence(
      process.Vtxtracks
)
		
process.validation = cms.Sequence(
    process.multiTrackValidator
)
	
process.p = cms.Path(
      process.cutsRecoTracks
    * process.T2V
    * process.FVTC
    * process.VT
    * process.validation
)
	
process.schedule = cms.Schedule(
      process.p
)
