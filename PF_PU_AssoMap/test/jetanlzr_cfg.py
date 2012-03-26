import FWCore.ParameterSet.Config as cms
import sys, os
	
datasample = ['QCD','DYToMuMu','GluGluToHToGG']
	
spectra = ['15to30','30to50','50to80','80to120','120to170','170to300','300to470','470to600','600to800','800to1000']
	
	
Outfile = "JetAnlzr"
TightSelection=True

argBegin = 2

if "jetanlzr_cfg.py" in str(sys.argv[0]):
    argBegin=1
		
for i in range(len(datasample)):
    if datasample[i] in str(sys.argv[argBegin]):
        Outfile+= "_" + str(datasample[i])	
		
for i in range(len(spectra)):
    if spectra[i] in str(sys.argv[argBegin]):
        Outfile+= "_Pt-" + str(spectra[i])	
	
if "False" in str(sys.argv[argBegin]):
    TightSelection = False
	    
Outfile += ".root"

print " Outfile set to " + Outfile

process = cms.Process("JetAnlzr")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/tmp/mgeisler/GluGluToHToGG_M-120_7TeV-pythia6.root'),
    fileNames = cms.untracked.vstring('file:/tmp/mgeisler/QCD_Pt-600to800.root'),
    #eventsToProcess = cms.untracked.VEventRange('1:93299'),
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)
 	
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)
		
### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_44_V13::All'

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

process.load("MGeisler.PF_PU_AssoMap.cuts_cff")
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
from MGeisler.TrackValidator.TrackingParticleSelection_cfi import *

########### track selection configuration ########
process.cutsRecoTracks.minRapidity = cms.double(-2.5)
process.cutsRecoTracks.maxRapidity = cms.double(2.5)	

### AssociationMap-specific includes		
from MGeisler.PF_PU_AssoMap.pf_pu_assomap_cfi import Tracks2Vertex
		
process.Tracks2Vertex1st = Tracks2Vertex.clone(
	AllSteps = cms.untracked.bool(False),
	VertexAssOneDim = cms.untracked.bool(False),
	VertexAssUseAbsDistance = cms.untracked.bool(True),
)
			 
process.Tracks2VertexAll = Tracks2Vertex.clone(
	AllSteps = cms.untracked.bool(True),
	VertexAssOneDim = cms.untracked.bool(False),
	VertexAssUseAbsDistance = cms.untracked.bool(True),
)		  
				  
process.FirstVertexTrackCollection1st = cms.EDProducer('FirstVertexTracks',
	  TrackCollection = cms.InputTag('cutsRecoTracks'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex1st'),
)		  
				       
process.FirstVertexTrackCollectionAll = cms.EDProducer('FirstVertexTracks',
	  TrackCollection = cms.InputTag('cutsRecoTracks'),
          VertexTrackAssociationMap = cms.InputTag('Tracks2VertexAll'),
)

process.PFCandidates1st = cms.EDProducer('PFCand_NoPU_WithAM',
	  PFCandidateCollection = cms.InputTag('particleFlow'),
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  VertexTrackAssociationMap = cms.InputTag('Tracks2Vertex1st'),
	  ConversionsCollection = cms.InputTag('allConversions'),
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	  IpMin = cms.double(0.3),
)

process.PFCandidatesAll1 = cms.EDProducer('PFCand_NoPU_WithAM',
	  PFCandidateCollection = cms.InputTag('particleFlow'),
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  VertexTrackAssociationMap = cms.InputTag('Tracks2VertexAll'),
	  ConversionsCollection = cms.InputTag('allConversions'),
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	  IpMin = cms.double(0.3),
)

process.PFCandidatesAll3 = cms.EDProducer('PFCand_NoPU_WithAM',
	  PFCandidateCollection = cms.InputTag('particleFlow'),
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  VertexTrackAssociationMap = cms.InputTag('Tracks2VertexAll'),
	  ConversionsCollection = cms.InputTag('allConversions'),
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	  IpMin = cms.double(0.6),
)

process.PFCandidatesAll5 = cms.EDProducer('PFCand_NoPU_WithAM',
	  PFCandidateCollection = cms.InputTag('particleFlow'),
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  VertexTrackAssociationMap = cms.InputTag('Tracks2VertexAll'),
	  ConversionsCollection = cms.InputTag('allConversions'),
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	  IpMin = cms.double(1.5),
)

process.PFCandidatesAll10 = cms.EDProducer('PFCand_NoPU_WithAM',
	  PFCandidateCollection = cms.InputTag('particleFlow'),
	  VertexCollection = cms.InputTag('offlinePrimaryVertices'),
	  VertexTrackAssociationMap = cms.InputTag('Tracks2VertexAll'),
	  ConversionsCollection = cms.InputTag('allConversions'),
	  V0KshortCollection = cms.InputTag('generalV0Candidates','Kshort'),
	  V0LambdaCollection = cms.InputTag('generalV0Candidates','Lambda'),
	  NIVertexCollection = cms.InputTag('particleFlowDisplacedVertex'),
	  IpMin = cms.double(3.),
)
				
process.trackValidator = cms.EDAnalyzer('TrackValidator',
	tcLabel = cms.VInputTag(cms.InputTag("cutsRecoTracks"),cms.InputTag("FirstVertexTrackCollection1st"),cms.InputTag("FirstVertexTrackCollectionAll")),
    	tcRefLabel = cms.InputTag("generalTracks"),
	pfLabel = cms.VInputTag(cms.InputTag("particleFlow"),cms.InputTag("PFCandidatesAll1"),cms.InputTag("PFCandidatesAll3"),cms.InputTag("PFCandidatesAll5"),cms.InputTag("PFCandidatesAll10"),cms.InputTag("PFCandidates1st")),
    	pfRefLabel = cms.InputTag("particleFlow"),
    	PULabel = cms.InputTag("addPileupInfo"),
    	TPLabel = cms.InputTag("mergedtruth","MergedTrackTruth"),
	ignoremissingtrackcollection=cms.bool(False),
	photonPtMin=cms.double(1.0),
	photonEtaMin=cms.double(-2.4),
	photonEtaMax=cms.double(2.4),
	photonLip=cms.double(30.),
	photonTip=cms.double(3.),
	UseLogPt=cms.bool(False),
	generalTpSelector = TrackingParticleSelectionGeneral,
)
	
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets	

process.kt4PFJets1st = kt4PFJets.clone(
	src = cms.InputTag("PFCandidates1st"),
	rParam = cms.double(0.4),
	doAreaFastjet = cms.bool(True),
	doRhoFastjet = cms.bool(True),
	Ghost_EtaMax = cms.double(6.5)
)	

process.kt4PFJetsAll1 = kt4PFJets.clone(
	src = cms.InputTag("PFCandidatesAll1"),
	rParam = cms.double(0.4),
	doAreaFastjet = cms.bool(True),
	doRhoFastjet = cms.bool(True),
	Ghost_EtaMax = cms.double(6.5)
)	

process.kt4PFJetsAll3 = kt4PFJets.clone(
	src = cms.InputTag("PFCandidatesAll3"),
	rParam = cms.double(0.4),
	doAreaFastjet = cms.bool(True),
	doRhoFastjet = cms.bool(True),
	Ghost_EtaMax = cms.double(6.5)
)	

process.kt4PFJetsAll5 = kt4PFJets.clone(
	src = cms.InputTag("PFCandidatesAll5"),
	rParam = cms.double(0.4),
	doAreaFastjet = cms.bool(True),
	doRhoFastjet = cms.bool(True),
	Ghost_EtaMax = cms.double(6.5)
)	

process.kt4PFJetsAll10 = kt4PFJets.clone(
	src = cms.InputTag("PFCandidatesAll10"),
	rParam = cms.double(0.4),
	doAreaFastjet = cms.bool(True),
	doRhoFastjet = cms.bool(True),
	Ghost_EtaMax = cms.double(6.5)
)

process.kt4GenJetsn = cms.EDProducer("CandViewNtpProducer",
	src = cms.InputTag("kt4GenJets"),
    	lazyParser = cms.untracked.bool(True),
    	prefix = cms.untracked.string(""),
    	eventInfo = cms.untracked.bool(True),
	variables = cms.VPSet(
		 cms.PSet(
		   tag = cms.untracked.string("pt"),
		   quantity = cms.untracked.string("pt")
		 ),
            	cms.PSet(
                   tag = cms.untracked.string("eta"),
                   quantity = cms.untracked.string("eta")
             	),
            	cms.PSet(
                   tag = cms.untracked.string("phi"),
                   quantity = cms.untracked.string("phi")
            	), 
	),     
)

process.kt4PFJetsAll1n = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJetsAll1"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.kt4PFJets1stn = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJets1st"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.kt4PFJetsAll3n = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJetsAll3"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.kt4PFJetsAll5n = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJetsAll5"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.kt4PFJetsAll10n = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJetsAll10"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.kt4PFJetsn = cms.EDProducer("CandViewNtpProducer", 
    src = cms.InputTag("kt4PFJets"),
    lazyParser = cms.untracked.bool(True),
    prefix = cms.untracked.string(""),
    eventInfo = cms.untracked.bool(True),
    variables = cms.VPSet(
        cms.PSet(
            tag = cms.untracked.string("pt"),
            quantity = cms.untracked.string("pt")
        ),
        cms.PSet(
            tag = cms.untracked.string("jetArea"),
            quantity = cms.untracked.string("jetArea")
        ),
        cms.PSet(
            tag = cms.untracked.string("eta"),
            quantity = cms.untracked.string("eta")
        ),
        cms.PSet(
            tag = cms.untracked.string("phi"),
            quantity = cms.untracked.string("phi")
        ),
    )
)

process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
    
process.jetanalyzer = cms.EDAnalyzer('JetAnlzr',
	genJets = cms.string("kt4GenJetsn"),
	recoJets = cms.vstring("kt4PFJets1stn","kt4PFJetsAll1n","kt4PFJetsAll3n","kt4PFJetsAll5n","kt4PFJetsAll10n","kt4PFJetsn"),
	PileUpInfo = cms.string("addPileupInfo"),
	JetCorrector = cms.string("kt4PFL1L2L3")
)


# paths	

if TightSelection:
    print " Tight Track selection is used"
    process.cutsRecoTracks.quality = cms.vstring('highPurity')
    process.cutsRecoTracks.tip = cms.double(3.)
    process.cutsRecoTracks.lip = cms.double(30.)
    process.cutsRecoTracks.ptMin = cms.double(1.)
else:
    print " Loose Track selection is used"
    process.cutsRecoTracks.quality = cms.vstring('highPurity','tight','loose')
    process.cutsRecoTracks.tip = cms.double(120.)
    process.cutsRecoTracks.lip = cms.double(280.)
    process.cutsRecoTracks.ptMin = cms.double(0.1)
    
    TrackingParticleSelectionGeneral.tip = cms.double(120.0)
    TrackingParticleSelectionGeneral.lip = cms.double(280.0)
    TrackingParticleSelectionGeneral.ptMin = cms.double(0.1)
    TrackingParticleSelectionGeneral.minHit = cms.int32(3)
    
    process.trackValidator.photonPtMin=cms.double(0.1)
    process.trackValidator.photonLip=cms.double(280.)
    process.trackValidator.photonTip=cms.double(120.)
    
### produce the association map
process.T2V = cms.Sequence(
      process.Tracks2Vertex1st
    + process.Tracks2VertexAll
)

### produce a collection of tracks associated to the first vertex
process.FVTC = cms.Sequence(
      process.FirstVertexTrackCollection1st
    + process.FirstVertexTrackCollectionAll
)

### produce a collection of PFCandidates associated to the first vertex
process.PFCand = cms.Sequence(
      process.PFCandidates1st
    + process.PFCandidatesAll1
    + process.PFCandidatesAll3
    + process.PFCandidatesAll5
    + process.PFCandidatesAll10
)

### do the efficiency and fake rate analyses for the charged particles
process.MTV = cms.Sequence(
      process.trackValidator
)

### produce a jet collection from the PFCandidates
process.PFJ = cms.Sequence(
      process.kt4PFJets1st
    + process.kt4PFJetsAll1
    + process.kt4PFJetsAll3
    + process.kt4PFJetsAll5
    + process.kt4PFJetsAll10
)

### produce the edm nTuples from the jet collection
process.nTs = cms.Sequence(
      process.kt4GenJetsn
    + process.kt4PFJets1stn
    + process.kt4PFJetsAll1n
    + process.kt4PFJetsAll3n
    + process.kt4PFJetsAll5n
    + process.kt4PFJetsAll10n
    + process.kt4PFJetsn
)
		
### do the jet analysis				  
process.JA = cms.Sequence(
     process.jetanalyzer
)

process.p = cms.Path(
      process.cutsRecoTracks    ### produce the track collection for the analysis of of the charged tracks
    * process.T2V		### produce the association map
    * process.FVTC		### produce a collection of tracks associated to the first vertex	
    * process.PFCand		### produce a collection of PFCandidates associated to the first vertex	
    * process.MTV		### do the efficiency and fake rate analyses for the particles 
    * process.PFJ		### produce a jet collection from the PFCandidates
    * process.nTs		### produce the edm nTuples from the jet collection
    * process.JA		### do the jet analysis	
)

process.schedule = cms.Schedule(
      process.p
)