
## set the outname of the analyzer

import sys

datasample = ['QCD','DYToMuMu','GluGluToHToGG']
	
spectra = ['15to30','30to50','50to80','80to120','120to170','170to300','300to470','470to600','600to800','800to1000']


OutfileAnlzr = "JetAnlzr"
    	
if len(sys.argv)>2:
		
    argBegin = 2

    if "jetanlzr" in str(sys.argv[0]):
        argBegin=1	
    for i in range(len(datasample)):
        if datasample[i] in str(sys.argv[argBegin]):
            OutfileAnlzr+= "_" + str(datasample[i])	
    for i in range(len(spectra)):
        if spectra[i] in str(sys.argv[argBegin]):
            OutfileAnlzr+= "_Pt-" + str(spectra[i])

Outfile= OutfileAnlzr+".root"
	

print " Outfile set to " + Outfile

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *
	
getattr(process,"outpath").remove(
    getattr(process,"out")
)

## Options and Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.source.fileNames = cms.untracked.vstring('file:/user/geisler/QCD_Pt-600to800_TuneZ2star_8TeV_PU_S7_START52_V9-v1_GEN-SIM-RECODEBUG.root')
		
process.maxEvents.input = cms.untracked.int32(5)
	
process.GlobalTag.globaltag = 'START52_V11C::All'
 	
process.TFileService = cms.Service('TFileService',
    fileName = cms.string(Outfile),
    closeFileFast = cms.untracked.bool(True),
)
		
# Next Configure PAT to use PF2PAT WITH CHS as reference
from PhysicsTools.PatAlgos.tools.pfTools import *
		
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix, jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']), pvCollection=cms.InputTag('selectedPrimaryVertexQuality'))
	
### includes
	
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

### validation-specific includes
process.load("Validation.Configuration.postValidation_cff")
process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
		
process.perfectTracks = cms.EDProducer('PerfectAssociation',
	PFCandidateCollection = cms.InputTag('particleFlow'),
	trackCollection = cms.InputTag('generalTracks'),
	genParticles = cms.InputTag('genParticles'),
    	TrackingParticles = cms.InputTag("mergedtruth","MergedTrackTruth"),
)	
	
process.selectedPrimaryVertexQuality = cms.EDFilter("VertexSelector",
   	src = cms.InputTag('offlinePrimaryVertices'),
	cut = cms.string("isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2."),
	filter = cms.bool(False),
)	
				    		
### PFCandidate AssociationMap-specific includes
from CommonTools.RecoUtils.pfcand_assomap_cfi import PFCandAssoMap
		
process.PFCand2VertexAM = PFCandAssoMap.clone(
          VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
	  UseBeamSpotCompatibility = cms.untracked.bool(False),
)
		
process.PFCand2VertexAM2 = PFCandAssoMap.clone(
          VertexCollection = cms.InputTag('selectedPrimaryVertexQuality'),
	  UseBeamSpotCompatibility = cms.untracked.bool(True),
)
		
### PFCandidateCollection-specific includes
from CommonTools.RecoUtils.pfcand_nopu_witham_cfi import FirstVertexPFCandidates
		
process.PFCand = FirstVertexPFCandidates.clone(
          VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAM'),
)
		
process.PFCand2 = FirstVertexPFCandidates.clone(
          VertexPFCandAssociationMap = cms.InputTag('PFCand2VertexAM2'),
)
	
### JetProducer-specific includes
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets	

process.ak5PFJetsNew = ak5PFJets.clone(
	src = cms.InputTag("PFCand")
)

process.ak5PFJetsNew2 = ak5PFJets.clone(
	src = cms.InputTag("PFCand2")
)

process.perfectJets = ak5PFJets.clone(
	src = cms.InputTag("perfectTracks")
)
		

########################## PF2PAT ############################
	
process.pfPileUpPFlow.Enable = True
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpPFlow.Vertices = cms.InputTag('selectedPrimaryVertexQuality')
process.pfJetsPFlow.doAreaFastjet = True
process.pfJetsPFlow.doRhoFastjet = False

process.ak5PFJetsPFlow = ak5PFJets.clone(
    rParam = cms.double(0.5),
    src = cms.InputTag("pfNoElectron"+postfix),
    doAreaFastjet = cms.bool(False),
    doRhoFastjet = cms.bool(False)
)


getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.ak5PFJetsPFlow 
)
		
########################## PF2PAT END ############################
				
	
### paths & sequences
		
##sequence to produce basic stuff
process.priorstuff = cms.Sequence(
	process.perfectTracks *
	process.selectedPrimaryVertexQuality 
)
		
##sequence up to the producer of the pf candidates
process.assomap = cms.Sequence(
	process.PFCand2VertexAM *
	process.PFCand2VertexAM2 *
	process.PFCand *
	process.PFCand2 
)
		
##sequence to produce  the jet collection
process.jets = cms.Sequence(
	process.ak5PFJetsNew *
	process.ak5PFJetsNew2 *
	process.perfectJets  
)
	
process.help = process.jets.copy()
		
##sequence to produce the chs jet as reference
from MGeisler.JetAnlzr.cfg_tools import RemoveModules
process.pf2pat = RemoveModules(process,getattr(process,"patPF2PATSequence"+postfix),"ak5PFJetsPFlow"+postfix)

##sequence to produce the jet corrections
corrections = ["L123Residual","L23Residual"]
from MGeisler.JetAnlzr.corrections_cfi import GetCorrSequence
process.help.__iadd__(process.ak5PFJetsPFlow)
process.corr = GetCorrSequence(process,process.help, corrections)

##sequence to produce the nTuple for the jet collections
from MGeisler.JetAnlzr.nTuple_cfi import GetNTupsSequence
process.help.__iadd__(process.corr)
process.ntup = GetNTupsSequence(process,process.help)

##sequence for the jet analyzers
from MGeisler.JetAnlzr.jAnlzr_cfi import GetJAnlzrSequence
process.janlzr = GetJAnlzrSequence(process,process.help, corrections)
  
process.p = cms.Path( 
	  process.priorstuff *
	  process.assomap *
	  process.jets *
	  process.pf2pat *
	  process.corr *
	  process.ntup *
	  process.janlzr
)
