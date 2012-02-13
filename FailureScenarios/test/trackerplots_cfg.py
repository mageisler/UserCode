import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("TrackerPlots")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()		 
)
		
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.Trackerplots = cms.EDAnalyzer("TrackerPlots",
    FileName = cms.untracked.string("RECO-Tracks_08")	
)

process.Trackhitcount = cms.EDAnalyzer("TrackHitCount",
    FileName = cms.untracked.string("RECO-Tracks_08")	
)

process.Ksmass = cms.EDAnalyzer('KsMassFit',
    FileName = cms.untracked.string("RECO-Tracks_08")	
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('../files/FS2012.root')
)

		
#Additional settings
abbr=''
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
    if abbr == "All":
        OutName="RECO-Tracks_All"
    process.Trackerplots.FileName = OutName
    process.Trackhitcount.FileName = OutName   
    process.Ksmass.FileName = OutName
    process.TFileService.fileName = "../files/FS2012_" + OutName + ".root"   
    print "\n scenario set to " + abbr
    print " output is ../pictures/"+OutName+"_TrackerMap.png and ../files/FS2012_" + OutName + ".root \n"
    if len(sys.argv) > 3:
        process.maxEvents.input= int(sys.argv[3])
        if len(sys.argv) > 4:
            if sys.argv[4] == "True":
	        InDir = '/store/user/mgeisler/FailureScenarios/RECO/Skim/'

InFile = InDir+'QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_'+abbr+'.root'

process.source.fileNames = cms.untracked.vstring(InFile)
print "Input file is " + InFile

process.p = cms.Path(process.Trackerplots+process.Trackhitcount+process.Ksmass)
