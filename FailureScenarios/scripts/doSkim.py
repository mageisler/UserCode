import os, fcntl, fcntl, select, sys, subprocess, math, array, getopt
from MGeisler.TrackerPlots.sourceFiles_cfi import *

print "\n This script should skim the failure scenario files to a usable size \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'n:d:'
keywords = ['number','dcache']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

NumberOfEvents=1
UseDCache="False"

for o,p in opts: 
    if o in ['-n','--number']:
        NumberOfEvents = p
    if o in ['-d','--dcache']:
        if p=="True":
	    UseDCache = "True"
        elif p=="False":
	    UseDCache = "False"
        else:
	    print "ERROR: Please set dCache use to True or False!"
	
if NumberOfEvents=="-1":
    NumberOfEvents=30000

print "Number of Events per failure scenario: " + str(NumberOfEvents)


#_____________________
#
# CMSRUN
#_____________________

input='/user/geisler/FS2011/FileLists/InputFileList.txt'
dictOfLists = ReadProdDetails(input)

for i, scenario in enumerate(dictOfLists["abbreviation"]):
    cmsRun='cmsRun ../test/skim_cfg.py ' + str(scenario) + ' ' + str(NumberOfEvents) + ' ' + UseDCache + ' &>../logs/SK_log_'+scenario+'.txt'
    print cmsRun
    os.system(cmsRun)
    if UseDCache == "True":
	FileName='QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_PU_S6-START44_SKIM_'+str(scenario)+'.root' 
	tmpdir='/user/geisler/FS2011/RECO/Tmp/'
	dcdir='/grid-srm:8443/pnfs/physik.rwth-aachen.de/cms/store/user/mgeisler/FailureScenarios/RECO/Skim/'
        dcrmComm = 'srmrm srm:/' + dcdir + FileName
	print dcrmComm
        os.system(dcrmComm)
        cpComm = 'srmcp file:///' + tmpdir + FileName + ' srm:/' + dcdir + FileName
	print cpComm
        os.system(cpComm)
	rmComm = 'rm ' + tmpdir + FileName
	print rmComm
        os.system(rmComm)

print ""