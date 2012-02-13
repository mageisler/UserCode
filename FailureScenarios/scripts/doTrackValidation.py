import os, fcntl, fcntl, select, sys, subprocess, math, array, getopt
from MGeisler.TrackerPlots.sourceFiles_cfi import *
from MGeisler.TrackerPlots.failureScenarioTools_cfi import *

print "\n This script should create the trackerplots from the failure scenarios \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'n:r:d:'
keywords = ['number','run','dcache']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

NumberOfEvents=-1
RunTrackValidator=False
UseDCache="False"

for o,p in opts: 
    if o in ['-n','--number']:
        NumberOfEvents = p
    if o in ['-r','--run']:
        if p=="True":
	    RunTrackValidator = True
        elif p=="False":
	    RunTrackValidator = False
        else:
	    print "ERROR: Please set run to True or False!"
    if o in ['-d','--dcache']:
        if p=="True":
	    UseDCache = "True"
        elif p=="False":
	    UseDCache = "False"
        else:
	    print "ERROR: Please set dCache use to True or False!"


#_____________________
#
# Check Input Parameters
#_____________________

if RunTrackValidator:
    print " Will run the TrackValidator first. This may take some time!"


#_____________________
#
# CMSRUN
#_____________________

input='/user/geisler/FS2011/FileLists/InputFileList.txt'
dictOfLists = ReadProdDetails(input)

if RunTrackValidator:
    print "Number of Events per failure scenario: " + str(NumberOfEvents)	
    for i, scenario in enumerate(dictOfLists["abbreviation"]):
        if NumberOfEvents == "-1":
            cmsRun='cmsRun ../test/trackvalidation_cfg.py ' + str(scenario) + ' ' + UseDCache + ' &>../logs/TV_log_'+scenario+'.txt'
        else:
            cmsRun='cmsRun ../test/trackvalidation_cfg.py ' + str(scenario) + ' ' + str(NumberOfEvents) + ' ' + UseDCache + ' &>../logs/TV_log_'+scenario+'.txt'
        print cmsRun
        os.system(cmsRun)
	
for i, scenario in enumerate(dictOfLists["abbreviation"]):
    ref = True
    print "\n Create tracking efficiency plots for scenario " + str(scenario)
    if str(scenario) == "Draft":
        ref = False
    CreateTrackingEfficiencyPlots(scenario,ref)