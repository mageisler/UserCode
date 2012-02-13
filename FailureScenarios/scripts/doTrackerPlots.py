import os, fcntl, fcntl, select, sys, subprocess, math, array, getopt
import FWCore.ParameterSet.Config as cms
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
CreateTrackerPlots=False
UseDCache="False"

for o,p in opts: 
    if o in ['-n','--number']:
        NumberOfEvents = p
    if o in ['-r','--run']:
        if p=="True":
	    CreateTrackerPlots = True
        elif p=="False":
	    CreateTrackerPlots = False
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

if CreateTrackerPlots:
    print " Will create the TrackerPlots first. This may take some time!"


#_____________________
#
# CMSRUN
#_____________________

input='/user/geisler/FS2011/FileLists/InputFileList.txt'
dictOfLists = ReadProdDetails(input)

if CreateTrackerPlots:
    print "Number of Events per failure scenario: " + str(NumberOfEvents)	
    for i, scenario in enumerate(dictOfLists["abbreviation"]):
        if NumberOfEvents == "-1":
            cmsRun='cmsRun ../test/trackerplots_cfg.py ' + str(scenario) + ' ' + UseDCache + ' &>../logs/TP_log_'+scenario+'.txt'
        else:
            cmsRun='cmsRun ../test/trackerplots_cfg.py ' + str(scenario) + ' ' + str(NumberOfEvents) + ' ' + UseDCache + ' &>../logs/TP_log_'+scenario+'.txt'
        print cmsRun
        os.system(cmsRun)
	
scenarios = []

for i, scenario in enumerate(dictOfLists["abbreviation"]):
    if not str(scenario) == "Draft":
        print "\n Create fraction of valid hits for scenario " + str(scenario)
        #CreateFractionOfValidHits(scenario)
    print " Create Ks mass distribution for scenario " + str(scenario)
    scenarios.append(scenario)
    #CreateMassFit(scenarios,"Ks")
    CreateMassFit(scenarios,"Lambda")