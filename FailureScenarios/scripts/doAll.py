import os, fcntl, fcntl, select, sys, subprocess, math, array, getopt
from MGeisler.FailureScenarios.sourceFiles_cfi import *
from MGeisler.FailureScenarios.failureScenarioTools_cfi import *

print "\n This script will run all neccessary steps for the failure scenarios \n" 


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'n:s:p:v:r:d:'
keywords = ['number','skim','plots','validation','run','dcache']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

NumberOfEvents=-1
Skim=False
TrackerPlots=False
TrackValidation=False
UseDCache="False"
Run="False"

for o,p in opts: 
    if o in ['-n','--number']:
        NumberOfEvents = p
    if o in ['-s','--skim']:
        if p=="True":
	    Skim = True
        elif p=="False":
	    Skim = False
        else:
	    print "ERROR: Please set skim to True or False!"
    if o in ['-p','--plots']:
        if p=="True":
	    TrackerPlots = True
        elif p=="False":
	    TrackerPlots = False
        else:
	    print "ERROR: Please set plots to True or False!"
    if o in ['-v','--validation']:
        if p=="True":
	    TrackValidation = True
        elif p=="False":
	    TrackValidation = False
        else:
	    print "ERROR: Please set validation to True or False!"
    if o in ['-r','--run']:
        if p=="True":
	    Run = "True"
        elif p=="False":
	    Run = "False"
        else:
	    print "ERROR: Please set validation to True or False!"
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


print " Number of Events per failure scenario: " + str(NumberOfEvents) + " \n"

if Skim:
    print " Will skim the input files first. This may take some time!"

if TrackerPlots:
    print " Will create the TrackerPlots. This may take some time!"

if TrackValidation:
    print " Will do the TrackValidation. This may take some time!"

if UseDCache == "True":
    print " Will save the output file at /store/..."
else:
    print " Will save the output file at /user/..."
    
print " "


#_____________________
#
# CMSRUN
#_____________________

if Skim:	
    pyRun = 'python doSkim.py -n ' + str(NumberOfEvents)  + ' -d ' + UseDCache + ' &>../logs/SK_log.txt'
    print pyRun
    os.system(pyRun)

if TrackerPlots:	
    pyRun = 'python doTrackerPlots.py -n ' + str(NumberOfEvents) + ' -r ' + Run + ' -d ' + UseDCache + ' &>../logs/TP_log.txt'
    print pyRun
    os.system(pyRun)

if TrackValidation:	
    pyRun = 'python doTrackValidation.py -n ' + str(NumberOfEvents) + ' -r ' + Run + ' -d ' + UseDCache + ' &>../logs/TV_log.txt'
    print pyRun
    os.system(pyRun)