#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT, os
from MGeisler.JetAnlzr.AnalysisTools import *


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:'
keywords = ['path']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)



path=""

for o,p in opts: 
    if o in ['-p','--path']:
        path = p


filelist = GetListOfFiles(path)
file_number = len(filelist)
#file_number = 1

if file_number<1:
    print "Too less files"
    sys.exit()
    
npu_number = len(npu_ranges)

eta_number = len(eta_ranges)

pt_number = len(pt_ranges)

#_____________________
#
# READ INPUT FILES
#_____________________

liste = ""
paras = ()

for file_ite in range(file_number):

    inputfile = file(filelist[file_ite])
    
    print "File " +str(file_ite+1)+ ": " + str(filelist[file_ite]) + " ...",
    
    for line in inputfile.readlines():  	     
        liste += line
	splittedline = line.split("\t")
	lastelement = len(splittedline)
	   
	if lastelement==4:
	    paras+=[float(splittedline[0]),float(splittedline[1]),float(splittedline[2]),float(splittedline[3])],

    print "done" 

    
## order is genPt - recoPt - genEta - npu 
    
for npu_ite in range(npu_number):
		    
    if (npu_ite==(pt_number-1)):
	npuMin = 0
	npuMax = 60
    else:
	npuMin = npu_ranges[npu_ite]
	npuMax = npu_ranges[npu_ite+1]
	
    print "Npu range from " + str(npuMin) + " to " +str(npuMax)
	
    for pt_ite in range(pt_number):
		    
        if (pt_ite==(pt_number-1)):
	    ptMin = 0.
	    ptMax = pt_ranges[pt_ite]
        else:
	    ptMin = pt_ranges[pt_ite]
	    ptMax = pt_ranges[pt_ite+1]
	
        print "GenPt range from " + str(ptMin) + " to " +str(ptMax)
     
        for eta_ite in range(eta_number):
		    
            if (eta_ite==(eta_number-1)):
 	        etaMin = -5.
		etaMax = 5.
	    else:		    
		etaMin = eta_ranges[eta_ite]
		etaMax = eta_ranges[eta_ite+1]
	
            print "GenEta range from " + str(etaMin) + " to " +str(etaMax)
		
	    sum_help = 0.
	    num_help = 0
		
	    for para_ite in range(len(paras)):
	        if ((paras[para_ite][0]>=ptMin) and (paras[para_ite][0]<=ptMax)):
	            if ((paras[para_ite][2]>=etaMin) and (paras[para_ite][2]<=etaMax)):
	                if ((paras[para_ite][3]>=npuMin) and (paras[para_ite][3]<=npuMax)):
			    sum_help+= paras[para_ite][1]/paras[para_ite][0]
			    num_help+= 1
		
	    if num_help==0:
		result = 0.
	    else:
	        result = sum_help*1./num_help
				
	    print result
	    print ""
	