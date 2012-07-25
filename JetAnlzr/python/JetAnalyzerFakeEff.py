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

if file_number<1:
    print "Too less files"
    sys.exit()
    
name = ''
    
for dataset in possibleDatasets:
    if dataset in filelist[0]:
	name = dataset
    

#_____________________
#
# READ INPUT FILES
#_____________________

dirnames = GetListOfSubdirectories(filelist[0],"","jetanalyzer")
dir_number = len(dirnames)
dir_number = 1

subdirnames = ()

for i in range(dir_number):
	
    subdirnames+=GetListOfSubdirectories(filelist[0],dirnames[i],"ak5PFJets"),
    
ass_number = len(subdirnames[0])

#print "Concentrating on "
    
#for ass_ite in range(2,ass_number):
    #for dir_ite in range(dir_number):
        #print subdirnames[dir_ite][ass_ite]
	
print ""

File_ref = [[]] * file_number

histo_sim_jet_PtcR_NpucR_EtacR_queue = ()
histo_assoc_PtcR_NpucR_EtacR_queue = ()

histo_reco_jet_PtcR_NpucR_EtacR_queue = ()
histo_assoc2_PtcR_NpucR_EtacR_queue = ()


histo_sim_jet_PtcR_Npu15_EtacR_queue = ()
histo_assoc_PtcR_Npu15_EtacR_queue = ()

histo_reco_jet_PtcR_Npu15_EtacR_queue = ()
histo_assoc2_PtcR_Npu15_EtacR_queue = ()


histo_sim_jet_PtcR_Npu30_EtacR_queue = ()
histo_assoc_PtcR_Npu30_EtacR_queue = ()

histo_reco_jet_PtcR_Npu30_EtacR_queue = ()
histo_assoc2_PtcR_Npu30_EtacR_queue = ()


histo_sim_jet_PtcR_Npu50_EtacR_queue = ()
histo_assoc_PtcR_Npu50_EtacR_queue = ()

histo_reco_jet_PtcR_Npu50_EtacR_queue = ()
histo_assoc2_PtcR_Npu50_EtacR_queue = ()

for file_ite in range(file_number):

    File_ref[file_ite] = ROOT.TFile.Open(filelist[file_ite])
    
    print "File " + str(file_ite+1)+ ": " + str(filelist[file_ite]) + " ...",

    histo_sim_jet_PtcR_NpucR_EtacR_ass = ()
    histo_assoc_PtcR_NpucR_EtacR_ass = ()

    histo_reco_jet_PtcR_NpucR_EtacR_ass = ()
    histo_assoc2_PtcR_NpucR_EtacR_ass = ()
    

    histo_sim_jet_PtcR_Npu15_EtacR_ass = ()
    histo_assoc_PtcR_Npu15_EtacR_ass = ()

    histo_reco_jet_PtcR_Npu15_EtacR_ass = ()
    histo_assoc2_PtcR_Npu15_EtacR_ass = ()
    

    histo_sim_jet_PtcR_Npu30_EtacR_ass = ()
    histo_assoc_PtcR_Npu30_EtacR_ass = ()

    histo_reco_jet_PtcR_Npu30_EtacR_ass = ()
    histo_assoc2_PtcR_Npu30_EtacR_ass = ()
    

    histo_sim_jet_PtcR_Npu50_EtacR_ass = ()
    histo_assoc_PtcR_Npu50_EtacR_ass = ()

    histo_reco_jet_PtcR_Npu50_EtacR_ass = ()
    histo_assoc2_PtcR_Npu50_EtacR_ass = ()
    
    for ass_ite in range(ass_number):

        histo_sim_jet_PtcR_NpucR_EtacR_ass_dir = () 
        histo_assoc_PtcR_NpucR_EtacR_ass_dir = ()

        histo_reco_jet_PtcR_NpucR_EtacR_ass_dir = ()
        histo_assoc2_PtcR_NpucR_EtacR_ass_dir = () 
	

        histo_sim_jet_PtcR_Npu15_EtacR_ass_dir = () 
        histo_assoc_PtcR_Npu15_EtacR_ass_dir = ()

        histo_reco_jet_PtcR_Npu15_EtacR_ass_dir = ()
        histo_assoc2_PtcR_Npu15_EtacR_ass_dir = ()  
	

        histo_sim_jet_PtcR_Npu30_EtacR_ass_dir = () 
        histo_assoc_PtcR_Npu30_EtacR_ass_dir = ()

        histo_reco_jet_PtcR_Npu30_EtacR_ass_dir = ()
        histo_assoc2_PtcR_Npu30_EtacR_ass_dir = () 
	

        histo_sim_jet_PtcR_Npu50_EtacR_ass_dir = () 
        histo_assoc_PtcR_Npu50_EtacR_ass_dir = ()

        histo_reco_jet_PtcR_Npu50_EtacR_ass_dir = ()
        histo_assoc2_PtcR_Npu50_EtacR_ass_dir = () 
	
	for dir_ite in range(dir_number):
		
	    help = random.random()
	    
	    sim_jet_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+"SimJetsVsPtVsEtaVsNpu")
	    assoc_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+"AssocJetsVsPtVsEtaVsNpu")
	    
	    reco_jet_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+"RecoJetsVsPtVsEtaVsNpu")
	    assoc2_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+"Assoc2JetsVsPtVsEtaVsNpu")
	    
	    histo_sim_jet_PtcR_NpucR_EtacR_ass_dir+= sim_jet_help.Project3D(str(help)+"xyo"),
	    histo_assoc_PtcR_NpucR_EtacR_ass_dir+= assoc_help.Project3D(str(help)+"xyo"),
	    
	    histo_reco_jet_PtcR_NpucR_EtacR_ass_dir+= reco_jet_help.Project3D(str(help)+"xyo"),
	    histo_assoc2_PtcR_NpucR_EtacR_ass_dir+= assoc2_help.Project3D(str(help)+"xyo"),
	    
	    ##Npu 15
	    sim_jet_help.GetZaxis().SetRange(1,16)
	    assoc_help.GetZaxis().SetRange(1,16)
	    histo_sim_jet_PtcR_Npu15_EtacR_ass_dir+= sim_jet_help.Project3D(str(help+1)+"xyo"),
	    histo_assoc_PtcR_Npu15_EtacR_ass_dir+= assoc_help.Project3D(str(help+1)+"xyo"),
	    
	    reco_jet_help.GetZaxis().SetRange(1,16)
	    assoc2_help.GetZaxis().SetRange(1,16)
	    histo_reco_jet_PtcR_Npu15_EtacR_ass_dir+= reco_jet_help.Project3D(str(help+1)+"xyo"),
	    histo_assoc2_PtcR_Npu15_EtacR_ass_dir+= assoc2_help.Project3D(str(help+1)+"xyo"),
	    
	    ##Npu 30
	    sim_jet_help.GetZaxis().SetRange(17,31)
	    assoc_help.GetZaxis().SetRange(17,31)
	    histo_sim_jet_PtcR_Npu30_EtacR_ass_dir+= sim_jet_help.Project3D(str(help+2)+"xyo"),
	    histo_assoc_PtcR_Npu30_EtacR_ass_dir+= assoc_help.Project3D(str(help+2)+"xyo"),
	    
	    reco_jet_help.GetZaxis().SetRange(17,31)
	    assoc2_help.GetZaxis().SetRange(17,31)
	    histo_reco_jet_PtcR_Npu30_EtacR_ass_dir+= reco_jet_help.Project3D(str(help+2)+"xyo"),
	    histo_assoc2_PtcR_Npu30_EtacR_ass_dir+= assoc2_help.Project3D(str(help+2)+"xyo"),
	    
	    ##Npu 50
	    sim_jet_help.GetZaxis().SetRange(32,50)
	    assoc_help.GetZaxis().SetRange(32,50)
	    histo_sim_jet_PtcR_Npu50_EtacR_ass_dir+= sim_jet_help.Project3D(str(help+3)+"xyo"),
	    histo_assoc_PtcR_Npu50_EtacR_ass_dir+= assoc_help.Project3D(str(help+3)+"xyo"),
	    
	    reco_jet_help.GetZaxis().SetRange(32,50)
	    assoc2_help.GetZaxis().SetRange(32,50)
	    histo_reco_jet_PtcR_Npu50_EtacR_ass_dir+= reco_jet_help.Project3D(str(help+3)+"xyo"),
	    histo_assoc2_PtcR_Npu50_EtacR_ass_dir+= assoc2_help.Project3D(str(help+3)+"xyo"),
	            
	    
        histo_sim_jet_PtcR_NpucR_EtacR_ass+= histo_sim_jet_PtcR_NpucR_EtacR_ass_dir,
        histo_assoc_PtcR_NpucR_EtacR_ass+= histo_assoc_PtcR_NpucR_EtacR_ass_dir,
	
        histo_reco_jet_PtcR_NpucR_EtacR_ass+= histo_reco_jet_PtcR_NpucR_EtacR_ass_dir,
        histo_assoc2_PtcR_NpucR_EtacR_ass+= histo_assoc2_PtcR_NpucR_EtacR_ass_dir,
	            
	    
        histo_sim_jet_PtcR_Npu15_EtacR_ass+= histo_sim_jet_PtcR_Npu15_EtacR_ass_dir,
        histo_assoc_PtcR_Npu15_EtacR_ass+= histo_assoc_PtcR_Npu15_EtacR_ass_dir,
	
        histo_reco_jet_PtcR_Npu15_EtacR_ass+= histo_reco_jet_PtcR_Npu15_EtacR_ass_dir,
        histo_assoc2_PtcR_Npu15_EtacR_ass+= histo_assoc2_PtcR_Npu15_EtacR_ass_dir,
	            
	    
        histo_sim_jet_PtcR_Npu30_EtacR_ass+= histo_sim_jet_PtcR_Npu30_EtacR_ass_dir,
        histo_assoc_PtcR_Npu30_EtacR_ass+= histo_assoc_PtcR_Npu30_EtacR_ass_dir,
	
        histo_reco_jet_PtcR_Npu30_EtacR_ass+= histo_reco_jet_PtcR_Npu30_EtacR_ass_dir,
        histo_assoc2_PtcR_Npu30_EtacR_ass+= histo_assoc2_PtcR_Npu30_EtacR_ass_dir,
	            
	    
        histo_sim_jet_PtcR_Npu50_EtacR_ass+= histo_sim_jet_PtcR_Npu50_EtacR_ass_dir,
        histo_assoc_PtcR_Npu50_EtacR_ass+= histo_assoc_PtcR_Npu50_EtacR_ass_dir,
	
        histo_reco_jet_PtcR_Npu50_EtacR_ass+= histo_reco_jet_PtcR_Npu50_EtacR_ass_dir,
        histo_assoc2_PtcR_Npu50_EtacR_ass+= histo_assoc2_PtcR_Npu50_EtacR_ass_dir,
	
	
    histo_sim_jet_PtcR_NpucR_EtacR_queue+= histo_sim_jet_PtcR_NpucR_EtacR_ass,
    histo_assoc_PtcR_NpucR_EtacR_queue+= histo_assoc_PtcR_NpucR_EtacR_ass,	
	
    histo_reco_jet_PtcR_NpucR_EtacR_queue+= histo_reco_jet_PtcR_NpucR_EtacR_ass,	
    histo_assoc2_PtcR_NpucR_EtacR_queue+= histo_assoc2_PtcR_NpucR_EtacR_ass,
	
	
    histo_sim_jet_PtcR_Npu15_EtacR_queue+= histo_sim_jet_PtcR_Npu15_EtacR_ass,
    histo_assoc_PtcR_Npu15_EtacR_queue+= histo_assoc_PtcR_Npu15_EtacR_ass,	
	
    histo_reco_jet_PtcR_Npu15_EtacR_queue+= histo_reco_jet_PtcR_Npu15_EtacR_ass,	
    histo_assoc2_PtcR_Npu15_EtacR_queue+= histo_assoc2_PtcR_Npu15_EtacR_ass,
	
	
    histo_sim_jet_PtcR_Npu30_EtacR_queue+= histo_sim_jet_PtcR_Npu30_EtacR_ass,
    histo_assoc_PtcR_Npu30_EtacR_queue+= histo_assoc_PtcR_Npu30_EtacR_ass,	
	
    histo_reco_jet_PtcR_Npu30_EtacR_queue+= histo_reco_jet_PtcR_Npu30_EtacR_ass,	
    histo_assoc2_PtcR_Npu30_EtacR_queue+= histo_assoc2_PtcR_Npu30_EtacR_ass,
	
	
    histo_sim_jet_PtcR_Npu50_EtacR_queue+= histo_sim_jet_PtcR_Npu50_EtacR_ass,
    histo_assoc_PtcR_Npu50_EtacR_queue+= histo_assoc_PtcR_Npu50_EtacR_ass,	
	
    histo_reco_jet_PtcR_Npu50_EtacR_queue+= histo_reco_jet_PtcR_Npu50_EtacR_ass,	
    histo_assoc2_PtcR_Npu50_EtacR_queue+= histo_assoc2_PtcR_Npu50_EtacR_ass,

    print "done"
    
print ""
	    
histo_sim_jet_PtcR_NpucR_EtacR_com = GetCombinedTH2FsWeight2(histo_sim_jet_PtcR_NpucR_EtacR_queue,filelist)
histo_assoc_PtcR_NpucR_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc_PtcR_NpucR_EtacR_queue,filelist)

histo_reco_jet_PtcR_NpucR_EtacR_com = GetCombinedTH2FsWeight2(histo_reco_jet_PtcR_NpucR_EtacR_queue,filelist)
histo_assoc2_PtcR_NpucR_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc2_PtcR_NpucR_EtacR_queue,filelist)
	    
##Npu 15	    
histo_sim_jet_PtcR_Npu15_EtacR_com = GetCombinedTH2FsWeight2(histo_sim_jet_PtcR_Npu15_EtacR_queue,filelist)
histo_assoc_PtcR_Npu15_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc_PtcR_Npu15_EtacR_queue,filelist)

histo_reco_jet_PtcR_Npu15_EtacR_com = GetCombinedTH2FsWeight2(histo_reco_jet_PtcR_Npu15_EtacR_queue,filelist)
histo_assoc2_PtcR_Npu15_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc2_PtcR_Npu15_EtacR_queue,filelist)
	    
##Npu 30	    
histo_sim_jet_PtcR_Npu30_EtacR_com = GetCombinedTH2FsWeight2(histo_sim_jet_PtcR_Npu30_EtacR_queue,filelist)
histo_assoc_PtcR_Npu30_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc_PtcR_Npu30_EtacR_queue,filelist)

histo_reco_jet_PtcR_Npu30_EtacR_com = GetCombinedTH2FsWeight2(histo_reco_jet_PtcR_Npu30_EtacR_queue,filelist)
histo_assoc2_PtcR_Npu30_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc2_PtcR_Npu30_EtacR_queue,filelist)
	    
##Npu 50	    
histo_sim_jet_PtcR_Npu50_EtacR_com = GetCombinedTH2FsWeight2(histo_sim_jet_PtcR_Npu50_EtacR_queue,filelist)
histo_assoc_PtcR_Npu50_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc_PtcR_Npu50_EtacR_queue,filelist)

histo_reco_jet_PtcR_Npu50_EtacR_com = GetCombinedTH2FsWeight2(histo_reco_jet_PtcR_Npu50_EtacR_queue,filelist)
histo_assoc2_PtcR_Npu50_EtacR_com = GetCombinedTH2FsWeight2(histo_assoc2_PtcR_Npu50_EtacR_queue,filelist)


## print the overall fractions
for ass_ite in range(ass_number):
	
    histo_effic_num_entries = histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][0].Integral() 
    histo_effic_denum_entries = histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][0].Integral()
	
    histo_fakerate_num_entries = histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][0].Integral() 
    histo_fakerate_denum_entries = histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][0].Integral()
      
    effic = histo_effic_num_entries *1./histo_effic_denum_entries
    fakerate = 1.- histo_fakerate_num_entries *1./histo_fakerate_denum_entries
    
    frac = effic *1./fakerate
      
    print ""  
    print "Version: " + subdirnames[0][ass_ite]
    print "Effic: " + str(effic)
    print "Fake rate: " + str(fakerate)
    print "Fraction: " + str(frac)
    print ""
    
print ""    
print "Creating 2D ratio"
    
effic_PtcR_NpucR_EtacR = ()
fakerate_PtcR_NpucR_EtacR = ()
    
effic_PtcR_Npu15_EtacR = ()
fakerate_PtcR_Npu15_EtacR = ()
    
effic_PtcR_Npu30_EtacR = ()
fakerate_PtcR_Npu30_EtacR = ()
    
effic_PtcR_Npu50_EtacR = ()
fakerate_PtcR_Npu50_EtacR = ()
    
for ass_ite in range(ass_number):
    
    effic_PtcR_NpucR_EtacR+= GetFractionTH2Fs(histo_assoc_PtcR_NpucR_EtacR_com[ass_ite],histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite],"effic"),
    fakerate_PtcR_NpucR_EtacR+=  GetFractionTH2Fs(histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite],histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite],"fakerate"),
    
    effic_PtcR_Npu15_EtacR+= GetFractionTH2Fs(histo_assoc_PtcR_Npu15_EtacR_com[ass_ite],histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite],"effic"),
    fakerate_PtcR_Npu15_EtacR+=  GetFractionTH2Fs(histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite],histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite],"fakerate"),
    
    effic_PtcR_Npu30_EtacR+= GetFractionTH2Fs(histo_assoc_PtcR_Npu30_EtacR_com[ass_ite],histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite],"effic"),
    fakerate_PtcR_Npu30_EtacR+=  GetFractionTH2Fs(histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite],histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite],"fakerate"),
    
    effic_PtcR_Npu50_EtacR+= GetFractionTH2Fs(histo_assoc_PtcR_Npu50_EtacR_com[ass_ite],histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite],"effic"),
    fakerate_PtcR_Npu50_EtacR+=  GetFractionTH2Fs(histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite],histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite],"fakerate"),
    
    
histo_assoc_NpucR_EtacR_Pt = ()
histo_sim_jet_NpucR_EtacR_Pt = ()
    
histo_assoc_NpucR_PtcR_Eta = ()
histo_sim_jet_NpucR_PtcR_Eta = ()   
    
histo_assoc2_NpucR_EtacR_Pt = ()
histo_reco_jet_NpucR_EtacR_Pt = ()
    
histo_assoc2_NpucR_PtcR_Eta = ()
histo_reco_jet_NpucR_PtcR_Eta = ()
    
##Npu 15   
histo_assoc_Npu15_EtacR_Pt = ()
histo_sim_jet_Npu15_EtacR_Pt = ()
    
histo_assoc_Npu15_PtcR_Eta = ()
histo_sim_jet_Npu15_PtcR_Eta = ()   
    
histo_assoc2_Npu15_EtacR_Pt = ()
histo_reco_jet_Npu15_EtacR_Pt = ()
    
histo_assoc2_Npu15_PtcR_Eta = ()
histo_reco_jet_Npu15_PtcR_Eta = ()
    
##Npu 30  
histo_assoc_Npu30_EtacR_Pt = ()
histo_sim_jet_Npu30_EtacR_Pt = ()
    
histo_assoc_Npu30_PtcR_Eta = ()
histo_sim_jet_Npu30_PtcR_Eta = ()   
    
histo_assoc2_Npu30_EtacR_Pt = ()
histo_reco_jet_Npu30_EtacR_Pt = ()
    
histo_assoc2_Npu30_PtcR_Eta = ()
histo_reco_jet_Npu30_PtcR_Eta = ()
    
##Npu 50   
histo_assoc_Npu50_EtacR_Pt = ()
histo_sim_jet_Npu50_EtacR_Pt = ()
    
histo_assoc_Npu50_PtcR_Eta = ()
histo_sim_jet_Npu50_PtcR_Eta = ()   
    
histo_assoc2_Npu50_EtacR_Pt = ()
histo_reco_jet_Npu50_EtacR_Pt = ()
    
histo_assoc2_Npu50_PtcR_Eta = ()
histo_reco_jet_Npu50_PtcR_Eta = ()
    
    
for ass_ite in range(ass_number):
    
    histo_assoc_NpucR_EtacR_Pt_ass = ()
    histo_sim_jet_NpucR_EtacR_Pt_ass = ()    
    
    histo_assoc_NpucR_PtcR_Eta_ass = ()
    histo_sim_jet_NpucR_PtcR_Eta_ass = ()    
    
    histo_assoc2_NpucR_EtacR_Pt_ass = ()
    histo_reco_jet_NpucR_EtacR_Pt_ass = ()    
    
    histo_assoc2_NpucR_PtcR_Eta_ass = ()
    histo_reco_jet_NpucR_PtcR_Eta_ass = ()
    
    
    histo_assoc_Npu15_EtacR_Pt_ass = ()
    histo_sim_jet_Npu15_EtacR_Pt_ass = ()    
    
    histo_assoc_Npu15_PtcR_Eta_ass = ()
    histo_sim_jet_Npu15_PtcR_Eta_ass = ()    
    
    histo_assoc2_Npu15_EtacR_Pt_ass = ()
    histo_reco_jet_Npu15_EtacR_Pt_ass = ()    
    
    histo_assoc2_Npu15_PtcR_Eta_ass = ()
    histo_reco_jet_Npu15_PtcR_Eta_ass = ()
    
    
    histo_assoc_Npu30_EtacR_Pt_ass = ()
    histo_sim_jet_Npu30_EtacR_Pt_ass = ()    
    
    histo_assoc_Npu30_PtcR_Eta_ass = ()
    histo_sim_jet_Npu30_PtcR_Eta_ass = ()    
    
    histo_assoc2_Npu30_EtacR_Pt_ass = ()
    histo_reco_jet_Npu30_EtacR_Pt_ass = ()    
    
    histo_assoc2_Npu30_PtcR_Eta_ass = ()
    histo_reco_jet_Npu30_PtcR_Eta_ass = ()
    
    
    histo_assoc_Npu50_EtacR_Pt_ass = ()
    histo_sim_jet_Npu50_EtacR_Pt_ass = ()    
    
    histo_assoc_Npu50_PtcR_Eta_ass = ()
    histo_sim_jet_Npu50_PtcR_Eta_ass = ()    
    
    histo_assoc2_Npu50_EtacR_Pt_ass = ()
    histo_reco_jet_Npu50_EtacR_Pt_ass = ()    
    
    histo_assoc2_Npu50_PtcR_Eta_ass = ()
    histo_reco_jet_Npu50_PtcR_Eta_ass = ()
    
	
    for dir_ite in range(1):
		
	help = random.random()
    
        histo_assoc_NpucR_EtacR_Pt_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_EtacR" + str(help)),
        histo_sim_jet_NpucR_EtacR_Pt_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_EtacR" + str(help)),
    
        histo_assoc_NpucR_EtacR_Pt_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Eta13" + str(help),19,33),
        histo_sim_jet_NpucR_EtacR_Pt_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Eta13" + str(help),19,33),
    
        histo_assoc_NpucR_EtacR_Pt_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Eta25" + str(help),13,18),
        histo_sim_jet_NpucR_EtacR_Pt_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Eta25" + str(help),13,18),
    
        histo_assoc_NpucR_EtacR_Pt_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Eta50" + str(help),1,12),
        histo_sim_jet_NpucR_EtacR_Pt_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Eta50" + str(help),1,12),
    
    
        histo_assoc_NpucR_PtcR_Eta_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_PtcR" + str(help)),
        histo_sim_jet_NpucR_PtcR_Eta_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_PtcR" + str(help)),
    
        histo_assoc_NpucR_PtcR_Eta_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Pt100" + str(help),1,5),
        histo_sim_jet_NpucR_PtcR_Eta_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Pt100" + str(help),1,5),
    
        histo_assoc_NpucR_PtcR_Eta_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Pt300" + str(help),6,15),
        histo_sim_jet_NpucR_PtcR_Eta_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Pt300" + str(help),6,15),
    
        histo_assoc_NpucR_PtcR_Eta_ass+= histo_assoc_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_NpucR_Pt1000" + str(help),16,50),
        histo_sim_jet_NpucR_PtcR_Eta_ass+=  histo_sim_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_NpucR_Pt1000" + str(help),16,50),
	
    
        histo_assoc2_NpucR_EtacR_Pt_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_EtacR" + str(help)),
        histo_reco_jet_NpucR_EtacR_Pt_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_EtacR" + str(help)),
    
        histo_assoc2_NpucR_EtacR_Pt_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Eta13" + str(help),19,33),
        histo_reco_jet_NpucR_EtacR_Pt_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Eta13" + str(help),19,33),
    
        histo_assoc2_NpucR_EtacR_Pt_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Eta25" + str(help),13,18),
        histo_reco_jet_NpucR_EtacR_Pt_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Eta25" + str(help),13,18),
    
        histo_assoc2_NpucR_EtacR_Pt_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Eta50" + str(help),1,12),
        histo_reco_jet_NpucR_EtacR_Pt_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Eta50" + str(help),1,12),
	
    
        histo_assoc2_NpucR_PtcR_Eta_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_PtcR" + str(help)),
        histo_reco_jet_NpucR_PtcR_Eta_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_PtcR" + str(help)),
	
        histo_assoc2_NpucR_PtcR_Eta_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Pt100" + str(help),1,5),
        histo_reco_jet_NpucR_PtcR_Eta_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Pt100" + str(help),1,5),
	
        histo_assoc2_NpucR_PtcR_Eta_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Pt300" + str(help),6,15),
        histo_reco_jet_NpucR_PtcR_Eta_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Pt300" + str(help),6,15),
	
        histo_assoc2_NpucR_PtcR_Eta_ass+= histo_assoc2_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_NpucR_Pt1000" + str(help),16,50),
        histo_reco_jet_NpucR_PtcR_Eta_ass+=  histo_reco_jet_PtcR_NpucR_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_NpucR_Pt1000" + str(help),16,50),
	
	
	##Npu 15
    
        histo_assoc_Npu15_EtacR_Pt_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_EtacR"),
        histo_sim_jet_Npu15_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_EtacR"),
    
        histo_assoc_Npu15_EtacR_Pt_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Eta13",19,33),
        histo_sim_jet_Npu15_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Eta13",19,33),
    
        histo_assoc_Npu15_EtacR_Pt_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Eta25",13,18),
        histo_sim_jet_Npu15_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Eta25",13,18),
    
        histo_assoc_Npu15_EtacR_Pt_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Eta50",1,12),
        histo_sim_jet_Npu15_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Eta50",1,12),
    
    
        histo_assoc_Npu15_PtcR_Eta_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_PtcR"),
        histo_sim_jet_Npu15_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_PtcR"),
    
        histo_assoc_Npu15_PtcR_Eta_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Pt100",1,5),
        histo_sim_jet_Npu15_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Pt100",1,5),
    
        histo_assoc_Npu15_PtcR_Eta_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Pt300",6,15),
        histo_sim_jet_Npu15_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Pt300",6,15),
    
        histo_assoc_Npu15_PtcR_Eta_ass+= histo_assoc_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu15_Pt1000",16,50),
        histo_sim_jet_Npu15_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu15_Pt1000",16,50),
	
    
        histo_assoc2_Npu15_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_EtacR"),
        histo_reco_jet_Npu15_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu15_EtacR"),
    
        histo_assoc2_Npu15_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Eta13",19,33),
        histo_reco_jet_Npu15_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu15_Eta13",19,33),
    
        histo_assoc2_Npu15_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Eta25",13,18),
        histo_reco_jet_Npu15_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu15_Eta25",13,18),
    
        histo_assoc2_Npu15_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Eta50",1,12),
        histo_reco_jet_Npu15_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu15_Eta50",1,12),
    
    
        histo_assoc2_Npu15_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_PtcR"),
        histo_reco_jet_Npu15_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_reco_Npu15_PtcR"),
    
        histo_assoc2_Npu15_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Pt100",1,5),
        histo_reco_jet_Npu15_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_reco_Npu15_Pt100",1,5),
    
        histo_assoc2_Npu15_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Pt300",6,15),
        histo_reco_jet_Npu15_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_reco_Npu15_Pt300",6,15),
    
        histo_assoc2_Npu15_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu15_Pt1000",16,50),
        histo_reco_jet_Npu15_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu15_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_reco_Npu15_Pt1000",16,50),
	
	##Npu 30
    
        histo_assoc_Npu30_EtacR_Pt_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_EtacR"),
        histo_sim_jet_Npu30_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_EtacR"),
    
        histo_assoc_Npu30_EtacR_Pt_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Eta13",19,33),
        histo_sim_jet_Npu30_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Eta13",19,33),
    
        histo_assoc_Npu30_EtacR_Pt_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Eta25",13,18),
        histo_sim_jet_Npu30_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Eta25",13,18),
    
        histo_assoc_Npu30_EtacR_Pt_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Eta50",1,12),
        histo_sim_jet_Npu30_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Eta50",1,12),
    
    
        histo_assoc_Npu30_PtcR_Eta_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_PtcR"),
        histo_sim_jet_Npu30_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_PtcR"),
    
        histo_assoc_Npu30_PtcR_Eta_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Pt100",1,5),
        histo_sim_jet_Npu30_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Pt100",1,5),
    
        histo_assoc_Npu30_PtcR_Eta_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Pt300",6,15),
        histo_sim_jet_Npu30_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Pt300",6,15),
    
        histo_assoc_Npu30_PtcR_Eta_ass+= histo_assoc_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu30_Pt1000",16,50),
        histo_sim_jet_Npu30_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu30_Pt1000",16,50),
	
    
        histo_assoc2_Npu30_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_EtacR"),
        histo_reco_jet_Npu30_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_EtacR"),
    
        histo_assoc2_Npu30_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Eta13",19,33),
        histo_reco_jet_Npu30_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Eta13",19,33),
    
        histo_assoc2_Npu30_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Eta25",13,18),
        histo_reco_jet_Npu30_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Eta25",13,18),
    
        histo_assoc2_Npu30_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Eta50",1,12),
        histo_reco_jet_Npu30_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Eta50",1,12),
	
    
        histo_assoc2_Npu30_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_PtcR"),
        histo_reco_jet_Npu30_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_PtcR"),
	
        histo_assoc2_Npu30_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Pt100",1,5),
        histo_reco_jet_Npu30_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Pt100",1,5),
	
        histo_assoc2_Npu30_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Pt300",6,15),
        histo_reco_jet_Npu30_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Pt300",6,15),
	
        histo_assoc2_Npu30_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu30_Pt1000",16,50),
        histo_reco_jet_Npu30_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu30_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu30_Pt1000",16,50),
	
	
	##Npu 50
    
        histo_assoc_Npu50_EtacR_Pt_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_EtacR"),
        histo_sim_jet_Npu50_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_EtacR"),
    
        histo_assoc_Npu50_EtacR_Pt_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Eta13",19,33),
        histo_sim_jet_Npu50_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Eta13",19,33),
    
        histo_assoc_Npu50_EtacR_Pt_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Eta25",13,18),
        histo_sim_jet_Npu50_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Eta25",13,18),
    
        histo_assoc_Npu50_EtacR_Pt_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Eta50",1,12),
        histo_sim_jet_Npu50_EtacR_Pt_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Eta50",1,12),
    
    
        histo_assoc_Npu50_PtcR_Eta_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_PtcR"),
        histo_sim_jet_Npu50_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_PtcR"),
    
        histo_assoc_Npu50_PtcR_Eta_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Pt100",1,5),
        histo_sim_jet_Npu50_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Pt100",1,5),
    
        histo_assoc_Npu50_PtcR_Eta_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Pt300",6,15),
        histo_sim_jet_Npu50_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Pt300",6,15),
    
        histo_assoc_Npu50_PtcR_Eta_ass+= histo_assoc_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "effic_Npu50_Pt1000",16,50),
        histo_sim_jet_Npu50_PtcR_Eta_ass+=  histo_sim_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "sim_jet_Npu50_Pt1000",16,50),
	
    
        histo_assoc2_Npu50_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_EtacR"),
        histo_reco_jet_Npu50_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_EtacR"),
    
        histo_assoc2_Npu50_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Eta13",19,33),
        histo_reco_jet_Npu50_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Eta13",19,33),
    
        histo_assoc2_Npu50_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Eta25",13,18),
        histo_reco_jet_Npu50_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Eta25",13,18),
    
        histo_assoc2_Npu50_EtacR_Pt_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Eta50",1,12),
        histo_reco_jet_Npu50_EtacR_Pt_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionY(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Eta50",1,12),
	
    
        histo_assoc2_Npu50_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_PtcR"),
        histo_reco_jet_Npu50_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_PtcR"),
	
        histo_assoc2_Npu50_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Pt100",1,5),
        histo_reco_jet_Npu50_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Pt100",1,5),
	
        histo_assoc2_Npu50_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Pt300",6,15),
        histo_reco_jet_Npu50_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Pt300",6,15),
	
        histo_assoc2_Npu50_PtcR_Eta_ass+= histo_assoc2_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "fakerate_Npu50_Pt1000",16,50),
        histo_reco_jet_Npu50_PtcR_Eta_ass+=  histo_reco_jet_PtcR_Npu50_EtacR_com[ass_ite][dir_ite].ProjectionX(subdirnames[dir_ite][ass_ite]+ "reco_jet_Npu50_Pt1000",16,50),
	
    
    histo_assoc_NpucR_EtacR_Pt+= histo_assoc_NpucR_EtacR_Pt_ass,
    histo_sim_jet_NpucR_EtacR_Pt+= histo_sim_jet_NpucR_EtacR_Pt_ass,   
    
    histo_assoc_NpucR_PtcR_Eta+= histo_assoc_NpucR_PtcR_Eta_ass,
    histo_sim_jet_NpucR_PtcR_Eta+= histo_sim_jet_NpucR_PtcR_Eta_ass,	
    
    histo_assoc2_NpucR_EtacR_Pt+= histo_assoc2_NpucR_EtacR_Pt_ass,
    histo_reco_jet_NpucR_EtacR_Pt+= histo_reco_jet_NpucR_EtacR_Pt_ass,    
    
    histo_assoc2_NpucR_PtcR_Eta+= histo_assoc2_NpucR_PtcR_Eta_ass,
    histo_reco_jet_NpucR_PtcR_Eta+= histo_reco_jet_NpucR_PtcR_Eta_ass,
    
    
    ##Npu 15
    
    histo_assoc_Npu15_EtacR_Pt+= histo_assoc_Npu15_EtacR_Pt_ass,
    histo_sim_jet_Npu15_EtacR_Pt+= histo_sim_jet_Npu15_EtacR_Pt_ass,    
    
    histo_assoc_Npu15_PtcR_Eta+= histo_assoc_Npu15_PtcR_Eta_ass,
    histo_sim_jet_Npu15_PtcR_Eta+= histo_sim_jet_Npu15_PtcR_Eta_ass,	
    
    histo_assoc2_Npu15_EtacR_Pt+= histo_assoc2_Npu15_EtacR_Pt_ass,
    histo_reco_jet_Npu15_EtacR_Pt+= histo_reco_jet_Npu15_EtacR_Pt_ass,    
    
    histo_assoc2_Npu15_PtcR_Eta+= histo_assoc2_Npu15_PtcR_Eta_ass,
    histo_reco_jet_Npu15_PtcR_Eta+= histo_reco_jet_Npu15_PtcR_Eta_ass,
    
    
    ##Npu 30
    
    histo_assoc_Npu30_EtacR_Pt+= histo_assoc_Npu30_EtacR_Pt_ass,
    histo_sim_jet_Npu30_EtacR_Pt+= histo_sim_jet_Npu30_EtacR_Pt_ass,    
    
    histo_assoc_Npu30_PtcR_Eta+= histo_assoc_Npu30_PtcR_Eta_ass,
    histo_sim_jet_Npu30_PtcR_Eta+= histo_sim_jet_Npu30_PtcR_Eta_ass,	
    
    histo_assoc2_Npu30_EtacR_Pt+= histo_assoc2_Npu30_EtacR_Pt_ass,
    histo_reco_jet_Npu30_EtacR_Pt+= histo_reco_jet_Npu30_EtacR_Pt_ass,    
    
    histo_assoc2_Npu30_PtcR_Eta+= histo_assoc2_Npu30_PtcR_Eta_ass,
    histo_reco_jet_Npu30_PtcR_Eta+= histo_reco_jet_Npu30_PtcR_Eta_ass,
    
    
    ##Npu 50
    
    histo_assoc_Npu50_EtacR_Pt+= histo_assoc_Npu50_EtacR_Pt_ass,
    histo_sim_jet_Npu50_EtacR_Pt+= histo_sim_jet_Npu50_EtacR_Pt_ass,    
    
    histo_assoc_Npu50_PtcR_Eta+= histo_assoc_Npu50_PtcR_Eta_ass,
    histo_sim_jet_Npu50_PtcR_Eta+= histo_sim_jet_Npu50_PtcR_Eta_ass,	
    
    histo_assoc2_Npu50_EtacR_Pt+= histo_assoc2_Npu50_EtacR_Pt_ass,
    histo_reco_jet_Npu50_EtacR_Pt+= histo_reco_jet_Npu50_EtacR_Pt_ass,    
    
    histo_assoc2_Npu50_PtcR_Eta+= histo_assoc2_Npu50_PtcR_Eta_ass,
    histo_reco_jet_Npu50_PtcR_Eta+= histo_reco_jet_Npu50_PtcR_Eta_ass,
    
    
print ""    
print "Creating 1D ratio" 
    
effic_1D_EtacR_NpucR_Pt = GetFractionTH1Fs(histo_assoc_NpucR_EtacR_Pt,histo_sim_jet_NpucR_EtacR_Pt,"effic")
effic_1D_PtcR_NpucR_Eta = GetFractionTH1Fs(histo_assoc_NpucR_PtcR_Eta,histo_sim_jet_NpucR_PtcR_Eta,"effic")
    
fakerate_1D_EtacR_NpucR_Pt = GetFractionTH1Fs(histo_assoc2_NpucR_EtacR_Pt,histo_reco_jet_NpucR_EtacR_Pt,"fakerate")
fakerate_1D_PtcR_NpucR_Eta = GetFractionTH1Fs(histo_assoc2_NpucR_PtcR_Eta,histo_reco_jet_NpucR_PtcR_Eta,"fakerate")
    
##Npu 00    

effic_1D_EtacR_Npu15_Pt = GetFractionTH1Fs(histo_assoc_Npu15_EtacR_Pt,histo_sim_jet_Npu15_EtacR_Pt,"effic")
effic_1D_PtcR_Npu15_Eta = GetFractionTH1Fs(histo_assoc_Npu15_PtcR_Eta,histo_sim_jet_Npu15_PtcR_Eta,"effic")
    
fakerate_1D_EtacR_Npu15_Pt = GetFractionTH1Fs(histo_assoc2_Npu15_EtacR_Pt,histo_reco_jet_Npu15_EtacR_Pt,"fakerate")
fakerate_1D_PtcR_Npu15_Eta = GetFractionTH1Fs(histo_assoc2_Npu15_PtcR_Eta,histo_reco_jet_Npu15_PtcR_Eta,"fakerate")
    
##Npu 15   

effic_1D_EtacR_Npu30_Pt = GetFractionTH1Fs(histo_assoc_Npu30_EtacR_Pt,histo_sim_jet_Npu30_EtacR_Pt,"effic")
effic_1D_PtcR_Npu30_Eta = GetFractionTH1Fs(histo_assoc_Npu30_PtcR_Eta,histo_sim_jet_Npu30_PtcR_Eta,"effic")
    
fakerate_1D_EtacR_Npu30_Pt = GetFractionTH1Fs(histo_assoc2_Npu30_EtacR_Pt,histo_reco_jet_Npu30_EtacR_Pt,"fakerate")
fakerate_1D_PtcR_Npu30_Eta = GetFractionTH1Fs(histo_assoc2_Npu30_PtcR_Eta,histo_reco_jet_Npu30_PtcR_Eta,"fakerate")
    
##Npu 50   

effic_1D_EtacR_Npu50_Pt = GetFractionTH1Fs(histo_assoc_Npu50_EtacR_Pt,histo_sim_jet_Npu50_EtacR_Pt,"effic")
effic_1D_PtcR_Npu50_Eta = GetFractionTH1Fs(histo_assoc_Npu50_PtcR_Eta,histo_sim_jet_Npu50_PtcR_Eta,"effic")
    
fakerate_1D_EtacR_Npu50_Pt = GetFractionTH1Fs(histo_assoc2_Npu50_EtacR_Pt,histo_reco_jet_Npu50_EtacR_Pt,"fakerate")
fakerate_1D_PtcR_Npu50_Eta = GetFractionTH1Fs(histo_assoc2_Npu50_PtcR_Eta,histo_reco_jet_Npu50_PtcR_Eta,"fakerate")


all_effic_1D_histos_Vs_eta = ()
all_effic_1D_histos_Vs_eta+= effic_1D_PtcR_NpucR_Eta,
all_effic_1D_histos_Vs_eta+= effic_1D_PtcR_Npu15_Eta,
all_effic_1D_histos_Vs_eta+= effic_1D_PtcR_Npu30_Eta,
all_effic_1D_histos_Vs_eta+= effic_1D_PtcR_Npu50_Eta,


all_fakerate_1D_histos_Vs_eta = ()
all_fakerate_1D_histos_Vs_eta+= fakerate_1D_PtcR_NpucR_Eta,
all_fakerate_1D_histos_Vs_eta+= fakerate_1D_PtcR_Npu15_Eta,
all_fakerate_1D_histos_Vs_eta+= fakerate_1D_PtcR_Npu30_Eta,
all_fakerate_1D_histos_Vs_eta+= fakerate_1D_PtcR_Npu50_Eta,


all_effic_1D_histos_Vs_pt = ()
all_effic_1D_histos_Vs_pt+= effic_1D_EtacR_NpucR_Pt,
all_effic_1D_histos_Vs_pt+= effic_1D_EtacR_Npu15_Pt,
all_effic_1D_histos_Vs_pt+= effic_1D_EtacR_Npu30_Pt,
all_effic_1D_histos_Vs_pt+= effic_1D_EtacR_Npu50_Pt,

all_fakerate_1D_histos_Vs_pt = ()
all_fakerate_1D_histos_Vs_pt+= fakerate_1D_EtacR_NpucR_Pt,
all_fakerate_1D_histos_Vs_pt+= fakerate_1D_EtacR_Npu15_Pt,
all_fakerate_1D_histos_Vs_pt+= fakerate_1D_EtacR_Npu30_Pt,
all_fakerate_1D_histos_Vs_pt+= fakerate_1D_EtacR_Npu50_Pt,

for ite1 in range(len(all_effic_1D_histos_Vs_eta)):
    for ite2 in range(len(all_effic_1D_histos_Vs_eta[0])):
        for ite3 in range(len(all_effic_1D_histos_Vs_eta[0][0])):
		
	    all_effic_1D_histos_Vs_eta[ite1][ite2][ite3].SetTitle("efficiency vs eta")
	    all_effic_1D_histos_Vs_eta[ite1][ite2][ite3].GetXaxis().SetTitle("#eta^{GEN}")
	    all_effic_1D_histos_Vs_eta[ite1][ite2][ite3].GetYaxis().SetTitle("efficiency")
		
	    all_fakerate_1D_histos_Vs_eta[ite1][ite2][ite3].SetTitle("fake rate vs eta")
	    all_fakerate_1D_histos_Vs_eta[ite1][ite2][ite3].GetXaxis().SetTitle("#eta^{GEN}")
	    all_fakerate_1D_histos_Vs_eta[ite1][ite2][ite3].GetYaxis().SetTitle("fake rate")
		
	    all_effic_1D_histos_Vs_pt[ite1][ite2][ite3].SetTitle("efficiency vs pt")
	    all_effic_1D_histos_Vs_pt[ite1][ite2][ite3].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
	    all_effic_1D_histos_Vs_pt[ite1][ite2][ite3].GetYaxis().SetTitle("efficiency")
		
	    all_fakerate_1D_histos_Vs_pt[ite1][ite2][ite3].SetTitle("fake rate vs pt")
	    all_fakerate_1D_histos_Vs_pt[ite1][ite2][ite3].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
	    all_fakerate_1D_histos_Vs_pt[ite1][ite2][ite3].GetYaxis().SetTitle("fake rate")
	    all_fakerate_1D_histos_Vs_pt[ite1][ite2][ite3].GetXaxis().SetRangeUser(0.,300.)


all_effic_PtcR_EtacR_NpucR = ()
all_effic_PtcR_EtacR_Npu15 = ()
all_effic_PtcR_EtacR_Npu30 = ()
all_effic_PtcR_EtacR_Npu50 = ()

all_fakerate_PtcR_EtacR_NpucR = ()
all_fakerate_PtcR_EtacR_Npu15 = ()
all_fakerate_PtcR_EtacR_Npu30 = ()
all_fakerate_PtcR_EtacR_Npu50 = ()

for ass_ite in range(ass_number):
    all_effic_PtcR_EtacR_NpucR+= effic_PtcR_NpucR_EtacR[ass_ite][0],
    all_effic_PtcR_EtacR_Npu15+= effic_PtcR_Npu15_EtacR[ass_ite][0],
    all_effic_PtcR_EtacR_Npu30+= effic_PtcR_Npu30_EtacR[ass_ite][0],
    all_effic_PtcR_EtacR_Npu50+= effic_PtcR_Npu50_EtacR[ass_ite][0],
    
    all_fakerate_PtcR_EtacR_NpucR+= fakerate_PtcR_NpucR_EtacR[ass_ite][0],
    all_fakerate_PtcR_EtacR_Npu15+= fakerate_PtcR_Npu15_EtacR[ass_ite][0],
    all_fakerate_PtcR_EtacR_Npu30+= fakerate_PtcR_Npu30_EtacR[ass_ite][0],
    all_fakerate_PtcR_EtacR_Npu50+= fakerate_PtcR_Npu50_EtacR[ass_ite][0],
    
 
all_effic_histo = ()
all_effic_histo+= all_effic_PtcR_EtacR_NpucR,
all_effic_histo+= all_effic_PtcR_EtacR_Npu15,
all_effic_histo+= all_effic_PtcR_EtacR_Npu30,
all_effic_histo+= all_effic_PtcR_EtacR_Npu50,
    
all_fakerate_histo = ()
all_fakerate_histo+= all_fakerate_PtcR_EtacR_NpucR,
all_fakerate_histo+= all_fakerate_PtcR_EtacR_Npu15,
all_fakerate_histo+= all_fakerate_PtcR_EtacR_Npu30,
all_fakerate_histo+= all_fakerate_PtcR_EtacR_Npu50,

print ""   

DrawTH2Fs(all_effic_histo,"JetAnlzr_"+name+"_Effic","COLZ")
DrawTH2Fs(all_fakerate_histo,"JetAnlzr_"+name+"_Fakerate","COLZ")

DrawTH1Fs(all_effic_1D_histos_Vs_eta,"JetAnlzr_"+name+"_Effic_1D_Eta",True,True,True)
DrawTH1Fs(all_fakerate_1D_histos_Vs_eta,"JetAnlzr_"+name+"_Fakerate_1D_Eta",True,True,True)

DrawTH1Fs(all_effic_1D_histos_Vs_pt,"JetAnlzr_"+name+"_Effic_1D_Pt",True,True,True)
DrawTH1Fs(all_fakerate_1D_histos_Vs_pt,"JetAnlzr_"+name+"_Fakerate_1D_Pt",True,True,True)


###for presentations
  
canv = ROOT.TCanvas("canvPres","Guck mal",10,10,1200,850)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetTitleX(0.58)
ROOT.gStyle.SetTitleY(0.93)
ROOT.gStyle.SetTitleW(1.)
ROOT.gStyle.SetOptStat(0)

canv.Divide(2,2)

ROOT.gPad.SetGrid()
ROOT.gPad.SetBorderMode(0)
ROOT.gPad.SetFillColor(0)
		    
canv.cd(1)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

effic_1D_PtcR_NpucR_Eta[0][0].SetTitle("efficiency vs eta")
effic_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetTitle("efficiency")
effic_1D_PtcR_NpucR_Eta[0][0].GetXaxis().SetTitle("#eta^{GEN}")

effic_1D_PtcR_NpucR_Eta[0][0].GetXaxis().CenterTitle()
effic_1D_PtcR_NpucR_Eta[0][0].GetYaxis().CenterTitle()
  
effic_1D_PtcR_NpucR_Eta[0][0].GetXaxis().SetTitleOffset(1.1)
effic_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetTitleOffset(1.1)
effic_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetRangeUser(0.4,1.05)
	
color = my_colorsP[0]
effic_1D_PtcR_NpucR_Eta[0][0].SetLineColor(color)
effic_1D_PtcR_NpucR_Eta[0][0].SetLineWidth(2)
effic_1D_PtcR_NpucR_Eta[0][0].SetMarkerStyle(20)
effic_1D_PtcR_NpucR_Eta[0][0].SetMarkerSize(0.8)
effic_1D_PtcR_NpucR_Eta[0][0].SetMarkerColor(color)

effic_1D_PtcR_NpucR_Eta[0][0].Draw()

etaNo1_shadow.Draw("same")
etaNo2_shadow.Draw("same")
etaEc1_shadow.Draw("same")
etaEc2_shadow.Draw("same")

effic_1D_PtcR_NpucR_Eta[0][0].Draw("same")
	
color = my_colorsP[1]
effic_1D_PtcR_NpucR_Eta[1][0].SetLineColor(color)
effic_1D_PtcR_NpucR_Eta[1][0].SetLineWidth(2)
effic_1D_PtcR_NpucR_Eta[1][0].SetMarkerStyle(20)
effic_1D_PtcR_NpucR_Eta[1][0].SetMarkerSize(0.8)
effic_1D_PtcR_NpucR_Eta[1][0].SetMarkerColor(color)

effic_1D_PtcR_NpucR_Eta[1][0].Draw("same")
	
color = my_colorsP[2]
effic_1D_PtcR_NpucR_Eta[2][0].SetLineColor(color)
effic_1D_PtcR_NpucR_Eta[2][0].SetLineWidth(2)
effic_1D_PtcR_NpucR_Eta[2][0].SetMarkerStyle(20)
effic_1D_PtcR_NpucR_Eta[2][0].SetMarkerSize(0.8)
effic_1D_PtcR_NpucR_Eta[2][0].SetMarkerColor(color)

effic_1D_PtcR_NpucR_Eta[2][0].Draw("same")
	
color = my_colorsP[3]
effic_1D_PtcR_NpucR_Eta[3][0].SetLineColor(color)
effic_1D_PtcR_NpucR_Eta[3][0].SetLineWidth(2)
effic_1D_PtcR_NpucR_Eta[3][0].SetMarkerStyle(20)
effic_1D_PtcR_NpucR_Eta[3][0].SetMarkerSize(0.8)
effic_1D_PtcR_NpucR_Eta[3][0].SetMarkerColor(color)

effic_1D_PtcR_NpucR_Eta[3][0].Draw("same")
		    
canv.cd(2)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)
ROOT.gPad.SetLogx()

effic_1D_EtacR_NpucR_Pt[0][0].SetTitle("efficiency vs p_{t}")
effic_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetTitle("efficiency")
effic_1D_EtacR_NpucR_Pt[0][0].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")

effic_1D_EtacR_NpucR_Pt[0][0].GetXaxis().CenterTitle()
effic_1D_EtacR_NpucR_Pt[0][0].GetYaxis().CenterTitle()
  
effic_1D_EtacR_NpucR_Pt[0][0].GetXaxis().SetTitleOffset(1.1)
effic_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetTitleOffset(1.1)
effic_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetRangeUser(0.6,1.05)
	
color = my_colorsP[0]
effic_1D_EtacR_NpucR_Pt[0][0].SetLineColor(color)
effic_1D_EtacR_NpucR_Pt[0][0].SetLineWidth(2)
effic_1D_EtacR_NpucR_Pt[0][0].SetMarkerStyle(20)
effic_1D_EtacR_NpucR_Pt[0][0].SetMarkerSize(0.8)
effic_1D_EtacR_NpucR_Pt[0][0].SetMarkerColor(color)

effic_1D_EtacR_NpucR_Pt[0][0].Draw()
	
color = my_colorsP[1]
effic_1D_EtacR_NpucR_Pt[1][0].SetLineColor(color)
effic_1D_EtacR_NpucR_Pt[1][0].SetLineWidth(2)
effic_1D_EtacR_NpucR_Pt[1][0].SetMarkerStyle(20)
effic_1D_EtacR_NpucR_Pt[1][0].SetMarkerSize(0.8)
effic_1D_EtacR_NpucR_Pt[1][0].SetMarkerColor(color)

effic_1D_EtacR_NpucR_Pt[1][0].Draw("same")
	
color = my_colorsP[2]
effic_1D_EtacR_NpucR_Pt[2][0].SetLineColor(color)
effic_1D_EtacR_NpucR_Pt[2][0].SetLineWidth(2)
effic_1D_EtacR_NpucR_Pt[2][0].SetMarkerStyle(20)
effic_1D_EtacR_NpucR_Pt[2][0].SetMarkerSize(0.8)
effic_1D_EtacR_NpucR_Pt[2][0].SetMarkerColor(color)

effic_1D_EtacR_NpucR_Pt[2][0].Draw("same")
	
color = my_colorsP[3]
effic_1D_EtacR_NpucR_Pt[3][0].SetLineColor(color)
effic_1D_EtacR_NpucR_Pt[3][0].SetLineWidth(2)
effic_1D_EtacR_NpucR_Pt[3][0].SetMarkerStyle(20)
effic_1D_EtacR_NpucR_Pt[3][0].SetMarkerSize(0.8)
effic_1D_EtacR_NpucR_Pt[3][0].SetMarkerColor(color)

effic_1D_EtacR_NpucR_Pt[3][0].Draw("same")
		    
canv.cd(3)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

fakerate_1D_PtcR_NpucR_Eta[0][0].SetTitle("fake rate vs eta")
fakerate_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetTitle("fake rate")
fakerate_1D_PtcR_NpucR_Eta[0][0].GetXaxis().SetTitle("#eta^{GEN}")

fakerate_1D_PtcR_NpucR_Eta[0][0].GetXaxis().CenterTitle()
fakerate_1D_PtcR_NpucR_Eta[0][0].GetYaxis().CenterTitle()
  
fakerate_1D_PtcR_NpucR_Eta[0][0].GetXaxis().SetTitleOffset(1.1)
fakerate_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetTitleOffset(1.1)
fakerate_1D_PtcR_NpucR_Eta[0][0].GetYaxis().SetRangeUser(0.4,1.05)
	
color = my_colorsP[0]
fakerate_1D_PtcR_NpucR_Eta[0][0].SetLineColor(color)
fakerate_1D_PtcR_NpucR_Eta[0][0].SetLineWidth(2)
fakerate_1D_PtcR_NpucR_Eta[0][0].SetMarkerStyle(20)
fakerate_1D_PtcR_NpucR_Eta[0][0].SetMarkerSize(0.8)
fakerate_1D_PtcR_NpucR_Eta[0][0].SetMarkerColor(color)

fakerate_1D_PtcR_NpucR_Eta[0][0].Draw()

etaNo1_shadow.Draw("same")
etaNo2_shadow.Draw("same")
etaEc1_shadow.Draw("same")
etaEc2_shadow.Draw("same")

fakerate_1D_PtcR_NpucR_Eta[0][0].Draw("same")
	
color = my_colorsP[1]
fakerate_1D_PtcR_NpucR_Eta[1][0].SetLineColor(color)
fakerate_1D_PtcR_NpucR_Eta[1][0].SetLineWidth(2)
fakerate_1D_PtcR_NpucR_Eta[1][0].SetMarkerStyle(20)
fakerate_1D_PtcR_NpucR_Eta[1][0].SetMarkerSize(0.8)
fakerate_1D_PtcR_NpucR_Eta[1][0].SetMarkerColor(color)

fakerate_1D_PtcR_NpucR_Eta[1][0].Draw("same")
	
color = my_colorsP[2]
fakerate_1D_PtcR_NpucR_Eta[2][0].SetLineColor(color)
fakerate_1D_PtcR_NpucR_Eta[2][0].SetLineWidth(2)
fakerate_1D_PtcR_NpucR_Eta[2][0].SetMarkerStyle(20)
fakerate_1D_PtcR_NpucR_Eta[2][0].SetMarkerSize(0.8)
fakerate_1D_PtcR_NpucR_Eta[2][0].SetMarkerColor(color)

fakerate_1D_PtcR_NpucR_Eta[2][0].Draw("same")
	
color = my_colorsP[3]
fakerate_1D_PtcR_NpucR_Eta[3][0].SetLineColor(color)
fakerate_1D_PtcR_NpucR_Eta[3][0].SetLineWidth(2)
fakerate_1D_PtcR_NpucR_Eta[3][0].SetMarkerStyle(20)
fakerate_1D_PtcR_NpucR_Eta[3][0].SetMarkerSize(0.8)
fakerate_1D_PtcR_NpucR_Eta[3][0].SetMarkerColor(color)

fakerate_1D_PtcR_NpucR_Eta[3][0].Draw("same")
		    
canv.cd(4)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)
ROOT.gPad.SetLogx()

fakerate_1D_EtacR_NpucR_Pt[0][0].SetTitle("fake rate vs p_{t}")
fakerate_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetTitle("fake rate")
fakerate_1D_EtacR_NpucR_Pt[0][0].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")

fakerate_1D_EtacR_NpucR_Pt[0][0].GetXaxis().CenterTitle()
fakerate_1D_EtacR_NpucR_Pt[0][0].GetYaxis().CenterTitle()
  
fakerate_1D_EtacR_NpucR_Pt[0][0].GetXaxis().SetTitleOffset(1.1)
fakerate_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetTitleOffset(1.1)
fakerate_1D_EtacR_NpucR_Pt[0][0].GetYaxis().SetRangeUser(0.,1.05)
	
color = my_colorsP[0]
fakerate_1D_EtacR_NpucR_Pt[0][0].SetLineColor(color)
fakerate_1D_EtacR_NpucR_Pt[0][0].SetLineWidth(2)
fakerate_1D_EtacR_NpucR_Pt[0][0].SetMarkerStyle(20)
fakerate_1D_EtacR_NpucR_Pt[0][0].SetMarkerSize(0.8)
fakerate_1D_EtacR_NpucR_Pt[0][0].SetMarkerColor(color)

fakerate_1D_EtacR_NpucR_Pt[0][0].Draw()
	
color = my_colorsP[1]
fakerate_1D_EtacR_NpucR_Pt[1][0].SetLineColor(color)
fakerate_1D_EtacR_NpucR_Pt[1][0].SetLineWidth(2)
fakerate_1D_EtacR_NpucR_Pt[1][0].SetMarkerStyle(20)
fakerate_1D_EtacR_NpucR_Pt[1][0].SetMarkerSize(0.8)
fakerate_1D_EtacR_NpucR_Pt[1][0].SetMarkerColor(color)

fakerate_1D_EtacR_NpucR_Pt[1][0].Draw("same")
	
color = my_colorsP[2]
fakerate_1D_EtacR_NpucR_Pt[2][0].SetLineColor(color)
fakerate_1D_EtacR_NpucR_Pt[2][0].SetLineWidth(2)
fakerate_1D_EtacR_NpucR_Pt[2][0].SetMarkerStyle(20)
fakerate_1D_EtacR_NpucR_Pt[2][0].SetMarkerSize(0.8)
fakerate_1D_EtacR_NpucR_Pt[2][0].SetMarkerColor(color)

fakerate_1D_EtacR_NpucR_Pt[2][0].Draw("same")
	
color = my_colorsP[3]
fakerate_1D_EtacR_NpucR_Pt[3][0].SetLineColor(color)
fakerate_1D_EtacR_NpucR_Pt[3][0].SetLineWidth(2)
fakerate_1D_EtacR_NpucR_Pt[3][0].SetMarkerStyle(20)
fakerate_1D_EtacR_NpucR_Pt[3][0].SetMarkerSize(0.8)
fakerate_1D_EtacR_NpucR_Pt[3][0].SetMarkerColor(color)

fakerate_1D_EtacR_NpucR_Pt[3][0].Draw("same")

canv.SaveAs("../pictures/JetAnlzr_FakeEff_"+name+"_Presentation.png")
canv.SaveAs("../root/JetAnlzr_FakeEff_"+name+"_Presentation.root")
canv.SaveAs("../pdfs/JetAnlzr_FakeEff_"+name+"_Presentation.pdf")
