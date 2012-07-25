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
#dir_number = 1

subdirnames = ()

for i in range(dir_number):
	
    subdirnames+=GetListOfSubdirectories(filelist[0],dirnames[i],"ak5PFJets"),
    
ass_number = len(subdirnames[0])
#ass_number = 1

npu_ranges = GetListOfKeyWords(filelist[0],subdirnames[0][0],"ResolutionVsPtVsEta_")
    
npu_number = len(npu_ranges)

eta_number = len(eta_ranges)

pt_number = len(pt_ranges)
	
print ""

File_ref = [[]] * file_number

##the order is file - ass - dir - npu - pt - eta

histo_resolution_queue = ()

for file_ite in range(file_number):

    File_ref[file_ite] = ROOT.TFile.Open(filelist[file_ite])
    
    print "File " +str(file_ite+1)+ ": " + str(filelist[file_ite]) + " ...",
    
    histo_resolution_ass = ()
    
    for ass_ite in range(ass_number):
	    
        histo_resolution_ass_dir = ()
	
	for dir_ite in range(dir_number):

            histo_resolution_ass_dir_npu = ()   
	
	    for npu_ite in range(npu_number):

                histo_resolution_ass_dir_npu_pt = () 
	    
	        resolution_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+ npu_ranges[npu_ite])
	    
	        resolution_help.SetName(subdirnames[dir_ite][ass_ite] + npu_ranges[npu_ite] + str(help))
	
	        for pt_ite in range(pt_number):

                    histo_resolution_ass_dir_npu_pt_eta = () 
		    
		    if (pt_ite==(pt_number-1)):
		        ptMin = 0
		        ptMax = -1
		    else:
		        ptMin = GetBinNumber(resolution_help.GetXaxis(),pt_ranges[pt_ite])
		        ptMax = GetBinNumber(resolution_help.GetXaxis(),pt_ranges[pt_ite+1])
 	
	            for eta_ite in range(eta_number):
		    
		        if (eta_ite==(eta_number-1)):
 		            etaMin = 0
		            etaMax = -1
		        else:		    
		            etaMin = GetBinNumber(resolution_help.GetYaxis(),eta_ranges[eta_ite])
		            etaMax = GetBinNumber(resolution_help.GetYaxis(),eta_ranges[eta_ite+1])
 		
	                help = random.random()
	    	    
	                histo_resolution_ass_dir_npu_pt_eta+= resolution_help.ProjectionZ(subdirnames[dir_ite][ass_ite] + npu_ranges[npu_ite] + str(help),ptMin,ptMax,etaMin,etaMax),
	            
	    
                    histo_resolution_ass_dir_npu_pt+= histo_resolution_ass_dir_npu_pt_eta,
		
                histo_resolution_ass_dir_npu+= histo_resolution_ass_dir_npu_pt,
	    
            histo_resolution_ass_dir+= histo_resolution_ass_dir_npu,
	    
        histo_resolution_ass+= histo_resolution_ass_dir,
	
    histo_resolution_queue+= histo_resolution_ass,

    print "done" 
    
##the former order is file - ass - dir - npu - pt - eta
##the new order is ass - dir - npu - pt - eta

print ""
	    
histo_resolution_com = GetCombinedTH1FsWeight5D(histo_resolution_queue,filelist)

print "Rearranging"
    
##the former order is ass - dir - npu - pt - eta
##the new order is npu - pt - eta - ass - dir 
    
histo_resolution = ()
    
for npu_ite in range(npu_number):
    
    histo_resolution_npu = ()  
	
    for pt_ite in range(pt_number):
	    
        histo_resolution_npu_pt = ()
     
        for eta_ite in range(eta_number):
    
            histo_resolution_npu_pt_eta = ()  
	
            for ass_ite in range(ass_number):
    
                histo_resolution_npu_pt_eta_ass = ()  
	
                for dir_ite in range(dir_number):
			
		    histo_resolution_npu_pt_eta_ass+= histo_resolution_com[ass_ite][dir_ite][npu_ite][pt_ite][eta_ite],
		    
		histo_resolution_npu_pt_eta+= histo_resolution_npu_pt_eta_ass,
		
            histo_resolution_npu_pt+= histo_resolution_npu_pt_eta,
	    
	histo_resolution_npu+= histo_resolution_npu_pt,
	
    histo_resolution+= histo_resolution_npu,

print ""
    
##Set title
    
for npu_ite in range(npu_number): 
	
    for pt_ite in range(pt_number):
     
        for eta_ite in range(eta_number): 
	
            for dir_ite in range(dir_number):

                title = dirnames[dir_ite] + " - Eta" + eta_namesP[eta_ite]
		    
		histo_resolution[npu_ite][pt_ite][eta_ite][0][dir_ite].SetTitle("p_{t}-resolution - " + title)
		histo_resolution[npu_ite][pt_ite][eta_ite][0][dir_ite].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
		histo_resolution[npu_ite][pt_ite][eta_ite][0][dir_ite].GetYaxis().SetTitle("# entries")
   
#for npu_ite in range(npu_number): 
	
    #for pt_ite in range(pt_number):
	    	
        #DrawTH1Fs(histo_resolution[npu_ite][pt_ite],"JetAnlzr_"+name+"_"+npu_ranges[npu_ite]+"_Pt"+pt_namesP[pt_ite],True,True,False,True)


###for presentations
##the new order is npu - pt - eta - ass - dir 
  
canv = ROOT.TCanvas("canvPres","Guck mal",10,10,1200,850)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetTitleX(0.58)
ROOT.gStyle.SetTitleY(0.93)
ROOT.gStyle.SetTitleW(1.)
ROOT.gStyle.SetOptStat("mr")

canv.Divide(4,2)

ROOT.gPad.SetGrid()
ROOT.gPad.SetBorderMode(0)
ROOT.gPad.SetFillColor(0)
	
histo_resolution[0][3][3][0][0].SetTitle("Whole Npu - whole eta - whole pt")
histo_resolution[2][3][3][0][0].SetTitle("15<Npu<30 - whole eta - whole pt")
histo_resolution[0][0][3][0][0].SetTitle("Whole Npu - whole eta - 30<pt<100")
histo_resolution[0][3][0][0][0].SetTitle("Whole Npu - 0<|eta|<1.3 - whole pt")

histo_resolution[0][3][3][0][1].SetTitle("Whole Npu - whole eta - whole pt")
histo_resolution[2][3][3][0][1].SetTitle("15<Npu<30 - whole eta - whole pt")
histo_resolution[0][0][3][0][1].SetTitle("Whole Npu - whole eta - 30<pt<100")
histo_resolution[0][3][0][0][1].SetTitle("Whole Npu - 0<|eta|<1.3 - whole pt")

### pad 1 ####

canv.cd(1)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][3][3][1][0].GetMaximum()
histo_resolution[0][3][3][0][0].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][3][3][0][0].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][3][3][0][0].GetYaxis().SetTitle("# entries")

histo_resolution[0][3][3][0][0].GetXaxis().CenterTitle()
histo_resolution[0][3][3][0][0].GetYaxis().CenterTitle()
  
histo_resolution[0][3][3][0][0].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][3][3][0][0].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][3][3][0][0].SetLineColor(color)
histo_resolution[0][3][3][0][0].SetLineWidth(2)
histo_resolution[0][3][3][0][0].SetMarkerStyle(20)
histo_resolution[0][3][3][0][0].SetMarkerSize(0.8)
histo_resolution[0][3][3][0][0].SetMarkerColor(color)
histo_resolution[0][3][3][0][0].SetStats(1)

histo_resolution[0][3][3][0][0].Draw()
moveStatBox(canv, histo_resolution[0][3][3][0][0], color, 0, 1, 0.959136434808)
	
color = my_colorsP[1]
histo_resolution[0][3][3][1][0].SetLineColor(color)
histo_resolution[0][3][3][1][0].SetLineWidth(2)
histo_resolution[0][3][3][1][0].SetMarkerStyle(20)
histo_resolution[0][3][3][1][0].SetMarkerSize(0.8)
histo_resolution[0][3][3][1][0].SetMarkerColor(color)
histo_resolution[0][3][3][1][0].SetStats(1)

histo_resolution[0][3][3][1][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][3][1][0], color, 1, 1, 1.01832507387)
	
color = my_colorsP[3]
histo_resolution[0][3][3][3][0].SetLineColor(color)
histo_resolution[0][3][3][3][0].SetLineWidth(2)
histo_resolution[0][3][3][3][0].SetMarkerStyle(20)
histo_resolution[0][3][3][3][0].SetMarkerSize(0.8)
histo_resolution[0][3][3][3][0].SetMarkerColor(color)
histo_resolution[0][3][3][3][0].SetStats(1)

histo_resolution[0][3][3][3][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][3][3][0], color, 2, 1, 1.02370260367)

### pad 2 ####

canv.cd(2)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[2][3][3][1][0].GetMaximum()
histo_resolution[2][3][3][0][0].SetMaximum(maximum + 0.1*maximum)

histo_resolution[2][3][3][0][0].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[2][3][3][0][0].GetYaxis().SetTitle("# entries")

histo_resolution[2][3][3][0][0].GetXaxis().CenterTitle()
histo_resolution[2][3][3][0][0].GetYaxis().CenterTitle()
  
histo_resolution[2][3][3][0][0].GetXaxis().SetTitleOffset(1.1)
histo_resolution[2][3][3][0][0].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[2][3][3][0][0].SetLineColor(color)
histo_resolution[2][3][3][0][0].SetLineWidth(2)
histo_resolution[2][3][3][0][0].SetMarkerStyle(20)
histo_resolution[2][3][3][0][0].SetMarkerSize(0.8)
histo_resolution[2][3][3][0][0].SetMarkerColor(color)
histo_resolution[2][3][3][0][0].SetStats(1)

histo_resolution[2][3][3][0][0].Draw()
moveStatBox(canv, histo_resolution[2][3][3][0][0], color, 0, 1, 0.950417971718)
	
color = my_colorsP[1]
histo_resolution[2][3][3][1][0].SetLineColor(color)
histo_resolution[2][3][3][1][0].SetLineWidth(2)
histo_resolution[2][3][3][1][0].SetMarkerStyle(20)
histo_resolution[2][3][3][1][0].SetMarkerSize(0.8)
histo_resolution[2][3][3][1][0].SetMarkerColor(color)
histo_resolution[2][3][3][1][0].SetStats(1)

histo_resolution[2][3][3][1][0].Draw("sames")
moveStatBox(canv, histo_resolution[2][3][3][1][0], color, 1, 1, 1.00393765777)
	
color = my_colorsP[3]
histo_resolution[2][3][3][3][0].SetLineColor(color)
histo_resolution[2][3][3][3][0].SetLineWidth(2)
histo_resolution[2][3][3][3][0].SetMarkerStyle(20)
histo_resolution[2][3][3][3][0].SetMarkerSize(0.8)
histo_resolution[2][3][3][3][0].SetMarkerColor(color)
histo_resolution[2][3][3][3][0].SetStats(1)

histo_resolution[2][3][3][3][0].Draw("sames")
moveStatBox(canv, histo_resolution[2][3][3][3][0], color, 2, 1, 1.00854925201)

### pad 3 ####

canv.cd(3)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][0][3][1][0].GetMaximum()
histo_resolution[0][0][3][0][0].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][0][3][0][0].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][0][3][0][0].GetYaxis().SetTitle("# entries")

histo_resolution[0][0][3][0][0].GetXaxis().CenterTitle()
histo_resolution[0][0][3][0][0].GetYaxis().CenterTitle()
  
histo_resolution[0][0][3][0][0].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][0][3][0][0].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][0][3][0][0].SetLineColor(color)
histo_resolution[0][0][3][0][0].SetLineWidth(2)
histo_resolution[0][0][3][0][0].SetMarkerStyle(20)
histo_resolution[0][0][3][0][0].SetMarkerSize(0.8)
histo_resolution[0][0][3][0][0].SetMarkerColor(color)
histo_resolution[0][0][3][0][0].SetStats(1)

histo_resolution[0][0][3][0][0].Draw()
moveStatBox(canv, histo_resolution[0][0][3][0][0], color, 0, 1, 0.993336005557)
	
color = my_colorsP[1]
histo_resolution[0][0][3][1][0].SetLineColor(color)
histo_resolution[0][0][3][1][0].SetLineWidth(2)
histo_resolution[0][0][3][1][0].SetMarkerStyle(20)
histo_resolution[0][0][3][1][0].SetMarkerSize(0.8)
histo_resolution[0][0][3][1][0].SetMarkerColor(color)
histo_resolution[0][0][3][1][0].SetStats(1)

histo_resolution[0][0][3][1][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][0][3][1][0], color, 1, 1, 1.0654565857)
	
color = my_colorsP[3]
histo_resolution[0][0][3][3][0].SetLineColor(color)
histo_resolution[0][0][3][3][0].SetLineWidth(2)
histo_resolution[0][0][3][3][0].SetMarkerStyle(20)
histo_resolution[0][0][3][3][0].SetMarkerSize(0.8)
histo_resolution[0][0][3][3][0].SetMarkerColor(color)
histo_resolution[0][0][3][3][0].SetStats(1)

histo_resolution[0][0][3][3][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][0][3][3][0], color, 2, 1, 1.07368008622)

### pad 4 ####

canv.cd(4)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][3][0][1][0].GetMaximum()
histo_resolution[0][3][0][0][0].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][3][0][0][0].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][3][0][0][0].GetYaxis().SetTitle("# entries")

histo_resolution[0][3][0][0][0].GetXaxis().CenterTitle()
histo_resolution[0][3][0][0][0].GetYaxis().CenterTitle()
  
histo_resolution[0][3][0][0][0].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][3][0][0][0].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][3][0][0][0].SetLineColor(color)
histo_resolution[0][3][0][0][0].SetLineWidth(2)
histo_resolution[0][3][0][0][0].SetMarkerStyle(20)
histo_resolution[0][3][0][0][0].SetMarkerSize(0.8)
histo_resolution[0][3][0][0][0].SetMarkerColor(color)
histo_resolution[0][3][0][0][0].SetStats(1)

histo_resolution[0][3][0][0][0].Draw()
moveStatBox(canv, histo_resolution[0][3][0][0][0], color, 0, 1, 0.956241582749)
	
color = my_colorsP[1]
histo_resolution[0][3][0][1][0].SetLineColor(color)
histo_resolution[0][3][0][1][0].SetLineWidth(2)
histo_resolution[0][3][0][1][0].SetMarkerStyle(20)
histo_resolution[0][3][0][1][0].SetMarkerSize(0.8)
histo_resolution[0][3][0][1][0].SetMarkerColor(color)
histo_resolution[0][3][0][1][0].SetStats(1)

histo_resolution[0][3][0][1][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][0][1][0], color, 1, 1, 1.00850373412)
	
color = my_colorsP[3]
histo_resolution[0][3][0][3][0].SetLineColor(color)
histo_resolution[0][3][0][3][0].SetLineWidth(2)
histo_resolution[0][3][0][3][0].SetMarkerStyle(20)
histo_resolution[0][3][0][3][0].SetMarkerSize(0.8)
histo_resolution[0][3][0][3][0].SetMarkerColor(color)
histo_resolution[0][3][0][3][0].SetStats(1)

histo_resolution[0][3][0][3][0].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][0][3][0], color, 2, 1, 1.01379419791)
	
### pad 5 ####

canv.cd(5)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][3][3][1][1].GetMaximum()
histo_resolution[0][3][3][0][1].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][3][3][0][1].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][3][3][0][1].GetYaxis().SetTitle("# entries")

histo_resolution[0][3][3][0][1].GetXaxis().CenterTitle()
histo_resolution[0][3][3][0][1].GetYaxis().CenterTitle()
  
histo_resolution[0][3][3][0][1].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][3][3][0][1].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][3][3][0][1].SetLineColor(color)
histo_resolution[0][3][3][0][1].SetLineWidth(2)
histo_resolution[0][3][3][0][1].SetMarkerStyle(20)
histo_resolution[0][3][3][0][1].SetMarkerSize(0.8)
histo_resolution[0][3][3][0][1].SetMarkerColor(color)
histo_resolution[0][3][3][0][1].SetStats(1)

histo_resolution[0][3][3][0][1].Draw()
moveStatBox(canv, histo_resolution[0][3][3][0][1], color, 0, 1, 0.970753019661)
	
color = my_colorsP[1]
histo_resolution[0][3][3][1][1].SetLineColor(color)
histo_resolution[0][3][3][1][1].SetLineWidth(2)
histo_resolution[0][3][3][1][1].SetMarkerStyle(20)
histo_resolution[0][3][3][1][1].SetMarkerSize(0.8)
histo_resolution[0][3][3][1][1].SetMarkerColor(color)
histo_resolution[0][3][3][1][1].SetStats(1)

histo_resolution[0][3][3][1][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][3][1][1], color, 1, 1, 1.03429483127)
	
color = my_colorsP[3]
histo_resolution[0][3][3][3][1].SetLineColor(color)
histo_resolution[0][3][3][3][1].SetLineWidth(2)
histo_resolution[0][3][3][3][1].SetMarkerStyle(20)
histo_resolution[0][3][3][3][1].SetMarkerSize(0.8)
histo_resolution[0][3][3][3][1].SetMarkerColor(color)
histo_resolution[0][3][3][3][1].SetStats(1)

histo_resolution[0][3][3][3][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][3][3][1], color, 2, 1, 1.03992302592)

### pad 6 ####

canv.cd(6)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[2][3][3][1][1].GetMaximum()
histo_resolution[2][3][3][0][1].SetMaximum(maximum + 0.1*maximum)

histo_resolution[2][3][3][0][1].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[2][3][3][0][1].GetYaxis().SetTitle("# entries")

histo_resolution[2][3][3][0][1].GetXaxis().CenterTitle()
histo_resolution[2][3][3][0][1].GetYaxis().CenterTitle()
  
histo_resolution[2][3][3][0][1].GetXaxis().SetTitleOffset(1.1)
histo_resolution[2][3][3][0][1].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[2][3][3][0][1].SetLineColor(color)
histo_resolution[2][3][3][0][1].SetLineWidth(2)
histo_resolution[2][3][3][0][1].SetMarkerStyle(20)
histo_resolution[2][3][3][0][1].SetMarkerSize(0.8)
histo_resolution[2][3][3][0][1].SetMarkerColor(color)
histo_resolution[2][3][3][0][1].SetStats(1)

histo_resolution[2][3][3][0][1].Draw()
moveStatBox(canv, histo_resolution[2][3][3][0][1], color, 0, 1, 0.972284369393)
	
color = my_colorsP[1]
histo_resolution[2][3][3][1][1].SetLineColor(color)
histo_resolution[2][3][3][1][1].SetLineWidth(2)
histo_resolution[2][3][3][1][1].SetMarkerStyle(20)
histo_resolution[2][3][3][1][1].SetMarkerSize(0.8)
histo_resolution[2][3][3][1][1].SetMarkerColor(color)
histo_resolution[2][3][3][1][1].SetStats(1)

histo_resolution[2][3][3][1][1].Draw("sames")
moveStatBox(canv, histo_resolution[2][3][3][1][1], color, 1, 1, 1.02969370237)
	
color = my_colorsP[3]
histo_resolution[2][3][3][3][1].SetLineColor(color)
histo_resolution[2][3][3][3][1].SetLineWidth(2)
histo_resolution[2][3][3][3][1].SetMarkerStyle(20)
histo_resolution[2][3][3][3][1].SetMarkerSize(0.8)
histo_resolution[2][3][3][3][1].SetMarkerColor(color)
histo_resolution[2][3][3][3][1].SetStats(1)

histo_resolution[2][3][3][3][1].Draw("sames")
moveStatBox(canv, histo_resolution[2][3][3][3][1], color, 2, 1, 1.03450303332)

### pad 7 ####

canv.cd(7)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][0][3][1][1].GetMaximum()
histo_resolution[0][0][3][0][1].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][0][3][0][1].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][0][3][0][1].GetYaxis().SetTitle("# entries")

histo_resolution[0][0][3][0][1].GetXaxis().CenterTitle()
histo_resolution[0][0][3][0][1].GetYaxis().CenterTitle()
  
histo_resolution[0][0][3][0][1].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][0][3][0][1].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][0][3][0][1].SetLineColor(color)
histo_resolution[0][0][3][0][1].SetLineWidth(2)
histo_resolution[0][0][3][0][1].SetMarkerStyle(20)
histo_resolution[0][0][3][0][1].SetMarkerSize(0.8)
histo_resolution[0][0][3][0][1].SetMarkerColor(color)
histo_resolution[0][0][3][0][1].SetStats(1)

histo_resolution[0][0][3][0][1].Draw()
moveStatBox(canv, histo_resolution[0][0][3][0][1], color, 0, 1, 0.96673067545)
	
color = my_colorsP[1]
histo_resolution[0][0][3][1][1].SetLineColor(color)
histo_resolution[0][0][3][1][1].SetLineWidth(2)
histo_resolution[0][0][3][1][1].SetMarkerStyle(20)
histo_resolution[0][0][3][1][1].SetMarkerSize(0.8)
histo_resolution[0][0][3][1][1].SetMarkerColor(color)
histo_resolution[0][0][3][1][1].SetStats(1)

histo_resolution[0][0][3][1][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][0][3][1][1], color, 1, 1, 1.04426948645)
	
color = my_colorsP[3]
histo_resolution[0][0][3][3][1].SetLineColor(color)
histo_resolution[0][0][3][3][1].SetLineWidth(2)
histo_resolution[0][0][3][3][1].SetMarkerStyle(20)
histo_resolution[0][0][3][3][1].SetMarkerSize(0.8)
histo_resolution[0][0][3][3][1].SetMarkerColor(color)
histo_resolution[0][0][3][3][1].SetStats(1)

histo_resolution[0][0][3][3][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][0][3][3][1], color, 2, 1, 1.05262166404)

### pad 8 ####

canv.cd(8)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

maximum= histo_resolution[0][3][0][1][1].GetMaximum()
histo_resolution[0][3][0][0][1].SetMaximum(maximum + 0.1*maximum)

histo_resolution[0][3][0][0][1].GetXaxis().SetTitle("(p_{t}^{RECO} - p_{t}^{GEN}) / p_{t}^{GEN}")
histo_resolution[0][3][0][0][1].GetYaxis().SetTitle("# entries")

histo_resolution[0][3][0][0][1].GetXaxis().CenterTitle()
histo_resolution[0][3][0][0][1].GetYaxis().CenterTitle()
  
histo_resolution[0][3][0][0][1].GetXaxis().SetTitleOffset(1.1)
histo_resolution[0][3][0][0][1].GetYaxis().SetTitleOffset(1.1)
	
color = my_colorsP[0]
histo_resolution[0][3][0][0][1].SetLineColor(color)
histo_resolution[0][3][0][0][1].SetLineWidth(2)
histo_resolution[0][3][0][0][1].SetMarkerStyle(20)
histo_resolution[0][3][0][0][1].SetMarkerSize(0.8)
histo_resolution[0][3][0][0][1].SetMarkerColor(color)
histo_resolution[0][3][0][0][1].SetStats(1)

histo_resolution[0][3][0][0][1].Draw()
moveStatBox(canv, histo_resolution[0][3][0][0][1], color, 0, 1, 0.971857048601)
	
color = my_colorsP[1]
histo_resolution[0][3][0][1][1].SetLineColor(color)
histo_resolution[0][3][0][1][1].SetLineWidth(2)
histo_resolution[0][3][0][1][1].SetMarkerStyle(20)
histo_resolution[0][3][0][1][1].SetMarkerSize(0.8)
histo_resolution[0][3][0][1][1].SetMarkerColor(color)
histo_resolution[0][3][0][1][1].SetStats(1)

histo_resolution[0][3][0][1][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][0][1][1], color, 1, 1, 1.02746259988)
	
color = my_colorsP[3]
histo_resolution[0][3][0][3][1].SetLineColor(color)
histo_resolution[0][3][0][3][1].SetLineWidth(2)
histo_resolution[0][3][0][3][1].SetMarkerStyle(20)
histo_resolution[0][3][0][3][1].SetMarkerSize(0.8)
histo_resolution[0][3][0][3][1].SetMarkerColor(color)
histo_resolution[0][3][0][3][1].SetStats(1)

histo_resolution[0][3][0][3][1].Draw("sames")
moveStatBox(canv, histo_resolution[0][3][0][3][1], color, 2, 1, 1.03297060199)	

canv.SaveAs("../pictures/JetAnlzr_"+name+"_PtResolution_Presentation.png")
canv.SaveAs("../root/JetAnlzr_"+name+"_PtResolution_Presentation.root")
canv.SaveAs("../pdfs/JetAnlzr_"+name+"_PtResolution_Presentation.pdf")
