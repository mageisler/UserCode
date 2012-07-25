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

npu_ranges = GetListOfKeyWords(filelist[0],subdirnames[0][0],"ResponseVsPtVsEta_")
    
npu_number = len(npu_ranges)

eta_number = len(eta_ranges)
	
print ""

File_ref = [[]] * file_number

##the order is file - ass - dir - npu - eta 

histo_response_queue = ()

for file_ite in range(file_number):

    File_ref[file_ite] = ROOT.TFile.Open(filelist[file_ite])
    
    print "File " +str(file_ite+1)+ ": " + str(filelist[file_ite]) + " ...",

    histo_response_ass = ()
    
    for ass_ite in range(ass_number):

        histo_response_ass_dir = ()     
	
	for dir_ite in range(dir_number):

            histo_response_ass_dir_npu = ()   
	
	    for npu_ite in range(npu_number):

                histo_response_ass_dir_npu_eta = () 
		
	        help = random.random()
	    
	        response_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+ npu_ranges[npu_ite])
	    
	        response_help.SetName(subdirnames[dir_ite][ass_ite] + npu_ranges[npu_ite] + str(help))
		
		histo_response_ass_dir_npu_eta+= response_help.Project3D("cR_zxo"),
	
	        for eta_ite in range(eta_number-1):	
	    
	            response_help.GetYaxis().SetRangeUser(eta_ranges[eta_ite],eta_ranges[eta_ite+1])	    
	            histo_response_ass_dir_npu_eta+= response_help.Project3D(str(eta_ranges[eta_ite+1]) + "_zxo"),
	            
	    
                histo_response_ass_dir_npu+= histo_response_ass_dir_npu_eta,
	    
            histo_response_ass_dir+= histo_response_ass_dir_npu,
	    
        histo_response_ass+= histo_response_ass_dir,
	
    histo_response_queue+= histo_response_ass,

    print "done" 
    
##the former order is file - ass - dir - npu - eta
##the new order is ass - dir - npu - eta

print ""
	    
histo_response_com = GetCombinedTH2FsWeight4D(histo_response_queue,filelist)

##the former order is ass - dir - npu - eta
##the new order is npu - eta - ass - dir 
    
mean_response = ()
width_response = ()

print "Profiling"
    
for npu_ite in range(npu_number):
    
    mean_response_npu = ()       
    width_response_npu = ()
	
    for eta_ite in range(eta_number):
    
        mean_response_npu_eta = ()       
        width_response_npu_eta = ()
	
        for ass_ite in range(ass_number):
    
            mean_response_npu_eta_ass = ()       
            width_response_npu_eta_ass = ()
	
            for dir_ite in range(dir_number):
	    
	        mean_response_npu_eta_ass+= GetMean2(histo_response_com[ass_ite][dir_ite][npu_ite][eta_ite]), 
	        width_response_npu_eta_ass+= GetWidth2(histo_response_com[ass_ite][dir_ite][npu_ite][eta_ite]),
    
            mean_response_npu_eta+= mean_response_npu_eta_ass,        
            width_response_npu_eta+= width_response_npu_eta_ass,
    
        mean_response_npu+= mean_response_npu_eta,        
        width_response_npu+= width_response_npu_eta,
    
    mean_response+= mean_response_npu,        
    width_response+= width_response_npu,
    
print "" 
	
for npu_ite in range(npu_number):
            
    for eta_ite in range(eta_number):
	    
        for dir_ite in range(dir_number):

            title = dirnames[dir_ite] + " - Eta" + eta_names[eta_ite]
		    
	    mean_response[npu_ite][eta_ite][0][dir_ite].SetTitle("p_{t}-response - MEAN - " + title)
	    mean_response[npu_ite][eta_ite][0][dir_ite].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
	    mean_response[npu_ite][eta_ite][0][dir_ite].GetYaxis().SetTitle("< p_{t}^{RECO} / p_{t}^{GEN} >") 
	    mean_response[npu_ite][eta_ite][0][dir_ite].GetYaxis().SetRangeUser(0.7,1.3) 	
			    
	    width_response[npu_ite][eta_ite][0][dir_ite].SetTitle("p_{t}-response - WIDTH - " + title)
	    width_response[npu_ite][eta_ite][0][dir_ite].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
	    width_response[npu_ite][eta_ite][0][dir_ite].GetYaxis().SetTitle("#sigma_{p_{t}^{RECO}/p_{t}^{GEN}}") 
	

#for npu_ite in range(npu_number):
    #DrawTH1Fs(mean_response[npu_ite],"JetAnlzr_"+name+"_"+npu_ranges[npu_ite].split("_")[0]+"_Pt_"+npu_ranges[npu_ite].split("_")[1]+"_Mean",False,True,True)
    #DrawTH1Fs(width_response[npu_ite],"JetAnlzr_"+name+"_"+npu_ranges[npu_ite].split("_")[0]+"_Pt_"+npu_ranges[npu_ite].split("_")[1]+"_Width",True,True,True)


###for presentations
##the new order is npu - eta - ass - dir 
  
canv = ROOT.TCanvas("canvPres","Guck mal",10,10,1600,800)
ROOT.gStyle.SetNumberContours(255)
ROOT.gStyle.SetTitleX(0.58)
ROOT.gStyle.SetTitleY(0.93)
ROOT.gStyle.SetTitleW(1.)
ROOT.gStyle.SetOptStat(0)

canv.Divide(4,2)

ROOT.gPad.SetGrid()
ROOT.gPad.SetBorderMode(0)
ROOT.gPad.SetFillColor(0)
		    

width_response[0][0][0][0].SetTitle("WIDTH vs pt - Whole npu range")
width_response[1][0][0][0].SetTitle("WIDTH vs pt - npu < 15.")
width_response[2][0][0][0].SetTitle("WIDTH vs pt - 15. <= npu < 30.")
width_response[3][0][0][0].SetTitle("WIDTH vs pt - 30. <= npu")
		   
for pt_ite in range(npu_number):
    canv.cd(pt_ite+1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetTopMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetFillColor(0)
    
    width_response[pt_ite][0][0][0].GetYaxis().SetTitle("#sigma_{p_{t}^{RECO}/p_{t}^{GEN}}")
    width_response[pt_ite][0][0][0].GetXaxis().SetTitle("#p_{t}^{GEN} / GeV")

    width_response[pt_ite][0][0][0].GetXaxis().CenterTitle()
    width_response[pt_ite][0][0][0].GetYaxis().CenterTitle()
  
    width_response[pt_ite][0][0][0].GetXaxis().SetTitleOffset(1.1)
    width_response[pt_ite][0][0][0].GetYaxis().SetTitleOffset(1.1)
    #width_response[pt_ite][0][0][0].GetYaxis().SetRangeUser(0.8,1.25)
	
    color = my_colorsP[0]
    width_response[pt_ite][0][0][0].SetLineColor(color)
    width_response[pt_ite][0][0][0].SetLineWidth(2)
    width_response[pt_ite][0][0][0].SetMarkerStyle(20)
    width_response[pt_ite][0][0][0].SetMarkerSize(0.8)
    width_response[pt_ite][0][0][0].SetMarkerColor(color)

    width_response[pt_ite][0][0][0].Draw()
	
    color = my_colorsP[1]
    width_response[pt_ite][0][1][0].SetLineColor(color)
    width_response[pt_ite][0][1][0].SetLineWidth(2)
    width_response[pt_ite][0][1][0].SetMarkerStyle(20)
    width_response[pt_ite][0][1][0].SetMarkerSize(0.8)
    width_response[pt_ite][0][1][0].SetMarkerColor(color)

    width_response[pt_ite][0][1][0].Draw("same")
	
    color = my_colorsP[3]
    width_response[pt_ite][0][3][0].SetLineColor(color)
    width_response[pt_ite][0][3][0].SetLineWidth(2)
    width_response[pt_ite][0][3][0].SetMarkerStyle(20)
    width_response[pt_ite][0][3][0].SetMarkerSize(0.8)
    width_response[pt_ite][0][3][0].SetMarkerColor(color)

    width_response[pt_ite][0][3][0].Draw("same")
		    	    

width_response[0][0][0][1].SetTitle("WIDTH vs pt - Whole npu range")
width_response[1][0][0][1].SetTitle("WIDTH vs pt - npu < 15.")
width_response[2][0][0][1].SetTitle("WIDTH vs pt - 15. <= npu < 30.")
width_response[3][0][0][1].SetTitle("WIDTH vs pt - 30. <= npu")
		   
for pt_ite in range(npu_number):
    canv.cd(pt_ite+5)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetTopMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetFillColor(0)
    
    width_response[pt_ite][0][0][1].GetYaxis().SetTitle("#sigma_{p_{t}^{RECO}/p_{t}^{GEN}}")
    width_response[pt_ite][0][0][1].GetXaxis().SetTitle("#p_{t}^{GEN} / GeV")
    
    width_response[pt_ite][0][0][1].GetXaxis().CenterTitle()
    width_response[pt_ite][0][0][1].GetYaxis().CenterTitle()
  
    width_response[pt_ite][0][0][1].GetXaxis().SetTitleOffset(1.1)
    width_response[pt_ite][0][0][1].GetYaxis().SetTitleOffset(1.1)
    #width_response[pt_ite][0][0][1].GetYaxis().SetRangeUser(0.8,1.25)
	
    color = my_colorsP[0]
    width_response[pt_ite][0][0][1].SetLineColor(color)
    width_response[pt_ite][0][0][1].SetLineWidth(2)
    width_response[pt_ite][0][0][1].SetMarkerStyle(20)
    width_response[pt_ite][0][0][1].SetMarkerSize(0.8)
    width_response[pt_ite][0][0][1].SetMarkerColor(color)

    width_response[pt_ite][0][0][1].Draw()
	
    color = my_colorsP[1]
    width_response[pt_ite][0][1][1].SetLineColor(color)
    width_response[pt_ite][0][1][1].SetLineWidth(2)
    width_response[pt_ite][0][1][1].SetMarkerStyle(20)
    width_response[pt_ite][0][1][1].SetMarkerSize(0.8)
    width_response[pt_ite][0][1][1].SetMarkerColor(color)

    width_response[pt_ite][0][1][1].Draw("same")
	
    color = my_colorsP[3]
    width_response[pt_ite][0][3][1].SetLineColor(color)
    width_response[pt_ite][0][3][1].SetLineWidth(2)
    width_response[pt_ite][0][3][1].SetMarkerStyle(20)
    width_response[pt_ite][0][3][1].SetMarkerSize(0.8)
    width_response[pt_ite][0][3][1].SetMarkerColor(color)

    width_response[pt_ite][0][3][1].Draw("same")

canv.SaveAs("../pictures/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedNpu_Width.png")
canv.SaveAs("../root/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedNpu_Width.root")
canv.SaveAs("../pdfs/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedNpu_Width.pdf")
