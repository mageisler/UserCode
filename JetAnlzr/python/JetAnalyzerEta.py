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

pt_number = len(pt_ranges)
	
print ""

File_ref = [[]] * file_number

##the order is file - ass - dir - npu - pt 

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

                histo_response_ass_dir_npu_pt = () 
		
	        help = random.random()
	    
	        response_help = File_ref[file_ite].Get(subdirnames[dir_ite][ass_ite]+ npu_ranges[npu_ite])
	    
	        response_help.SetName(subdirnames[dir_ite][ass_ite] + npu_ranges[npu_ite] + str(help))
		
		histo_response_ass_dir_npu_pt+= response_help.Project3D("cR_zyo"),
	
	        for pt_ite in range(pt_number-1):	
	    
	            response_help.GetXaxis().SetRangeUser(pt_ranges[pt_ite],pt_ranges[pt_ite+1])	    
	            histo_response_ass_dir_npu_pt+= response_help.Project3D(str(pt_ranges[pt_ite+1]) + "_zyo"),
	            
	    
                histo_response_ass_dir_npu+= histo_response_ass_dir_npu_pt,
	    
            histo_response_ass_dir+= histo_response_ass_dir_npu,
	    
        histo_response_ass+= histo_response_ass_dir,
	
    histo_response_queue+= histo_response_ass,

    print "done" 
    
##the former order is file - ass - dir - npu - pt
##the new order is ass - dir - npu - pt

print ""
	    
histo_response_com = GetCombinedTH2FsWeight4D(histo_response_queue,filelist)

##the former order is ass - dir - npu - pt
##the new order is npu - pt - ass - dir 
    
mean_response = ()
width_response = ()

print "Profiling"
    
for npu_ite in range(npu_number):
    
    mean_response_npu = ()       
    width_response_npu = ()
	
    for pt_ite in range(pt_number):
    
        mean_response_npu_pt = ()       
        width_response_npu_pt = ()
	
        for ass_ite in range(ass_number):
    
            mean_response_npu_pt_ass = ()       
            width_response_npu_pt_ass = ()
	
            for dir_ite in range(dir_number):
	    
	        mean_response_npu_pt_ass+= GetMean2(histo_response_com[ass_ite][dir_ite][npu_ite][pt_ite]), 
	        width_response_npu_pt_ass+= GetWidth2(histo_response_com[ass_ite][dir_ite][npu_ite][pt_ite]),
    
            mean_response_npu_pt+= mean_response_npu_pt_ass,        
            width_response_npu_pt+= width_response_npu_pt_ass,
    
        mean_response_npu+= mean_response_npu_pt,        
        width_response_npu+= width_response_npu_pt,
    
    mean_response+= mean_response_npu,        
    width_response+= width_response_npu,
    
print "" 
	
for npu_ite in range(npu_number):
            
    for pt_ite in range(pt_number):
	    
        for dir_ite in range(dir_number):

            title = dirnames[dir_ite] + " - pt" + pt_names[pt_ite]
		    
	    mean_response[npu_ite][pt_ite][0][dir_ite].SetTitle("p_{t}-response - MEAN - " + title)
	    mean_response[npu_ite][pt_ite][0][dir_ite].GetXaxis().SetTitle("#eta^{GEN}")
	    mean_response[npu_ite][pt_ite][0][dir_ite].GetYaxis().SetTitle("< p_{t}^{RECO} / p_{t}^{GEN} >") 
	    mean_response[npu_ite][pt_ite][0][dir_ite].GetYaxis().SetRangeUser(0.7,1.3) 	
			    
	    width_response[npu_ite][pt_ite][0][dir_ite].SetTitle("p_{t}-response - WIDTH - " + title)
	    width_response[npu_ite][pt_ite][0][dir_ite].GetXaxis().SetTitle("#eta^{GEN}")
	    width_response[npu_ite][pt_ite][0][dir_ite].GetYaxis().SetTitle("#sigma_{p_{t}^{RECO}/p_{t}^{GEN}}") 
	

#for npu_ite in range(npu_number):
    #DrawTH1Fs(mean_response[npu_ite],"JetAnlzr_"+name+"_"+npu_ranges[npu_ite].split("_")[0]+"_Eta_"+npu_ranges[npu_ite].split("_")[1]+"_Mean",False,True,True)
    #DrawTH1Fs(width_response[npu_ite],"JetAnlzr_"+name+"_"+npu_ranges[npu_ite].split("_")[0]+"_Eta_"+npu_ranges[npu_ite].split("_")[1]+"_Width",True,True,True)


###for presentations
##the new order is npu - pt - ass - dir 
  
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
		    

mean_response[0][0][0][0].SetTitle("MEAN vs eta - Whole pt range")
mean_response[0][1][0][0].SetTitle("MEAN vs eta - 30. <= pt <= 100.")
mean_response[0][2][0][0].SetTitle("MEAN vs eta - 100. <= pt <= 300.")
mean_response[0][3][0][0].SetTitle("MEAN vs eta - 300. <= pt <= 1000.")
		   
for pt_ite in range(pt_number):
    canv.cd(pt_ite+1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetTopMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetFillColor(0)
    
    mean_response[0][pt_ite][0][0].GetYaxis().SetTitle("< p_{t}^{RECO} / p_{t}^{GEN} >")
    mean_response[0][pt_ite][0][0].GetXaxis().SetTitle("#eta^{GEN}")

    mean_response[0][pt_ite][0][0].GetXaxis().CenterTitle()
    mean_response[0][pt_ite][0][0].GetYaxis().CenterTitle()
  
    mean_response[0][pt_ite][0][0].GetXaxis().SetTitleOffset(1.1)
    mean_response[0][pt_ite][0][0].GetYaxis().SetTitleOffset(1.1)
    mean_response[0][pt_ite][0][0].GetYaxis().SetRangeUser(0.8,1.25)
	
    color = my_colorsP[0]
    mean_response[0][pt_ite][0][0].SetLineColor(color)
    mean_response[0][pt_ite][0][0].SetLineWidth(2)
    mean_response[0][pt_ite][0][0].SetMarkerStyle(20)
    mean_response[0][pt_ite][0][0].SetMarkerSize(0.8)
    mean_response[0][pt_ite][0][0].SetMarkerColor(color)

    mean_response[0][pt_ite][0][0].Draw()
	
    color = my_colorsP[1]
    mean_response[0][pt_ite][1][0].SetLineColor(color)
    mean_response[0][pt_ite][1][0].SetLineWidth(2)
    mean_response[0][pt_ite][1][0].SetMarkerStyle(20)
    mean_response[0][pt_ite][1][0].SetMarkerSize(0.8)
    mean_response[0][pt_ite][1][0].SetMarkerColor(color)

    mean_response[0][pt_ite][1][0].Draw("same")
	
    color = my_colorsP[3]
    mean_response[0][pt_ite][3][0].SetLineColor(color)
    mean_response[0][pt_ite][3][0].SetLineWidth(2)
    mean_response[0][pt_ite][3][0].SetMarkerStyle(20)
    mean_response[0][pt_ite][3][0].SetMarkerSize(0.8)
    mean_response[0][pt_ite][3][0].SetMarkerColor(color)

    mean_response[0][pt_ite][3][0].Draw("same")
		    

mean_response[0][0][0][1].SetTitle("MEAN vs eta - Whole pt range")
mean_response[0][1][0][1].SetTitle("MEAN vs eta - 30. <= pt <= 100.")
mean_response[0][2][0][1].SetTitle("MEAN vs eta - 100. <= pt <= 300.")
mean_response[0][3][0][1].SetTitle("MEAN vs eta - 300. <= pt <= 1000.")
		   
for pt_ite in range(pt_number):
    canv.cd(pt_ite+5)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetTopMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetFillColor(0)
    
    mean_response[0][pt_ite][0][1].GetYaxis().SetTitle("< p_{t}^{RECO} / p_{t}^{GEN} >")
    mean_response[0][pt_ite][0][1].GetXaxis().SetTitle("#eta^{GEN}")
    
    mean_response[0][pt_ite][0][1].GetXaxis().CenterTitle()
    mean_response[0][pt_ite][0][1].GetYaxis().CenterTitle()
  
    mean_response[0][pt_ite][0][1].GetXaxis().SetTitleOffset(1.1)
    mean_response[0][pt_ite][0][1].GetYaxis().SetTitleOffset(1.1)
    mean_response[0][pt_ite][0][1].GetYaxis().SetRangeUser(0.8,1.25)
	
    color = my_colorsP[0]
    mean_response[0][pt_ite][0][1].SetLineColor(color)
    mean_response[0][pt_ite][0][1].SetLineWidth(2)
    mean_response[0][pt_ite][0][1].SetMarkerStyle(20)
    mean_response[0][pt_ite][0][1].SetMarkerSize(0.8)
    mean_response[0][pt_ite][0][1].SetMarkerColor(color)

    mean_response[0][pt_ite][0][1].Draw()
	
    color = my_colorsP[1]
    mean_response[0][pt_ite][1][1].SetLineColor(color)
    mean_response[0][pt_ite][1][1].SetLineWidth(2)
    mean_response[0][pt_ite][1][1].SetMarkerStyle(20)
    mean_response[0][pt_ite][1][1].SetMarkerSize(0.8)
    mean_response[0][pt_ite][1][1].SetMarkerColor(color)

    mean_response[0][pt_ite][1][1].Draw("same")
	
    color = my_colorsP[3]
    mean_response[0][pt_ite][3][1].SetLineColor(color)
    mean_response[0][pt_ite][3][1].SetLineWidth(2)
    mean_response[0][pt_ite][3][1].SetMarkerStyle(20)
    mean_response[0][pt_ite][3][1].SetMarkerSize(0.8)
    mean_response[0][pt_ite][3][1].SetMarkerColor(color)

    mean_response[0][pt_ite][3][1].Draw("same")

canv.SaveAs("../pictures/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedPt_Mean.png")
canv.SaveAs("../root/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedPt_Mean.root")
canv.SaveAs("../pdfs/JetAnlzr_"+name+"_PtResponse_Presentation_ChangedPt_Mean.pdf")
