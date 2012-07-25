#!/usr/bin/env python

import sys, math, array, getopt, random, ROOT
from MGeisler.JetAnlzr.AnalysisTools import *


#_____________________
#
# READ INPUT PARAMETERS
#_____________________

letters = 'p:'
keywords = ['path']
opts, extraparams = getopt.getopt(sys.argv[1:], letters, keywords)

selection = "LooseTP"


path=""

for o,p in opts: 
    if o in ['-p','--path']:
        path = p
	
if "Tight" in path:
    selection = "TightTP"


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
    
x_axis_list = ("_eta","_pt","_npu")
x_axis_number = len(x_axis_list)
    

#_____________________
#
# READ INPUT FILES
#_____________________

dirnames = GetListOfSubdirectories(filelist[0],"trackValidatorWO/","Track")
ass_number = len(dirnames)

histo_num_simul_tracks_queue = ()
histo_num_assoc_queue = ()

histo_num_reco_tracks_queue = ()
histo_num_assoc2_queue = ()

histo_num_removed_reco_PU_queue = ()
histo_num_track_reco_PU_queue = ()

histo_num_removed_reco_signal_queue = ()
histo_num_track_reco_signal_queue = ()

histo_effic2D_num_queue = ()
histo_effic2D_denum_queue = ()

histo_fakerate2D_num_queue = ()
histo_fakerate2D_denum_queue = ()

File_ref = [[]] * file_number
	
for file_ite in range(file_number):

    File_ref[file_ite] = ROOT.TFile.Open(filelist[file_ite])
    
    print "File " +str(file_ite)+ ": " + str(filelist[file_ite]) + " ...",
    
    histo_num_simul_tracks_ass = ()
    histo_num_assoc_ass = ()
    
    histo_num_reco_tracks_ass = ()
    histo_num_assoc2_ass = ()
    
    histo_num_removed_reco_PU_ass = ()
    histo_num_track_reco_PU_ass = ()
    
    histo_num_removed_reco_signal_ass = ()
    histo_num_track_reco_signal_ass = ()
    
    histo_effic2D_num_ass = ()
    histo_effic2D_denum_ass = ()
    
    histo_fakerate2D_num_ass = ()
    histo_fakerate2D_denum_ass = ()
    
    for ass_ite in range(ass_number):
    
        histo_num_simul_tracks_ass_axis = ()
        histo_num_assoc_ass_axis = ()
    
        histo_num_reco_tracks_ass_axis = ()
        histo_num_assoc2_ass_axis = ()
    
        histo_num_removed_reco_PU_ass_axis = ()
        histo_num_track_reco_PU_ass_axis = ()
    
        histo_num_removed_reco_signal_ass_axis = ()
        histo_num_track_reco_signal_ass_axis = ()
	    
        for x_axis_ite in range(x_axis_number):
		
	    help = int(random.random()*10000)
		 
            histo_num_simul_tracks_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_track_simul"+x_axis_list[x_axis_ite]),
            histo_num_assoc_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_assoc(simToReco)"+x_axis_list[x_axis_ite]),
	    
	    histo_num_simul_tracks_ass_axis[x_axis_ite].SetName(histo_num_simul_tracks_ass_axis[x_axis_ite].GetName() + str(help))
	    histo_num_assoc_ass_axis[x_axis_ite].SetName(histo_num_assoc_ass_axis[x_axis_ite].GetName() + str(help))
	    
            histo_num_reco_tracks_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_track_reco"+x_axis_list[x_axis_ite]),
            histo_num_assoc2_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_assoc(recoToSim)"+x_axis_list[x_axis_ite]),
	    
	    histo_num_reco_tracks_ass_axis[x_axis_ite].SetName(histo_num_reco_tracks_ass_axis[x_axis_ite].GetName() + str(help))
	    histo_num_assoc2_ass_axis[x_axis_ite].SetName(histo_num_assoc2_ass_axis[x_axis_ite].GetName() + str(help))
		 
            histo_num_removed_reco_PU_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_removed_reco_PU"+x_axis_list[x_axis_ite]),
            histo_num_track_reco_PU_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_track_reco_PU"+x_axis_list[x_axis_ite]),
	    
	    histo_num_removed_reco_PU_ass_axis[x_axis_ite].SetName(histo_num_removed_reco_PU_ass_axis[x_axis_ite].GetName() + str(help))
	    histo_num_track_reco_PU_ass_axis[x_axis_ite].SetName(histo_num_track_reco_PU_ass_axis[x_axis_ite].GetName() + str(help))
		 
            histo_num_removed_reco_signal_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_removed_reco_signal"+x_axis_list[x_axis_ite]),
            histo_num_track_reco_signal_ass_axis+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_track_reco_signal"+x_axis_list[x_axis_ite]),
	    
	    histo_num_removed_reco_signal_ass_axis[x_axis_ite].SetName(histo_num_removed_reco_signal_ass_axis[x_axis_ite].GetName() + str(help))
	    histo_num_track_reco_signal_ass_axis[x_axis_ite].SetName(histo_num_track_reco_signal_ass_axis[x_axis_ite].GetName() + str(help))
	    
	histo_num_simul_tracks_ass+=histo_num_simul_tracks_ass_axis,
	histo_num_assoc_ass+=histo_num_assoc_ass_axis,
	    
	histo_num_reco_tracks_ass+=histo_num_reco_tracks_ass_axis,
	histo_num_assoc2_ass+=histo_num_assoc2_ass_axis,
	    
	histo_num_removed_reco_PU_ass+=histo_num_removed_reco_PU_ass_axis,
	histo_num_track_reco_PU_ass+=histo_num_track_reco_PU_ass_axis,
	    
	histo_num_removed_reco_signal_ass+=histo_num_removed_reco_signal_ass_axis,
	histo_num_track_reco_signal_ass+=histo_num_track_reco_signal_ass_axis,
		
	help = int(random.random()*10000)
		 
        histo_effic2D_num_ass+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_assoc(simToReco)_npu_Contr"),
        histo_effic2D_denum_ass+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_simul_tracks_npu_Contr"),
	    
	histo_effic2D_num_ass[ass_ite].SetName(histo_effic2D_num_ass[ass_ite].GetName() + str(help))
	histo_effic2D_denum_ass[ass_ite].SetName(histo_effic2D_denum_ass[ass_ite].GetName() + str(help))
	
        histo_fakerate2D_num_ass+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_assoc(recoToSim)_npu_Contr"),
        histo_fakerate2D_denum_ass+= File_ref[file_ite].Get(dirnames[ass_ite]+"num_reco_tracks_npu_Contr"),
	    
	histo_fakerate2D_num_ass[ass_ite].SetName(histo_fakerate2D_num_ass[ass_ite].GetName() + str(help))
	histo_fakerate2D_denum_ass[ass_ite].SetName(histo_fakerate2D_denum_ass[ass_ite].GetName() + str(help))
	
    histo_num_simul_tracks_queue+=histo_num_simul_tracks_ass,
    histo_num_assoc_queue+=histo_num_assoc_ass,
	
    histo_num_reco_tracks_queue+=histo_num_reco_tracks_ass,
    histo_num_assoc2_queue+=histo_num_assoc2_ass,
	
    histo_num_removed_reco_PU_queue+=histo_num_removed_reco_PU_ass,
    histo_num_track_reco_PU_queue+=histo_num_track_reco_PU_ass,
	
    histo_num_removed_reco_signal_queue+=histo_num_removed_reco_signal_ass,
    histo_num_track_reco_signal_queue+=histo_num_track_reco_signal_ass,
	
    histo_effic2D_num_queue+=histo_effic2D_num_ass,
    histo_effic2D_denum_queue+=histo_effic2D_denum_ass,
    
    histo_fakerate2D_num_queue+=histo_fakerate2D_num_ass,
    histo_fakerate2D_denum_queue+=histo_fakerate2D_denum_ass,

    print "done"
    
print ""

histo_num_simul_tracks_com = GetCombinedTH1FsWeight(histo_num_simul_tracks_queue,filelist)
histo_num_assoc_com = GetCombinedTH1FsWeight(histo_num_assoc_queue,filelist)

histo_num_reco_tracks_com = GetCombinedTH1FsWeight(histo_num_reco_tracks_queue,filelist)
histo_num_assoc2_com = GetCombinedTH1FsWeight(histo_num_assoc2_queue,filelist)

histo_num_removed_reco_PU_com = GetCombinedTH1FsWeight(histo_num_removed_reco_PU_queue,filelist)
histo_num_track_reco_PU_com = GetCombinedTH1FsWeight(histo_num_track_reco_PU_queue,filelist)

histo_num_removed_reco_signal_com = GetCombinedTH1FsWeight(histo_num_removed_reco_signal_queue,filelist)
histo_num_track_reco_signal_com = GetCombinedTH1FsWeight(histo_num_track_reco_signal_queue,filelist)

histo_effic2D_num_com = GetCombinedTH2FsWeight(histo_effic2D_num_queue,filelist)
histo_effic2D_denum_com = GetCombinedTH2FsWeight(histo_effic2D_denum_queue,filelist)

histo_fakerate2D_num_com = GetCombinedTH2FsWeight(histo_fakerate2D_num_queue,filelist)
histo_fakerate2D_denum_com = GetCombinedTH2FsWeight(histo_fakerate2D_denum_queue,filelist)

## print the overall fractions
for ass_ite in range(ass_number):
	
    histo_effic_num_entries = histo_num_assoc_com[ass_ite][0].Integral() 
    histo_effic_denum_entries = histo_num_simul_tracks_com[ass_ite][0].Integral()
	
    histo_fakerate_num_entries = histo_num_assoc2_com[ass_ite][0].Integral() 
    histo_fakerate_denum_entries = histo_num_reco_tracks_com[ass_ite][0].Integral()
	
    histo_PU_effic_num_entries = histo_num_removed_reco_PU_com[ass_ite][0].Integral() 
    histo_PU_effic_denum_entries = histo_num_track_reco_PU_com[ass_ite][0].Integral()
	
    histo_PU_fakerate_num_entries = histo_num_removed_reco_signal_com[ass_ite][0].Integral() 
    histo_PU_fakerate_denum_entries = histo_num_track_reco_signal_com[ass_ite][0].Integral()
      
    effic = histo_effic_num_entries *1./histo_effic_denum_entries
    fakerate = 1.- histo_fakerate_num_entries *1./histo_fakerate_denum_entries
      
    PU_effic = histo_PU_effic_num_entries *1./histo_PU_effic_denum_entries    
    PU_fakerate = histo_PU_fakerate_num_entries *1./histo_PU_fakerate_denum_entries
    
    frac = effic *1./fakerate
    PU_frac = PU_effic *1./PU_fakerate
      
    print ""  
    print "Version: " + dirnames[ass_ite]
    print "Effic: " + str(effic)
    print "Fake rate: " + str(fakerate)
    print ""
    print "PU Effic: " + str(PU_effic)
    print "PU fake rate: " + str(PU_fakerate)
    print ""
    print "Fraction: " + str(frac)
    print "PU fraction: " + str(PU_frac)
    print ""
	
print "Calculating fractions...",

histo_effic = GetFractionTH1Fs(histo_num_assoc_com,histo_num_simul_tracks_com,"effic")
histo_fakerate = GetFractionTH1Fs(histo_num_assoc2_com,histo_num_reco_tracks_com,"fakerate")

histo_PU_effic = GetFractionTH1Fs(histo_num_removed_reco_PU_com,histo_num_track_reco_PU_com,"effic")
histo_PU_fakerate = GetFractionTH1Fs(histo_num_removed_reco_signal_com,histo_num_track_reco_signal_com,"effic")

histo_effic2D = GetFractionTH2Fs(histo_effic2D_num_com,histo_effic2D_denum_com,"effic")
histo_fakerate2D = GetFractionTH2Fs(histo_fakerate2D_num_com,histo_fakerate2D_denum_com,"fakerate")

print "done"
print ""

## create and draw the canvases

histo_effic[0][0].SetTitle("efficiency vs #eta")
histo_effic[0][0].GetXaxis().SetTitle("#eta^{GEN}")
histo_effic[0][0].GetYaxis().SetTitle("efficiency")
histo_effic[0][0].GetYaxis().SetRangeUser(0.,1.05)
histo_effic[0][1].SetTitle("efficiency vs pt")
histo_effic[0][1].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
histo_effic[0][1].GetYaxis().SetTitle("efficiency")
histo_effic[0][1].GetYaxis().SetRangeUser(0.,1.05)
histo_effic[0][2].SetTitle("efficiency vs npu")
histo_effic[0][2].GetXaxis().SetTitle("number of pileup interactions")
histo_effic[0][2].GetYaxis().SetTitle("efficiency")
histo_effic[0][2].GetYaxis().SetRangeUser(0.,1.05)

histo_fakerate[0][0].SetTitle("fake rate vs #eta")
histo_fakerate[0][0].GetXaxis().SetTitle("#eta^{GEN}")
histo_fakerate[0][0].GetYaxis().SetTitle("fake rate")
histo_fakerate[0][0].GetYaxis().SetRangeUser(0.,1.05)
histo_fakerate[0][1].SetTitle("fake rate vs pt")
histo_fakerate[0][1].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
histo_fakerate[0][1].GetYaxis().SetTitle("fake rate")
histo_fakerate[0][1].GetYaxis().SetRangeUser(0.,1.05)
histo_fakerate[0][2].SetTitle("fake rate vs npu")
histo_fakerate[0][2].GetXaxis().SetTitle("number of pileup interactions")
histo_fakerate[0][2].GetYaxis().SetTitle("fake rate")
histo_fakerate[0][2].GetYaxis().SetRangeUser(0.,1.05)
	 
all_histos = ()
all_histos+=histo_effic,
all_histos+=histo_fakerate,

#DrawTH1Fs(all_histos,"FakeEff_JetAnlzr_" + name + "_" + str(selection),False,False,True)


##PU
histo_PU_effic[0][0].SetTitle("PU efficiency vs #eta")
histo_PU_effic[0][0].GetXaxis().SetTitle("#eta^{GEN}")
histo_PU_effic[0][0].GetYaxis().SetTitle("PU efficiency")
histo_PU_effic[0][0].GetYaxis().SetRangeUser(0.,1.05)
histo_PU_effic[0][1].SetTitle("PU efficiency vs pt")
histo_PU_effic[0][1].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
histo_PU_effic[0][1].GetYaxis().SetTitle("PU efficiency")
histo_PU_effic[0][1].GetYaxis().SetRangeUser(0.,1.05)
histo_PU_effic[0][2].SetTitle("PU efficiency vs npu")
histo_PU_effic[0][2].GetXaxis().SetTitle("number of pileup interactions")
histo_PU_effic[0][2].GetYaxis().SetTitle("PU efficiency")
histo_PU_effic[0][2].GetYaxis().SetRangeUser(0.,1.05)

histo_PU_fakerate[0][0].SetTitle("PU fake rate vs #eta")
histo_PU_fakerate[0][0].GetXaxis().SetTitle("#eta^{GEN}")
histo_PU_fakerate[0][0].GetYaxis().SetTitle("PU fake rate")
histo_PU_fakerate[0][0].GetYaxis().SetRangeUser(0.,1.05)
histo_PU_fakerate[0][1].SetTitle("PU fake rate vs pt")
histo_PU_fakerate[0][1].GetXaxis().SetTitle("p_{t}^{GEN} / GeV")
histo_PU_fakerate[0][1].GetYaxis().SetTitle("PU fake rate")
histo_PU_fakerate[0][1].GetYaxis().SetRangeUser(0.,1.05)
histo_PU_fakerate[0][2].SetTitle("PU fake rate vs npu")
histo_PU_fakerate[0][2].GetXaxis().SetTitle("number of pileup interactions")
histo_PU_fakerate[0][2].GetYaxis().SetTitle("PU fake rate")
histo_PU_fakerate[0][2].GetYaxis().SetRangeUser(0.,1.05)

all_PU_histos = ()
all_PU_histos+=histo_PU_effic,
all_PU_histos+=histo_PU_fakerate,

#DrawTH1Fs(all_PU_histos,"FakeEff_JetAnlzr_PU_" + name + "_" + str(selection),False,False,True)
  
  
##2D
histo_effic2D[0].SetTitle("efficiency vs npu vs contribution")
histo_effic2D[0].GetXaxis().SetTitle("number of pileup interactions")
histo_effic2D[0].GetYaxis().SetTitle("contribution")
histo_effic2D[0].GetZaxis().SetTitle("efficiency")
histo_effic2D[0].GetZaxis().SetRangeUser(0.,1.05)
histo_effic2D[1].SetTitle("efficiency vs npu vs contribution")
histo_effic2D[1].GetXaxis().SetTitle("number of pileup interactions")
histo_effic2D[1].GetYaxis().SetTitle("contribution")
histo_effic2D[1].GetZaxis().SetTitle("efficiency")
histo_effic2D[1].GetZaxis().SetRangeUser(0.,1.05)
histo_effic2D[2].SetTitle("efficiency vs npu vs contribution")
histo_effic2D[2].GetXaxis().SetTitle("number of pileup interactions")
histo_effic2D[2].GetYaxis().SetTitle("contribution")
histo_effic2D[2].GetZaxis().SetTitle("efficiency")
histo_effic2D[2].GetZaxis().SetRangeUser(0.,1.05)

histo_fakerate2D[0].SetTitle("fake rate vs npu vs contribution")
histo_fakerate2D[0].GetXaxis().SetTitle("number of pileup interactions")
histo_fakerate2D[0].GetYaxis().SetTitle("contribution")
histo_fakerate2D[0].GetZaxis().SetTitle("fake rate")
histo_fakerate2D[0].GetZaxis().SetRangeUser(0.,1.05)
histo_fakerate2D[1].SetTitle("fake rate vs npu vs contribution")
histo_fakerate2D[1].GetXaxis().SetTitle("number of pileup interactions")
histo_fakerate2D[1].GetYaxis().SetTitle("contribution")
histo_fakerate2D[1].GetZaxis().SetTitle("fake rate")
histo_fakerate2D[1].GetZaxis().SetRangeUser(0.,1.05)
histo_fakerate2D[2].SetTitle("fake rate vs npu vs contribution")
histo_fakerate2D[2].GetXaxis().SetTitle("number of pileup interactions")
histo_fakerate2D[2].GetYaxis().SetTitle("contribution")
histo_fakerate2D[2].GetZaxis().SetTitle("fake rate")
histo_fakerate2D[2].GetZaxis().SetRangeUser(0.,1.05)

all_2D_histos = ()
all_2D_histos+=histo_effic2D,
all_2D_histos+=histo_fakerate2D,

#DrawTH2Fs(all_2D_histos,"FakeEff_JetAnlzr_2D_" + name + "_" + str(selection),"COLZ")


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

histo_effic[0][0].GetXaxis().CenterTitle()
histo_effic[0][0].GetYaxis().CenterTitle()
  
histo_effic[0][0].GetXaxis().SetTitleOffset(1.1)
histo_effic[0][0].GetYaxis().SetTitleOffset(1.1)
histo_effic[0][0].GetYaxis().SetRangeUser(0.3,0.7)
	
color = my_colors[0]
histo_effic[0][0].SetLineColor(color)
histo_effic[0][0].SetLineWidth(2)
histo_effic[0][0].SetMarkerStyle(20)
histo_effic[0][0].SetMarkerSize(0.8)
histo_effic[0][0].SetMarkerColor(color)

histo_effic[0][0].Draw()


color = my_colors[1]
histo_effic[1][0].SetLineColor(color)
histo_effic[1][0].SetLineWidth(2)		
histo_effic[1][0].SetMarkerStyle(20)
histo_effic[1][0].SetMarkerSize(0.8)
histo_effic[1][0].SetMarkerColor(color)

histo_effic[1][0].Draw("same")


color = my_colors[2]
histo_effic[2][0].SetLineColor(color)
histo_effic[2][0].SetLineWidth(2)		
histo_effic[2][0].SetMarkerStyle(20)
histo_effic[2][0].SetMarkerSize(0.8)
histo_effic[2][0].SetMarkerColor(color)

histo_effic[2][0].Draw("same")


color = my_colors[3]
histo_effic[3][0].SetLineColor(color)
histo_effic[3][0].SetLineWidth(2)		
histo_effic[3][0].SetMarkerStyle(20)
histo_effic[3][0].SetMarkerSize(0.8)
histo_effic[3][0].SetMarkerColor(color)

histo_effic[3][0].Draw("same")
		    
canv.cd(2)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

histo_PU_effic[0][0].GetXaxis().CenterTitle()
histo_PU_effic[0][0].GetYaxis().CenterTitle()
  
histo_PU_effic[0][0].GetXaxis().SetTitleOffset(1.1)
histo_PU_effic[0][0].GetYaxis().SetTitleOffset(1.1)
histo_PU_effic[0][0].GetYaxis().SetRangeUser(0.5,1.05)
	
color = my_colors[0]
histo_PU_effic[0][0].SetLineColor(color)
histo_PU_effic[0][0].SetLineWidth(2)
histo_PU_effic[0][0].SetMarkerStyle(20)
histo_PU_effic[0][0].SetMarkerSize(0.8)
histo_PU_effic[0][0].SetMarkerColor(color)

histo_PU_effic[0][0].Draw()


color = my_colors[1]
histo_PU_effic[1][0].SetLineColor(color)
histo_PU_effic[1][0].SetLineWidth(2)		
histo_PU_effic[1][0].SetMarkerStyle(20)
histo_PU_effic[1][0].SetMarkerSize(0.8)
histo_PU_effic[1][0].SetMarkerColor(color)

histo_PU_effic[1][0].Draw("same")


color = my_colors[2]
histo_PU_effic[2][0].SetLineColor(color)
histo_PU_effic[2][0].SetLineWidth(2)		
histo_PU_effic[2][0].SetMarkerStyle(20)
histo_PU_effic[2][0].SetMarkerSize(0.8)
histo_PU_effic[2][0].SetMarkerColor(color)

histo_PU_effic[2][0].Draw("same")


color = my_colors[3]
histo_PU_effic[3][0].SetLineColor(color)
histo_PU_effic[3][0].SetLineWidth(2)		
histo_PU_effic[3][0].SetMarkerStyle(20)
histo_PU_effic[3][0].SetMarkerSize(0.8)
histo_PU_effic[3][0].SetMarkerColor(color)

histo_PU_effic[3][0].Draw("same")
		    
canv.cd(3)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

histo_fakerate[0][0].GetXaxis().CenterTitle()
histo_fakerate[0][0].GetYaxis().CenterTitle()
  
histo_fakerate[0][0].GetXaxis().SetTitleOffset(1.1)
histo_fakerate[0][0].GetYaxis().SetTitleOffset(1.1)
histo_fakerate[0][0].GetYaxis().SetRangeUser(0.,1.05)
	
color = my_colors[0]
histo_fakerate[0][0].SetLineColor(color)
histo_fakerate[0][0].SetLineWidth(2)
histo_fakerate[0][0].SetMarkerStyle(20)
histo_fakerate[0][0].SetMarkerSize(0.8)
histo_fakerate[0][0].SetMarkerColor(color)

histo_fakerate[0][0].Draw()


color = my_colors[1]
histo_fakerate[1][0].SetLineColor(color)
histo_fakerate[1][0].SetLineWidth(2)		
histo_fakerate[1][0].SetMarkerStyle(20)
histo_fakerate[1][0].SetMarkerSize(0.8)
histo_fakerate[1][0].SetMarkerColor(color)

histo_fakerate[1][0].Draw("same")


color = my_colors[2]
histo_fakerate[2][0].SetLineColor(color)
histo_fakerate[2][0].SetLineWidth(2)		
histo_fakerate[2][0].SetMarkerStyle(20)
histo_fakerate[2][0].SetMarkerSize(0.8)
histo_fakerate[2][0].SetMarkerColor(color)

histo_fakerate[2][0].Draw("same")


color = my_colors[3]
histo_fakerate[3][0].SetLineColor(color)
histo_fakerate[3][0].SetLineWidth(2)		
histo_fakerate[3][0].SetMarkerStyle(20)
histo_fakerate[3][0].SetMarkerSize(0.8)
histo_fakerate[3][0].SetMarkerColor(color)

histo_fakerate[3][0].Draw("same")
		    
canv.cd(4)
ROOT.gPad.SetGrid()
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetLeftMargin(0.2)
ROOT.gPad.SetTopMargin(0.15)
ROOT.gPad.SetBottomMargin(0.15)
ROOT.gPad.SetFillColor(0)

histo_PU_fakerate[0][0].GetXaxis().CenterTitle()
histo_PU_fakerate[0][0].GetYaxis().CenterTitle()
  
histo_PU_fakerate[0][0].GetXaxis().SetTitleOffset(1.1)
histo_PU_fakerate[0][0].GetYaxis().SetTitleOffset(1.1)
histo_PU_fakerate[0][0].GetYaxis().SetRangeUser(0.,0.3)
	
color = my_colors[0]
histo_PU_fakerate[0][0].SetLineColor(color)
histo_PU_fakerate[0][0].SetLineWidth(2)
histo_PU_fakerate[0][0].SetMarkerStyle(20)
histo_PU_fakerate[0][0].SetMarkerSize(0.8)
histo_PU_fakerate[0][0].SetMarkerColor(color)

histo_PU_fakerate[0][0].Draw()


color = my_colors[1]
histo_PU_fakerate[1][0].SetLineColor(color)
histo_PU_fakerate[1][0].SetLineWidth(2)		
histo_PU_fakerate[1][0].SetMarkerStyle(20)
histo_PU_fakerate[1][0].SetMarkerSize(0.8)
histo_PU_fakerate[1][0].SetMarkerColor(color)

histo_PU_fakerate[1][0].Draw("same")


color = my_colors[2]
histo_PU_fakerate[2][0].SetLineColor(color)
histo_PU_fakerate[2][0].SetLineWidth(2)		
histo_PU_fakerate[2][0].SetMarkerStyle(20)
histo_PU_fakerate[2][0].SetMarkerSize(0.8)
histo_PU_fakerate[2][0].SetMarkerColor(color)

histo_PU_fakerate[2][0].Draw("same")


color = my_colors[3]
histo_PU_fakerate[3][0].SetLineColor(color)
histo_PU_fakerate[3][0].SetLineWidth(2)		
histo_PU_fakerate[3][0].SetMarkerStyle(20)
histo_PU_fakerate[3][0].SetMarkerSize(0.8)
histo_PU_fakerate[3][0].SetMarkerColor(color)

histo_PU_fakerate[3][0].Draw("same")

canv.SaveAs("../pictures/FakeEff_Presentation" + name + "_WithOutWeight.png")
canv.SaveAs("../root/FakeEff_Presentation" + name + "_WithOutWeight.root")
