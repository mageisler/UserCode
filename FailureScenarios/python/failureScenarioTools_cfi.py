#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms
import sys, os, math, array, ROOT
from MGeisler.TrackerPlots.sourceFiles_cfi import *


ROOT.gROOT.Macro( os.path.expanduser( '~/Phd-Study/CMSSW/Helpers/SetPrettyStyle.C' ) )

input='/user/geisler/FS2011/FileLists/InputFileList.txt'
dictOfLists = ReadProdDetails(input)

my_colors = dictOfLists["colour"]

errGraphs = [[] for i in range(5)]
efficRatios = [[] for i in range(4)]
numRatios = [[] for i in range(4)]

def moveStatBox1(canvas,histo,color,step):

	canvas.Update()
	statBox = histo.GetListOfFunctions().FindObject("stats")
	
	statBoxSpacing = 0.005
	
	statBox.SetLineWidth(2)
	statBox.SetLineColor(color)

	size_y = 0.06
	
	statBox.SetX1NDC(0.70)
	statBox.SetX2NDC(0.95)
	statBox.SetY2NDC(0.88 - step*(size_y+statBoxSpacing))
	statBox.SetY1NDC(0.88 - (step+1)*size_y-step*statBoxSpacing)

	statBox.Draw()

def moveStatBox2(canvas,histo,color,step):

	canvas.Update()
	statBox = histo.GetListOfFunctions().FindObject("stats")
	
	statBoxSpacing = 0.005
	
	statBox.SetLineWidth(2)
	statBox.SetLineColor(color)

	size_y = 0.13
	
	statBox.SetX1NDC(0.70)
	statBox.SetX2NDC(0.95)
	statBox.SetY2NDC(0.88 - step*(size_y+statBoxSpacing))
	statBox.SetY1NDC(0.88 - (step+1)*size_y-step*statBoxSpacing)

	statBox.Draw()

def GetRatio(num,denum):
    
    bin_num = num.GetNbinsX()
    x_min = num.GetBinLowEdge(1)
    x_max = num.GetBinLowEdge(bin_num+1)
	
    histo_ratio = ROOT.TH1F("histo_ratio","; #eta; ratio",bin_num,x_min,x_max)		
    histo_ratio.Sumw2()
    
    for x_ite in range(1,bin_num+1):
        if not denum.GetBinContent(x_ite)==0:
	    N = num.GetBinContent(x_ite)
	    D = denum.GetBinContent(x_ite)
	    frac = N*1./D
	    err_N = num.GetBinError(x_ite)
	    err_D = denum.GetBinError(x_ite)
	    err = math.sqrt((err_N*err_N/(D*D)) +(N*N*err_D*err_D/(D*D*D*D)))
	    histo_ratio.SetBinContent(x_ite,frac)
	    histo_ratio.SetBinError(x_ite,err)
	else:
	    histo_ratio.SetBinContent(x_ite,0.)
	    histo_ratio.SetBinError(x_ite,0.01)
	    
    return histo_ratio
       
def CreateFractionOfValidHits(scenario):
    
    ref_file = ROOT.TFile.Open("../files/FS2012_RECO-Tracks_Draft.root")
    fs_file = ROOT.TFile.Open("../files/FS2012_RECO-Tracks_" + str(scenario) + ".root")
    
    ref_prof = ref_file.Get("Trackhitcount/RECO-Tracks_Draft_TrackerHitProfile")
    fs_prof = fs_file.Get("Trackhitcount/RECO-Tracks_" + str(scenario) + "_TrackerHitProfile")
    
    x_bins_num = fs_prof.GetNbinsX()
    x_bins_low = fs_prof.GetBinLowEdge(1)
    x_bins_high = fs_prof.GetBinLowEdge(x_bins_num+1)
    
    y_bins_num = fs_prof.GetNbinsY()
    y_bins_low = fs_prof.GetYaxis().GetBinLowEdge(1)
    y_bins_high = fs_prof.GetYaxis().GetBinLowEdge(y_bins_num+1)
  
    frac_prof = ROOT.TH2F("frac_prof"+str(scenario),"Fraction of layers hit vs #phi vs #eta; #phi; #eta",x_bins_num,x_bins_low,x_bins_high,y_bins_num,y_bins_low,y_bins_high)
    
    print " Created empty frac_prof" +str(scenario)
    
    for x_ite in range(1,x_bins_num+1):
        for y_ite in range(1,y_bins_num+1):
	    if not ref_prof.GetBinContent(x_ite,y_ite)==0:
		frac = fs_prof.GetBinContent(x_ite,y_ite)/ref_prof.GetBinContent(x_ite,y_ite)
	        frac_prof.SetBinContent(x_ite,y_ite,frac)
    
    print " Filled frac_prof" +str(scenario)
  
    frac_canv = ROOT.TCanvas("frac_canv","Guck mal",10,10,1200,850) 
    ROOT.SetRootPalette(2)
    
    frac_canv.Divide(1,1)
  
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetBorderMode(0)
    ROOT.gPad.SetFillColor(0)
   
    frac_canv.cd(1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetRightMargin(0.12)
    ROOT.gPad.SetLeftMargin(0.08)
    ROOT.gPad.SetFillColor(0)

    frac_prof.GetXaxis().CenterTitle()
    frac_prof.GetYaxis().CenterTitle()
    frac_prof.GetZaxis().CenterTitle()
    #frac_prof.GetZaxis().SetRangeUser(0.4,1.1)
    frac_prof.SetFillColor(0)
    frac_prof.SetMaximum(1.)    
    
    frac_prof.GetXaxis().SetTitleOffset(1.1)
    frac_prof.GetYaxis().SetTitleOffset(0.8)
    frac_prof.GetZaxis().SetTitleOffset(1.08)
    frac_prof.SetStats(0)
    frac_prof.Draw("COLZ")
    
    frac_canv.SaveAs("../pictures/RECO-Tracks_" + str(scenario) + "_TrackerHitFraction.png")
    
def CreateTrackingEfficiencyPlots(scenario,WithRef):
	
    fileName = "../files/FS2012_RECO-Tracks_" + str(scenario) + "_TV.root"
    fs_file = ROOT.TFile.Open(fileName)
    
    effic_cR = fs_file.Get("DQMData/Tracking/Track/cutsReco_AssociatorByHits/effic")
    effic_vT = fs_file.Get("DQMData/Tracking/Track/Vtxtracks_AssociatorByHits/effic")
    effic_aM = fs_file.Get("DQMData/Tracking/Track/FirstVertexTrackCollection_AssociatorByHits/effic")
    
    num_cR = fs_file.Get("DQMData/Tracking/Track/cutsReco_AssociatorByHits/num_reco_eta")
    
    print " Got all effic plots from " + fileName
    
    if WithRef:
	    	    
	scen_num = 1
	if str(scenario) == "FS08":
	    scen_num = 2
	elif str(scenario) == "FS09":
	    scen_num = 3
	elif str(scenario) == "FS10":
	    scen_num = 4
	
	
        ref_fileName = "../files/FS2012_RECO-Tracks_Draft_TV.root"
        ref_file = ROOT.TFile.Open(ref_fileName)
    
        ref_effic_cR = ref_file.Get("DQMData/Tracking/Track/cutsReco_AssociatorByHits/effic")
        ref_effic_vT = ref_file.Get("DQMData/Tracking/Track/Vtxtracks_AssociatorByHits/effic")
        ref_effic_aM = ref_file.Get("DQMData/Tracking/Track/FirstVertexTrackCollection_AssociatorByHits/effic")
    
        ref_num_cR = ref_file.Get("DQMData/Tracking/Track/cutsReco_AssociatorByHits/num_reco_eta")
    
        print " Got the reference plots from " + ref_fileName
	
	efficRatios[scen_num-1] = GetRatio(effic_cR,ref_effic_cR)
	numRatios[scen_num-1] = GetRatio(num_cR,ref_num_cR)
	
    ROOT.SetRootPalette(1)
  
    effic_canv = ROOT.TCanvas("effic_canv","Guck mal",10,10,1200,850)
    num_canv = ROOT.TCanvas("num_canv","Guck mal",10,10,1200,850)
  
    if WithRef:
    
        effic_canv.Divide(1, 2, 0.01, 0.0)
   
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetBorderMode(0)
        ROOT.gPad.SetFillColor(0)
   
        effic_canv.cd(1)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetPad(0.0, 0.25, 1.0, 1.0)
        ROOT.gPad.SetTopMargin(0.1)
        ROOT.gPad.SetLeftMargin(0.16) #0.13
        ROOT.gPad.SetRightMargin(0.04) #0.05
        ROOT.gPad.SetFillColor(0)

        ref_effic_cR.GetXaxis().SetTitle("#eta")
        ref_effic_cR.GetYaxis().SetTitle("efficiency")
        ref_effic_cR.GetXaxis().CenterTitle()
        ref_effic_cR.GetYaxis().CenterTitle()
        ref_effic_cR.GetYaxis().SetRangeUser(0.5,1.05)
        ref_effic_cR.GetXaxis().SetTitleOffset(1.1)
        ref_effic_cR.GetYaxis().SetTitleOffset(1.0)
		
        ref_effic_cR.SetMarkerStyle(20)
        ref_effic_cR.SetMarkerColor(1)
        ref_effic_cR.SetLineColor(1)
        ref_effic_cR.SetLineWidth(2)
        ref_effic_cR.SetStats(0)
  
        ref_effic_vT.SetMarkerStyle(24)
        ref_effic_vT.SetMarkerColor(1)
        ref_effic_vT.SetLineColor(1)
        ref_effic_vT.SetLineWidth(2)
        ref_effic_vT.SetStats(0)
  
        ref_effic_aM.SetMarkerStyle(25)
        ref_effic_aM.SetMarkerColor(1)
        ref_effic_aM.SetLineColor(1)
        ref_effic_aM.SetLineWidth(2)
        ref_effic_aM.SetStats(0)
	
        ref_effic_cR.Draw()
        #ref_effic_vT.Draw("same")
        #ref_effic_aM.Draw("same")
    
        effic_cR.SetMarkerStyle(20)
        effic_cR.SetMarkerColor(int(my_colors[scen_num]))
        effic_cR.SetLineColor(int(my_colors[scen_num]))
        effic_cR.SetLineWidth(2)
	
        effic_vT.SetMarkerStyle(24)
        effic_vT.SetMarkerColor(int(my_colors[scen_num]))
        effic_vT.SetLineColor(int(my_colors[scen_num]))
        effic_vT.SetLineWidth(2)
	
        effic_aM.SetMarkerStyle(25)
        effic_aM.SetMarkerColor(int(my_colors[scen_num]))
        effic_aM.SetLineColor(int(my_colors[scen_num]))
        effic_aM.SetLineWidth(2)
		
        effic_cR.Draw("same")
        #effic_vT.Draw("same")
        #effic_aM.Draw("same")	
   
        effic_canv.cd(2)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.25)
        ROOT.gPad.SetBottomMargin(0.375)
        ROOT.gPad.SetLeftMargin(0.16) #0.13
        ROOT.gPad.SetRightMargin(0.04) #0.05
        ROOT.gPad.SetFillColor(0)
    
        efficRatios[0].GetYaxis().SetTitleSize(0.1) #0.11
        efficRatios[0].GetYaxis().SetTitleOffset(0.3) #0.55
        efficRatios[0].GetYaxis().SetLabelSize(0.1)
	        
        efficRatios[0].GetXaxis().SetTitleSize(0.16)
        efficRatios[0].GetXaxis().SetLabelSize(0.16)

        efficRatios[0].GetXaxis().CenterTitle()
        efficRatios[0].GetYaxis().CenterTitle()
        efficRatios[0].GetXaxis().SetTitleOffset(1.1)
        efficRatios[0].GetYaxis().SetRangeUser(0.8,1.05)
    
        efficRatios[0].SetLineColor(int(my_colors[1]))
        efficRatios[0].SetLineWidth(2)
        efficRatios[0].SetStats(0)
        efficRatios[0].Draw()
    
        for i in range(1,scen_num):
    
            efficRatios[i].SetLineColor(int(my_colors[i+1]))
            efficRatios[i].SetLineWidth(2)
            efficRatios[i].SetStats(0)
            efficRatios[i].Draw("same")
	
	##################################
	
	num_canv.Divide(1, 2, 0.01, 0.0)
   
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetBorderMode(0)
        ROOT.gPad.SetFillColor(0)
   
        num_canv.cd(1)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetPad(0.0, 0.25, 1.0, 1.0)
        ROOT.gPad.SetTopMargin(0.1)
        ROOT.gPad.SetLeftMargin(0.16) #0.13
        ROOT.gPad.SetRightMargin(0.04) #0.05
        ROOT.gPad.SetFillColor(0)
        ROOT.gStyle.SetOptStat("i")

        ref_num_cR.GetXaxis().SetTitle("#eta")
        ref_num_cR.GetYaxis().SetTitle("number of tracks")
        ref_num_cR.GetXaxis().CenterTitle()
        ref_num_cR.GetYaxis().CenterTitle()
        ref_num_cR.GetXaxis().SetTitleOffset(1.1)
        ref_num_cR.GetYaxis().SetTitleOffset(1.0)
		
        ref_num_cR.SetMarkerStyle(20)
        ref_num_cR.SetMarkerColor(1)
        ref_num_cR.SetLineColor(1)
        ref_num_cR.SetLineWidth(2)
  
        ref_num_cR.Draw()
        moveStatBox1(num_canv,ref_num_cR,int(my_colors[0]),0)
    
        num_cR.SetMarkerStyle(20)
        num_cR.SetMarkerColor(int(my_colors[scen_num]))
        num_cR.SetLineColor(int(my_colors[scen_num]))
        num_cR.SetLineWidth(2)
		
        num_cR.Draw("sames")	
        moveStatBox1(num_canv,num_cR,int(my_colors[scen_num]),1)
   
        num_canv.cd(2)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.25)
        ROOT.gPad.SetBottomMargin(0.375)
        ROOT.gPad.SetLeftMargin(0.16) #0.13
        ROOT.gPad.SetRightMargin(0.04) #0.05
        ROOT.gPad.SetFillColor(0)
    
        numRatios[0].GetYaxis().SetTitleSize(0.1) #0.11
        numRatios[0].GetYaxis().SetTitleOffset(0.3) #0.55
        numRatios[0].GetYaxis().SetLabelSize(0.1)
	        
        numRatios[0].GetXaxis().SetTitleSize(0.16)
        numRatios[0].GetXaxis().SetLabelSize(0.16)

        numRatios[0].GetXaxis().CenterTitle()
        numRatios[0].GetYaxis().CenterTitle()
        numRatios[0].GetXaxis().SetTitleOffset(1.1)
        numRatios[0].GetYaxis().SetRangeUser(0.8,1.05)
    
        numRatios[0].SetLineColor(int(my_colors[1]))
        numRatios[0].SetLineWidth(2)
        numRatios[0].SetStats(0)
        numRatios[0].Draw()
    
        for i in range(1,scen_num):
    
            numRatios[i].SetLineColor(int(my_colors[i+1]))
            numRatios[i].SetLineWidth(2)
            numRatios[i].Draw("same")
	    		
    else:  
    
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetBorderMode(0)
        ROOT.gPad.SetFillColor(0)
        effic_canv.Divide(1,1)
   
        effic_canv.cd(1)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetRightMargin(0.12)
        ROOT.gPad.SetLeftMargin(0.08)
        ROOT.gPad.SetFillColor(0)

        effic_cR.GetXaxis().SetTitle("#eta")
        effic_cR.GetYaxis().SetTitle("efficiency")
        effic_cR.GetXaxis().CenterTitle()
        effic_cR.GetYaxis().CenterTitle()
        effic_cR.GetYaxis().SetRangeUser(0.6,1.05)
        effic_cR.GetXaxis().SetTitleOffset(1.1)
        effic_cR.GetYaxis().SetTitleOffset(1.0)
    
        effic_cR.SetMarkerColor(1)
        effic_cR.SetLineColor(1)
        effic_cR.SetLineWidth(2)
        effic_cR.SetStats(0)
   
        effic_vT.SetMarkerColor(1)
        effic_vT.SetLineColor(1)
        effic_vT.SetLineWidth(2)
        effic_vT.SetStats(0)
  
        effic_aM.SetMarkerColor(1)
        effic_aM.SetLineColor(1)
        effic_aM.SetLineWidth(2)
        effic_aM.SetStats(0)
	
        effic_cR.SetMarkerStyle(20)
        effic_vT.SetMarkerStyle(24)
        effic_aM.SetMarkerStyle(25)
	
        effic_cR.Draw()
        #effic_vT.Draw("same")
        #effic_aM.Draw("same")
	
	############################
    
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetBorderMode(0)
        ROOT.gPad.SetFillColor(0)
        num_canv.Divide(1,1)
   
        num_canv.cd(1)
        ROOT.gPad.SetGrid()
        ROOT.gPad.SetRightMargin(0.12)
        ROOT.gPad.SetLeftMargin(0.08)
        ROOT.gPad.SetFillColor(0)
        ROOT.gStyle.SetOptStat("i")

        num_cR.GetXaxis().SetTitle("#eta")
        num_cR.GetYaxis().SetTitle("number of tracks")
        num_cR.GetXaxis().CenterTitle()
        num_cR.GetYaxis().CenterTitle()
        num_cR.GetXaxis().SetTitleOffset(1.1)
        num_cR.GetYaxis().SetTitleOffset(1.0)
    
        num_cR.SetMarkerColor(1)
        num_cR.SetLineColor(1)
        num_cR.SetLineWidth(2)
	
        num_cR.SetMarkerStyle(20)
	
        num_cR.Draw()
	moveStatBox1(num_canv,num_cR,int(my_colors[0]),0)
      
    effic_canv.SaveAs("../pictures/RECO-Tracks_" + str(scenario) + "_TrackValidation_effic.png")   
    num_canv.SaveAs("../pictures/RECO-Tracks_" + str(scenario) + "_TrackValidation_num.png")  
               
def CreateMassFit(scenarios,V0):
	
    num_of_scen = len(scenarios)
    
    finscen = scenarios[num_of_scen-1]
    
    print "Will create plot for scenario " + finscen
    
    files =[[] for i in range(num_of_scen)]
    distrs =[[] for i in range(num_of_scen)]
    fits =[[] for i in range(num_of_scen)]
    
    if str(V0)=="Ks":
        gausFit = ROOT.TF1("gausFit","TMath::Gaus(x,[0],[1])*[2]+[3]",400.,600.)
    elif str(V0)=="Lambda":
        gausFit = ROOT.TF1("gausFit","TMath::Gaus(x,[0],[1])*[2]+[3]+[4]*x",1.05,1.2)
    
    for i,scenario in enumerate(scenarios):
	    
	if i==0 or i==(num_of_scen-1):
    
            files[i] = ROOT.TFile.Open("../files/FS2012_RECO-Tracks_" + str(scenario) + ".root") 
            distrs[i] = files[i].Get("Ksmass/RECO-Tracks_" + str(scenario) + "_" + str(V0) + "MassDistribution")
    
            print "\n Opened ../files/FS2012_RECO-Tracks_" + str(scenario) + ".root"
	    
            if str(V0)=="Ks":
                gausFit.SetParameters(distrs[i].GetMean(),distrs[i].GetRMSError(),3.,250.)
            elif str(V0)=="Lambda":
                gausFit.SetParameters(distrs[i].GetMean(),distrs[i].GetRMSError(),3.5,10.,0.5)           
	    
            if str(V0)=="Ks":
    	        distrs[i].Fit("gausFit","0","I",450.,550.)
            elif str(V0)=="Lambda":
    	        distrs[i].Fit("gausFit","0","I",1.09,1.14)
    	    fits[i] = distrs[i].GetFunction("gausFit")
	    

    dummy = ROOT.TGraphErrors(1)
    
    x = fits[num_of_scen-1].GetParameter(0)
    x_err = abs(fits[num_of_scen-1].GetParameter(1))
    
    dummy.SetPoint(1,x,num_of_scen)
    dummy.SetPointError(1,x_err,0.2)
     
    if str(V0)=="Ks":
        dummy.GetXaxis().SetTitle("K_{s} mass / MeV")
    elif str(V0)=="Lambda":
        dummy.GetXaxis().SetTitle("#Lambda mass / GeV")
	
    dummy.GetYaxis().SetTitle("Failure scenario")
	    
    errGraphs[num_of_scen-1] = dummy
  
    distr_canv = ROOT.TCanvas("distr_canv","Guck mal",10,10,1200,850)
    ROOT.SetRootPalette(1)
    ROOT.gStyle.SetOptFit(0)
    
    distr_canv.Divide(1, 2, 0.01, 0.0)
  
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetBorderMode(0)
    ROOT.gPad.SetFillColor(0)
   
    distr_canv.cd(1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetPad(0.0, 0.25, 1.0, 1.0)
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.16) #0.13
    ROOT.gPad.SetRightMargin(0.04) #0.05
    ROOT.gPad.SetFillColor(0)
    ROOT.gStyle.SetOptStat(0)

    distrs[0].GetXaxis().CenterTitle()
    distrs[0].GetYaxis().CenterTitle()
    distrs[0].GetXaxis().SetTitleOffset(1.1)
    distrs[0].GetYaxis().SetTitleOffset(1.0)
    
    distrs[0].SetLineColor(int(my_colors[0]))
    distrs[0].SetLineWidth(2)
    distrs[0].Draw()
    
    statboxR = ROOT.TPaveStats(0.2,0.7,0.45,0.85,"brNDC")	

    statboxR.SetLineWidth(2)
    statboxR.SetLineColor(int(my_colors[0]))
    statboxR.SetFillColor(0)
    statboxR.SetBorderSize(1)
    statboxR.SetOptFit(1)
    
    statboxR.AddText("Mean  = %3.4f" % fits[0].GetParameter(0))
    statboxR.AddText("Error = %3.4f" % abs(fits[0].GetParameter(1)))
    statboxR.AddText("Entries = %5.f" % distrs[0].GetEntries())
    
    statboxR.Draw()
   
    fits[0].SetLineColor(int(my_colors[0]))
    fits[0].SetLineWidth(2)
    fits[0].Draw("same")
    
    if num_of_scen > 1:
    
        distrs[num_of_scen-1].SetLineColor(int(my_colors[num_of_scen-1]))
        distrs[num_of_scen-1].SetLineWidth(2)
        distrs[num_of_scen-1].Draw("same")
    
        statboxF = ROOT.TPaveStats(0.7,0.7,0.95,0.85,"brNDC")	

        statboxF.SetLineWidth(2)
        statboxF.SetLineColor(int(my_colors[num_of_scen-1]))
        statboxF.SetFillColor(0)
        statboxF.SetBorderSize(1)
        statboxF.SetOptFit(1)
    
        statboxF.AddText("Mean  = %.4f" % fits[num_of_scen-1].GetParameter(0))
        statboxF.AddText("Error = %.4f" % abs(fits[num_of_scen-1].GetParameter(1)))
        statboxF.AddText("Entries = %5.f" % distrs[num_of_scen-1].GetEntries())
      
        statboxF.Draw()
      
        fits[num_of_scen-1].SetLineColor(int(my_colors[num_of_scen-1]))
        fits[num_of_scen-1].SetLineWidth(2)
        fits[num_of_scen-1].Draw("same")
	
     
    if str(V0)=="Ks":
        dummyH = ROOT.TH2F("dummyH","; K_{s} mass / MeV; FS",100,400.,600.,5,0.5,5.5)
    elif str(V0)=="Lambda":
        dummyH = ROOT.TH2F("dummyH","; #Lambda mass / GeV; FS",100,1.05,1.2,5,0.5,5.5)
	
   
    distr_canv.cd(2)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.25)
    ROOT.gPad.SetBottomMargin(0.375)
    ROOT.gPad.SetLeftMargin(0.16) #0.13
    ROOT.gPad.SetRightMargin(0.04) #0.05
    ROOT.gPad.SetFillColor(0)
    
    dummyH.GetYaxis().SetTitleSize(0.16) #0.11
    dummyH.GetYaxis().SetTitleOffset(0.3) #0.55
    dummyH.GetYaxis().SetLabelSize(0.16)		
        
    dummyH.GetXaxis().SetTitleSize(0.16)
    dummyH.GetXaxis().SetLabelSize(0.16)

    dummyH.GetXaxis().CenterTitle()
    dummyH.GetYaxis().CenterTitle()
    dummyH.GetXaxis().SetTitleOffset(1.1)
    dummyH.GetYaxis().SetRangeUser(0.,6.)
    
    dummyH.GetYaxis().SetBinLabel(1,"Draft")
    dummyH.GetYaxis().SetBinLabel(2,"FS01")
    dummyH.GetYaxis().SetBinLabel(3,"FS08")
    dummyH.GetYaxis().SetBinLabel(4,"FS09")
    dummyH.GetYaxis().SetBinLabel(5,"FS10")	
    
    dummyH.GetYaxis().SetRangeUser(0.,6.)
    
    dummyH.Draw()
    
    for i in range(0,num_of_scen):
    
        errGraphs[i].SetLineColor(int(my_colors[i]))
        errGraphs[i].SetLineWidth(2)
        errGraphs[i].Draw("P")
    
     
    distr_canv.SaveAs("../pictures/RECO-Tracks_" + str(scenario) + "_" + str(V0) + "MassFit.png")