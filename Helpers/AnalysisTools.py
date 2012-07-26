#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms
import subprocess, ROOT, random, math, os

ROOT.gROOT.Macro( os.path.expanduser( '~/Phd-Study/CMSSW/Helpers/SetPrettyStyle.C' ) )

possibleDatasets = ["QCD","DYToMuMu","TT"]

spectra = ["15to30","30to50","50to80","80to120","120to170","170to300","300to470","470to600","600to800","800to1000","15to3000"]

weights = {"15to3000":1.,"15to30":8.159e08,"30to50":5.312e07,"50to80":6.359e06,"80to120":7.843e05,"120to170":1.151e05,"170to300":2.426e04,"300to470":1.168e04,"470to600":7.022e01,"600to800":1.555e01,"800to1000":1.844e00}

my_colors = [1,600,634,887,857,418,432,807]
#my_colorsP = [634,887,1]
my_colorsP = [634,887,1,857]

colorsMap = {1:"Black",600:"Blue",634:"Red",887:"Violet",857:"Azure",418:"Green",432:"Cyan",807:"Orange"}


pt_ranges = [30.,100.,300.,1000.]
pt_names = ["cR","30.to100.","100.to300.","300.to1000."]
pt_namesP = ["30to100","100to300","300to1000","cR"]

eta_ranges = [-5.,-2.5,-1.3,1.3]
eta_names = ["cR","2.5to5","1.3to2.5","0.to1.3"]
eta_namesP = ["2.5to5","1.3to2.5","0.to1.3","cR"]

npu_ranges = [0,15,30,50]
npu_names = ["cR","0to15","15to30","30to50"]
npu_namesP = ["0to15","15to30","30to50","cR"]


etaNo1_shadow = ROOT.TPave(-4.9,0.42,-2.6,1.03)
etaNo1_shadow.SetFillColor(921)
etaNo1_shadow.SetBorderSize(0)

etaNo2_shadow = ROOT.TPave(2.6,0.42,4.9,1.03)
etaNo2_shadow.SetFillColor(921)
etaNo2_shadow.SetBorderSize(0)

etaEc1_shadow = ROOT.TPave(-2.6,0.42,-1.3,1.03)
etaEc1_shadow.SetFillColor(920)
etaEc1_shadow.SetBorderSize(0)

etaEc2_shadow = ROOT.TPave(1.3,0.42,2.6,1.03)
etaEc2_shadow.SetFillColor(920)
etaEc2_shadow.SetBorderSize(0)


def GetBinNumber(axis,value):
  
    axis_bin_number = axis.GetNbins()
    axis_min = axis.GetBinLowEdge(1)
    axis_max = axis.GetBinUpEdge(axis_bin_number)
    
    step_size = (axis_max - axis_min)*1./axis_bin_number
    
    bin_num = int((value-axis_min)*1./step_size)
    
    if bin_num>0:
        return bin_num
    else:
        return 0


def GetMedian(histo):
  
    bin_number = histo.GetNbinsX()
    
    mean = histo.GetEntries()*1./2.
    
    entries_help = 0
    idx = 0
    
    while entries_help<mean:
        entries_help+= histo.GetBinContent(idx)
	idx+= 1
	
    idx-= 1
    
    return histo.GetBinCenter(idx)


def moveStatBox(canvas, histo, color, step, dim=1, rmsred=0.):

	canvas.Update()
	
	statBox = histo.GetListOfFunctions().FindObject("stats")
	   	
	if dim==2:    
            histo.GetListOfFunctions().Remove(statBox)
            histo.SetStats(0)
	    statBox.Clear()
	    statBox.AddText("Correlation  = %1.4f" % histo.GetCorrelationFactor())
	
	    statBox.SetX1NDC(0.4)
	    statBox.SetX2NDC(0.8)
	    statBox.SetY1NDC(0.7)
	    statBox.SetY2NDC(0.83)
	else:    
            histo.GetListOfFunctions().Remove(statBox)
            histo.SetStats(0)
	    statBox.Clear()
	    statBox.AddText("Mean  = %1.4f" % histo.GetMean())
	    statBox.AddText("RMS  = %1.4f" % histo.GetRMS())
	    if not float(rmsred)==0.:
		result = float(histo.GetRMS())/float(rmsred)
		print str(histo.GetRMS()) + " / " + str(rmsred) + " = " +str(result)
		statBox.AddText("RMS*  = %1.4f" % result)
	    		
	    statBoxSpacing = 0.01
	    size_y = 0.13	
	
	    statBox.SetX1NDC(0.65)
	    statBox.SetX2NDC(0.95)
	    statBox.SetY2NDC(0.83 - step*(size_y+statBoxSpacing))
	    statBox.SetY1NDC(0.83 - (step+1)*size_y-step*statBoxSpacing)
	
	statBox.SetLineWidth(2)	    
	statBox.SetLineColor(color)
        statBox.SetBorderSize(1)
	
	statBox.Draw()

def GetListOfFiles(path):
	
	ListOfNames = []
	
	p = subprocess.Popen(["ls -l "+path], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	
	for line in p.stdout.readlines():
                filename = line.strip().split()
		if len(filename)>8:
                	#print filename[8]
			name=filename[8]
		        ListOfNames.append(path+name)
			
		
	print "" 				
	print  str(len(ListOfNames)) + " files were found"
	
	return ListOfNames

def GetListOfFilesWithEnding(path,ending):
	
	ListOfNames = []
	
	p = subprocess.Popen(["ls -l "+path+" | grep "+ending], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	
	for line in p.stdout.readlines():
                filename = line.strip().split()
		if len(filename)>8:
                	#print filename[8]
			name=filename[8]
		        ListOfNames.append(path+name)
			
		
	print "" 				
	print  str(len(ListOfNames)) + " files were found"
	
	return ListOfNames


def GetListOfSubdirectories(filename,dirname,tagname):
		
	print ""
	
	ListOfDirecs = []
	
	f = ROOT.TFile.Open(filename)
	
	if not f.cd(dirname):
	    print dirname + " does not exist in " + filename
	    return ListOfDirecs

	
	f.cd(dirname)
	
	HelpList = ROOT.gDirectory.GetListOfKeys()
	
	for i in range(len(HelpList)):
	    if tagname in HelpList[i].GetName():
	        ListOfDirecs.append(dirname + HelpList[i].GetName() + "/")
		
        if len(ListOfDirecs) == 0:
	    print "No Subdirectories with " + tagname + " were found"
	    print "Returning empty subdirectory list"
        else:
            print str(len(ListOfDirecs)) + " subdirectories were found:"
	    for i in range(len(ListOfDirecs)):
	        print " " + ListOfDirecs[i]
		
	print ""
	
	return ListOfDirecs


def GetListOfKeyWords(filename,dirname,tagname):
		
	print ""
	
	ListOfKeys = []
	
	f = ROOT.TFile.Open(filename)
	
	if not f.cd(dirname):
	    print dirname + " does not exist in " + filename
	    return ListOfKeys

	
	f.cd(dirname)
	
	HelpList = ROOT.gDirectory.GetListOfKeys()
	
	for i in range(len(HelpList)):
	    if tagname in HelpList[i].GetName():
	        ListOfKeys.append(HelpList[i].GetName())
		
        if len(ListOfKeys) == 0:
	    print "No keywords with " + tagname + " were found"
	    print "Returning empty subdirectory list"
        else:
            print str(len(ListOfKeys)) + " keywords were found:"
	    for i in range(len(ListOfKeys)):
	        print " " + ListOfKeys[i]
		
	print ""
	
	return ListOfKeys
	        
		
def GetCombinedTH1FsWeight(histos,names):
	
	print "Combining TH1Fs with weight..."

        histo_output = ()
	
	weight = 1.
    
        for ass_ite in range(len(histos[0])):
    
            op_histo_ass_axis = ()
	    
            for x_axis_ite in range(len(histos[0][0])):
  
  	        x_bin_number = histos[0][ass_ite][x_axis_ite].GetNbinsX()
  	        x_min = histos[0][ass_ite][x_axis_ite].GetXaxis().GetBinLowEdge(1)
  	        x_max = histos[0][ass_ite][x_axis_ite].GetXaxis().GetBinUpEdge(x_bin_number)		
		
	        op_name = "combined_" + histos[0][ass_ite][x_axis_ite].GetName()
		op_title =  histos[0][ass_ite][x_axis_ite].GetTitle()
  	   
	        op_histo = ROOT.TH1F(op_name,op_title,x_bin_number,x_min,x_max)

  	        op_histo.GetXaxis().SetTitle(histos[0][ass_ite][x_axis_ite].GetXaxis().GetTitle())
  	        op_histo.GetYaxis().SetTitle(histos[0][ass_ite][x_axis_ite].GetYaxis().GetTitle())
		
		entries_help = 0
	    
                for file_ite in range(len(histos)):
	    
	            name = names[file_ite]
	
                    for spectrum in spectra:
		        if spectrum in name:
	                    weight = weights[spectrum]
			    break
	    
                    for bin_ite in range(x_bin_number):
			entries_help+= histos[file_ite][ass_ite][x_axis_ite].GetBinContent(bin_ite+1)
	                addOn = histos[file_ite][ass_ite][x_axis_ite].GetBinContent(bin_ite+1)*weight
		        new = op_histo.GetBinContent(bin_ite+1) + addOn
		        op_histo.SetBinContent(bin_ite+1,new)
			
	        op_histo.SetEntries(entries_help)			
	    
  		op_histo_ass_axis+=op_histo,
		 
            histo_output+=op_histo_ass_axis,
		    
		    		    
	return  histo_output 
	        
		
def GetCombinedTH1FsWeight5D(histos,names):
	
	print "Combining TH1Fs with weight..."
        print ""

        histo_output = ()
	
	weight = 1.
	
	dim1 = len(histos[0])
	dim2 = len(histos[0][0])
	dim3 = len(histos[0][0][0])
	dim4 = len(histos[0][0][0][0])
	dim5 = len(histos[0][0][0][0][0])
	
	overall = (dim1*dim2*dim3*dim4*dim5)*1./20.
	step = overall
	counter = 0
    
        for ite1 in range(dim1):
    
            op_histo_1 = ()
	    
            for ite2 in range(dim2):
    
                op_histo_2 = ()
	    
                for ite3 in range(dim3):
    
                    op_histo_3 = ()
	    
                    for ite4 in range(dim4):
    
                        op_histo_4 = ()
	    
                        for ite5 in range(dim5):
  
  	                    x_bin_number = histos[0][ite1][ite2][ite3][ite4][ite5].GetNbinsX()
  	                    x_min = histos[0][ite1][ite2][ite3][ite4][ite5].GetXaxis().GetBinLowEdge(1)
  	                    x_max = histos[0][ite1][ite2][ite3][ite4][ite5].GetXaxis().GetBinUpEdge(x_bin_number)
		
	                    op_name = "combined_" + histos[0][ite1][ite2][ite3][ite4][ite5].GetName()
  	   
	                    op_histo = ROOT.TH1F(op_name,"",x_bin_number,x_min,x_max)
		
	                    entries_help = 0
	    
                            for file_ite in range(len(histos)):
	    
	                        name = names[file_ite]
        	
                                for spectrum in spectra:
		                    if spectrum in name:
	                                weight = weights[spectrum]
	    
                                for x_bin_ite in range(x_bin_number):
			            entries_help+= histos[file_ite][ite1][ite2][ite3][ite4][ite5].GetBinContent(x_bin_ite+1)
   	                            addOn = histos[file_ite][ite1][ite2][ite3][ite4][ite5].GetBinContent(x_bin_ite+1)*weight
		                    new = op_histo.GetBinContent(x_bin_ite+1) + addOn
		                    op_histo.SetBinContent(x_bin_ite+1,new)		
			
	                    op_histo.SetEntries(entries_help)			
	    				    
  		            op_histo_4+=op_histo,
			    
			    counter+= 1
  				
			    if(counter>overall):
			        print " = ",
				overall+= step	
	    				    
  		        op_histo_3+=op_histo_4,			
	    				    
  		    op_histo_2+=op_histo_3,		
	    				    
  		op_histo_1+=op_histo_2,	
		 
            histo_output+=op_histo_1,
		    
 	print ""
		    		    
	return  histo_output
	        
		
def GetCombinedTH2FsWeight(histos,names):
	
	print "Combining TH2Fs with weight..."

        histo_output = ()
	
	weight = 1.
    
        for ass_ite in range(len(histos[0])):
  
  	    x_bin_number = histos[0][ass_ite].GetNbinsX()
  	    x_min = histos[0][ass_ite].GetXaxis().GetBinLowEdge(1)
  	    x_max = histos[0][ass_ite].GetXaxis().GetBinUpEdge(x_bin_number)
  
  	    y_bin_number = histos[0][ass_ite].GetNbinsY()
  	    y_min = histos[0][ass_ite].GetYaxis().GetBinLowEdge(1)
  	    y_max = histos[0][ass_ite].GetYaxis().GetBinUpEdge(y_bin_number)
		
	    op_name = "combined_" + histos[0][ass_ite].GetName()
            op_title =  histos[0][ass_ite].GetTitle()
  	   
	    op_histo = ROOT.TH2F(op_name,op_title,x_bin_number,x_min,x_max,y_bin_number,y_min,y_max)

  	    op_histo.GetXaxis().SetTitle(histos[0][ass_ite].GetXaxis().GetTitle())
  	    op_histo.GetYaxis().SetTitle(histos[0][ass_ite].GetYaxis().GetTitle())
  	    op_histo.GetZaxis().SetTitle(histos[0][ass_ite].GetZaxis().GetTitle())
		
	    entries_help = 0
	    
            for file_ite in range(len(histos)):
	    
	        name = names[file_ite]
	
                for spectrum in spectra:
		    if spectrum in name:
	                weight = weights[spectrum]
			break
	    
                for x_bin_ite in range(x_bin_number):
                    for y_bin_ite in range(y_bin_number):
			entries_help+= histos[file_ite][ass_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1)
   	                addOn = histos[file_ite][ass_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1)*weight
		        new = op_histo.GetBinContent(x_bin_ite+1,y_bin_ite+1) + addOn
		        op_histo.SetBinContent(x_bin_ite+1,y_bin_ite+1,new)
			
			
	    op_histo.SetEntries(entries_help)				    
		 
            histo_output+=op_histo,
		    
		    		    
	return  histo_output
	        
		
def GetCombinedTH2FsWeight2(histos,names):
	
	print "Combining TH2Fs with weight..."

        histo_output = ()
	
	weight = 1.
    
        for ass_ite in range(len(histos[0])):
    
            op_histo_dir_ass = ()
	    
            for dir_ite in range(len(histos[0][0])):
  
  	        x_bin_number = histos[0][ass_ite][dir_ite].GetNbinsX()
  	        x_min = histos[0][ass_ite][dir_ite].GetXaxis().GetBinLowEdge(1)
  	        x_max = histos[0][ass_ite][dir_ite].GetXaxis().GetBinUpEdge(x_bin_number)
    
  	        y_bin_number = histos[0][ass_ite][dir_ite].GetNbinsY()
  	        y_min = histos[0][ass_ite][dir_ite].GetYaxis().GetBinLowEdge(1)
  	        y_max = histos[0][ass_ite][dir_ite].GetYaxis().GetBinUpEdge(y_bin_number)
		
	        op_name = "combined_" + histos[0][ass_ite][dir_ite].GetName()
                op_title =  histos[0][ass_ite][dir_ite].GetName()
  	   
	        op_histo = ROOT.TH2F(op_name,op_title,x_bin_number,x_min,x_max,y_bin_number,y_min,y_max)

  	        op_histo.GetXaxis().SetTitle(histos[0][ass_ite][dir_ite].GetXaxis().GetTitle())
  	        op_histo.GetYaxis().SetTitle(histos[0][ass_ite][dir_ite].GetYaxis().GetTitle())
  	        op_histo.GetZaxis().SetTitle(histos[0][ass_ite][dir_ite].GetZaxis().GetTitle())
		
	        entries_help = 0
	    
                for file_ite in range(len(histos)):
	    
	            name = names[file_ite]
	
                    for spectrum in spectra:
		        if spectrum in name:
	                    weight = weights[spectrum]
	    
                    for x_bin_ite in range(x_bin_number):
                        for y_bin_ite in range(y_bin_number):
			    entries_help+= histos[file_ite][ass_ite][dir_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1)
   	                    addOn = histos[file_ite][ass_ite][dir_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1)*weight
		            new = op_histo.GetBinContent(x_bin_ite+1,y_bin_ite+1) + addOn
		            op_histo.SetBinContent(x_bin_ite+1,y_bin_ite+1,new)		
			
	        op_histo.SetEntries(entries_help)			
	    				    
  		op_histo_dir_ass+=op_histo,
		 
            histo_output+=op_histo_dir_ass,
		    
		    		    
	return  histo_output
	        
		
def GetCombinedTH2FsWeight4D(histos,names):
	
	print "Combining TH2Fs with weight..."
        print ""

        histo_output = ()
	
	weight = 1.
	
	dim1 = len(histos[0])
	dim2 = len(histos[0][0])
	dim3 = len(histos[0][0][0])
	dim4 = len(histos[0][0][0][0])
	
	overall = (dim1*dim2*dim3*dim4)*1./20.
	step = overall
	counter = 0
    
        for ite1 in range(len(histos[0])):
    
            op_histo_1 = ()
	    
            for ite2 in range(len(histos[0][0])):
    
                op_histo_2 = ()
	    
                for ite3 in range(len(histos[0][0][0])):
    
                    op_histo_3 = ()
	    
                    for ite4 in range(len(histos[0][0][0][0])):
  
  	                x_bin_number = histos[0][ite1][ite2][ite3][ite4].GetNbinsX()
  	                x_min = histos[0][ite1][ite2][ite3][ite4].GetXaxis().GetBinLowEdge(1)
  	                x_max = histos[0][ite1][ite2][ite3][ite4].GetXaxis().GetBinUpEdge(x_bin_number)
    
  	                y_bin_number = histos[0][ite1][ite2][ite3][ite4].GetNbinsY()
  	                y_min = histos[0][ite1][ite2][ite3][ite4].GetYaxis().GetBinLowEdge(1)
  	                y_max = histos[0][ite1][ite2][ite3][ite4].GetYaxis().GetBinUpEdge(y_bin_number)
		
	                op_name = "combined_" + histos[0][ite1][ite2][ite3][ite4].GetName()
                        op_title =  histos[0][ite1][ite2][ite3][ite4].GetName()
  	   
	                op_histo = ROOT.TH2F(op_name,op_title,x_bin_number,x_min,x_max,y_bin_number,y_min,y_max)

  	                op_histo.GetXaxis().SetTitle(histos[0][ite1][ite2][ite3][ite4].GetXaxis().GetTitle())
  	                op_histo.GetYaxis().SetTitle(histos[0][ite1][ite2][ite3][ite4].GetYaxis().GetTitle())
  	                op_histo.GetZaxis().SetTitle(histos[0][ite1][ite2][ite3][ite4].GetZaxis().GetTitle())
		
	                entries_help = 0
	    
                        for file_ite in range(len(histos)):
	    
	                    name = names[file_ite]
        	
                            for spectrum in spectra:
		                if spectrum in name:
	                            weight = weights[spectrum]
	    
                            for x_bin_ite in range(x_bin_number):
                                for y_bin_ite in range(y_bin_number):
			          entries_help+= histos[file_ite][ite1][ite2][ite3][ite4].GetBinContent(x_bin_ite+1,y_bin_ite+1)
   	                          addOn = histos[file_ite][ite1][ite2][ite3][ite4].GetBinContent(x_bin_ite+1,y_bin_ite+1)*weight
		                  new = op_histo.GetBinContent(x_bin_ite+1,y_bin_ite+1) + addOn
		                  op_histo.SetBinContent(x_bin_ite+1,y_bin_ite+1,new)		
			
	                op_histo.SetEntries(entries_help)			
	    				    
  		        op_histo_3+=op_histo,
			    
			counter+= 1
  			
			if(counter>overall):
			    print " = ",
			    overall+= step			
	    				    
  		    op_histo_2+=op_histo_3,	
	    				    
  		op_histo_1+=op_histo_2,	
		 
            histo_output+=op_histo_1,
		    
	print ""
				    
	return  histo_output
	        
		
def GetCombinedTH3FsWeight(histos,names):
	
	print "Combining TH3Fs with weight..."
	
	#if len(histos)==1:
	    #return histos[0]

        histo_output = ()
	
	weight = 1.
    
        for ass_ite in range(len(histos[0])):
  
  	    x_bin_number = histos[0][ass_ite].GetNbinsX()
  	    x_min = histos[0][ass_ite].GetXaxis().GetBinLowEdge(1)
  	    x_max = histos[0][ass_ite].GetXaxis().GetBinUpEdge(x_bin_number)
  
  	    y_bin_number = histos[0][ass_ite].GetNbinsY()
  	    y_min = histos[0][ass_ite].GetYaxis().GetBinLowEdge(1)
  	    y_max = histos[0][ass_ite].GetYaxis().GetBinUpEdge(y_bin_number)
  
  	    z_bin_number = histos[0][ass_ite].GetNbinsZ()
  	    z_min = histos[0][ass_ite].GetZaxis().GetBinLowEdge(1)
  	    z_max = histos[0][ass_ite].GetZaxis().GetBinUpEdge(z_bin_number)
		
	    help = random.random()
            op_name = "combined" + str(help)
	    op_title =  histos[0][ass_ite].GetTitle()
  	   
	    op_histo = ROOT.TH3F(op_name,op_title, x_bin_number,x_min,x_max,y_bin_number,y_min,y_max,z_bin_number,z_min,z_max)

  	    op_histo.GetXaxis().SetTitle(histos[0][ass_ite].GetXaxis().GetTitle())
  	    op_histo.GetYaxis().SetTitle(histos[0][ass_ite].GetYaxis().GetTitle())
  	    op_histo.GetZaxis().SetTitle(histos[0][ass_ite].GetZaxis().GetTitle())
		
            entries_help = 0
	    
            for file_ite in range(len(histos)):
	    
	        name = names[file_ite]
	
                for spectrum in spectra:
		    if spectrum in name:
	                weight = weights[spectrum]
	    
	    
                for x_bin_ite in range(x_bin_number):
                    for y_bin_ite in range(y_bin_number):
                        for z_bin_ite in range(z_bin_number):
			    entries_help+= histos[file_ite][ass_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1)
	                    addOn = histos[file_ite][ass_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1)*weight
		            new = op_histo.GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1) + addOn
		            op_histo.SetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1,new)
			    	
	    op_histo.SetEntries(entries_help)
		 
            histo_output+=op_histo,
		    
		    		    
	return  histo_output 
	        
		
def GetCombinedTH3FsWeight2(histos,names):
	
	print "Combining TH3Fs with weight..."
	
	if len(histos)==1:
	    return histos[0]

        histo_output = ()
	
	weight = 1.
    
        for ass_ite in range(len(histos[0])):
    
            op_histo_dir_ass = ()
	    
            for dir_ite in range(len(histos[0][0])):
  
  	        x_bin_number = histos[0][ass_ite][dir_ite].GetNbinsX()
  	        x_min = histos[0][ass_ite][dir_ite].GetXaxis().GetBinLowEdge(1)
  	        x_max = histos[0][ass_ite][dir_ite].GetXaxis().GetBinUpEdge(x_bin_number)
  
  	        y_bin_number = histos[0][ass_ite][dir_ite].GetNbinsY()
  	        y_min = histos[0][ass_ite][dir_ite].GetYaxis().GetBinLowEdge(1)
  	        y_max = histos[0][ass_ite][dir_ite].GetYaxis().GetBinUpEdge(y_bin_number)
  
  	        z_bin_number = histos[0][ass_ite][dir_ite].GetNbinsZ()
  	        z_min = histos[0][ass_ite][dir_ite].GetZaxis().GetBinLowEdge(1)
  	        z_max = histos[0][ass_ite][dir_ite].GetZaxis().GetBinUpEdge(z_bin_number)
		
	        help = random.random()
		op_name = "combined" + str(help)
		op_title =  histos[0][dir_ite][dir_ite].GetTitle()
  	   
	        op_histo = ROOT.TH3F(op_name,op_title, x_bin_number,x_min,x_max,y_bin_number,y_min,y_max,z_bin_number,z_min,z_max)

  	        op_histo.GetXaxis().SetTitle(histos[0][ass_ite][dir_ite].GetXaxis().GetTitle())
  	        op_histo.GetYaxis().SetTitle(histos[0][ass_ite][dir_ite].GetYaxis().GetTitle())
  	        op_histo.GetZaxis().SetTitle(histos[0][ass_ite][dir_ite].GetZaxis().GetTitle())
		
		entries_help = 0
	    
                for file_ite in range(len(histos)):
	    
	            name = names[file_ite]
	
                    for spectrum in spectra:
		        if spectrum in name:
	                    weight = weights[spectrum]
	    
	    
                    for x_bin_ite in range(x_bin_number):
                        for y_bin_ite in range(y_bin_number):
                            for z_bin_ite in range(z_bin_number):
			        entries_help+= histos[file_ite][ass_ite][dir_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1)
	                        addOn = histos[file_ite][ass_ite][dir_ite].GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1)*weight
		                new = op_histo.GetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1) + addOn
		                op_histo.SetBinContent(x_bin_ite+1,y_bin_ite+1,z_bin_ite+1,new)
			
			
	        op_histo.SetEntries(entries_help)			
	    				    
  		op_histo_dir_ass+=op_histo,
		 
            histo_output+=op_histo_dir_ass,
		    
		    		    
	return  histo_output  

	        
		
def GetMean(histo):	
  
  	x_bin_number = histo.GetNbinsX()
  	x_min = histo.GetXaxis().GetBinLowEdge(1)
  	x_max = histo.GetXaxis().GetBinUpEdge(x_bin_number)
		
	help = random.random()
	op_name = str(histo.GetName()) + "_mean_" + str(help)
	op_title = str(histo.GetTitle()) + " - MEAN"
  	   
	op_histo = ROOT.TH1F(op_name,op_title,x_bin_number,x_min,x_max)
	
	for x_bin_ite in range(1,x_bin_number+1):	
	     histo.GetXaxis().SetRange(x_bin_ite,x_bin_ite)
	     op_histo.SetBinContent(x_bin_ite,histo.Project3D("zo").GetMean())
	     op_histo.SetBinError(x_bin_ite,0.0001)
	     	     
	return op_histo 

	        
		
def GetMean2(histo):	
  
  	x_bin_number = histo.GetNbinsX()
  	x_min = histo.GetXaxis().GetBinLowEdge(1)
  	x_max = histo.GetXaxis().GetBinUpEdge(x_bin_number)
		
	op_name = str(histo.GetName()) + "_mean"
	op_title = str(histo.GetTitle()) + " - MEAN"
  	   
	op_histo = ROOT.TH1F(op_name,op_title,x_bin_number,x_min,x_max)
	
	for x_bin_ite in range(1,x_bin_number+1):	
	     op_histo.SetBinContent(x_bin_ite,histo.ProjectionY(op_name+str(x_bin_ite),x_bin_ite,x_bin_ite,"o").GetMean())
	     op_histo.SetBinError(x_bin_ite,0.0001)
	     	     
	return op_histo

	        
		
def GetWidth(histo):	
  
  	x_bin_number = histo.GetNbinsX()
  	x_min = histo.GetXaxis().GetBinLowEdge(1)
  	x_max = histo.GetXaxis().GetBinUpEdge(x_bin_number)
		
	help = random.random()
	op_name = str(histo.GetName()) + "_width_" + str(help)
  	   
	op_histo = ROOT.TH1F(op_name,"",x_bin_number,x_min,x_max)	
	
	for x_bin_ite in range(1,x_bin_number+1):	
	     histo.GetXaxis().SetRange(x_bin_ite,x_bin_ite)
	     op_histo.SetBinContent(x_bin_ite,histo.Project3D("zo").GetRMS())
	     op_histo.SetBinError(x_bin_ite,0.0001)
	     	     
	return op_histo

	        
		
def GetWidth2(histo):	
  
  	x_bin_number = histo.GetNbinsX()
  	x_min = histo.GetXaxis().GetBinLowEdge(1)
  	x_max = histo.GetXaxis().GetBinUpEdge(x_bin_number)
		
	op_name = str(histo.GetName()) + "_width"
  	   
	op_histo = ROOT.TH1F(op_name,"",x_bin_number,x_min,x_max)	
	
	for x_bin_ite in range(1,x_bin_number+1):		
	     op_histo.SetBinContent(x_bin_ite,histo.ProjectionY(op_name+str(x_bin_ite),x_bin_ite,x_bin_ite,"o").GetRMS())
	     op_histo.SetBinError(x_bin_ite,0.0001)
	     	     
	return op_histo
	        
		
def GetFractionTH1Fs(num,den,kind):

        histo_output = ()
    
        for ass_ite in range(len(num)):
    
            op_histo_ass_axis = ()
	    
            for x_axis_ite in range(len(num[0])):
  
  	        x_bin_number = num[ass_ite][x_axis_ite].GetNbinsX()
  	        x_min = num[ass_ite][x_axis_ite].GetXaxis().GetBinLowEdge(1)
  	        x_max = num[ass_ite][x_axis_ite].GetXaxis().GetBinUpEdge(x_bin_number)
		
	        op_name = str(num[ass_ite][x_axis_ite].GetName()) + "_fraction"
	        op_title = str(num[ass_ite][x_axis_ite].GetName())
  	   
	        op_histo = ROOT.TH1F(op_name,op_title,x_bin_number,x_min,x_max)

  	        op_histo.GetXaxis().SetTitle(num[ass_ite][x_axis_ite].GetXaxis().GetTitle())

                for bin_ite in range(1,x_bin_number+1):
 	    
	            val = 0.
	            err = 0.
	
	            if den[ass_ite][x_axis_ite].GetBinContent(bin_ite)!=0:
	                if kind=="effic":
	                    val = num[ass_ite][x_axis_ite].GetBinContent(bin_ite)*1./den[ass_ite][x_axis_ite].GetBinContent(bin_ite)
		            err = math.sqrt( val*(1-val)/den[ass_ite][x_axis_ite].GetBinContent(bin_ite))
	                elif kind=="fakerate":
  	                    val = 1.- (num[ass_ite][x_axis_ite].GetBinContent(bin_ite)*1./den[ass_ite][x_axis_ite].GetBinContent(bin_ite))
		            err = math.sqrt( val*(1-val)/den[ass_ite][x_axis_ite].GetBinContent(bin_ite))
	                elif kind=="npu":
	                    val = num[ass_ite][x_axis_ite].GetBinContent(bin_ite)*1./den[ass_ite][x_axis_ite].GetBinContent(bin_ite)
			    err = 0.0001
			    
                    op_histo.SetBinContent(bin_ite, val)
                    op_histo.SetBinError(bin_ite, err)			
	    
  		op_histo_ass_axis+=op_histo,
		 
            histo_output+=op_histo_ass_axis,
		    
		    		    
	return  histo_output 
	        
		
def GetXYFraction(histos):

        histo_output = ()
	
	for ass_ite in range(len(histos)):
		
	    op_name = str(histos[ass_ite].GetName()) + "_xyfraction"
	
	    op_histo = ROOT.TH1F(op_name,op_name,50,0.,2.)
	
	    for x_bin_ite in range(1,histos[ass_ite].GetNbinsX()+1):
	        for y_bin_ite in range(1,histos[ass_ite].GetNbinsY()+1):
		    x_value = histos[ass_ite].GetBinCenter(x_bin_ite)
		    y_value = histos[ass_ite].GetBinCenter(y_bin_ite)
		    value = x_value*1./y_value
		    weight = histos[ass_ite].GetBinContent(x_bin_ite,y_bin_ite)
		    op_histo.Fill(value,weight)
		
	    histo_output+=op_histo,
	
	return histo_output
	        
		
def GetTh3FromTProfile3D(histo,min_bin,max_bin,name=""):
	
	op_name = str(name + histo.GetName()) + "_" + str(min_bin) + "to" + str(max_bin)
  
  	x_bin_number = histo.GetNbinsX()
  	x_min = histo.GetXaxis().GetBinLowEdge(1)
  	x_max = histo.GetXaxis().GetBinUpEdge(x_bin_number)
  
  	y_bin_number = histo.GetNbinsY()
  	y_min = histo.GetYaxis().GetBinLowEdge(1)
  	y_max = histo.GetYaxis().GetBinUpEdge(y_bin_number)
	
	if min_bin<0:
	    min_bin=0
	    
	if max_bin>histo.GetNbinsZ():
	    max_bin=histo.GetNbinsZ()
  
  	z_bin_number = 100
  	z_min = 0.
  	z_max = 2.
	
	op_histo = ROOT.TH3F(op_name,op_name, x_bin_number,x_min,x_max,y_bin_number,y_min,y_max,z_bin_number,z_min,z_max)
	
	for x_bin_ite in range(1,x_bin_number+1):
	    for y_bin_ite in range(1,y_bin_number+1):
		z_val_num = 0.
		z_val_denum = 0.
		entries = 0.
	        for z_bin_ite in range(min_bin,max_bin+1):
		    if(histo.GetBinContent(x_bin_ite,y_bin_ite,z_bin_ite)>0.):
			err = histo.GetBinError(x_bin_ite,y_bin_ite,z_bin_ite)
			if err==0.:
			    err=1.
			z_val_num+= histo.GetBinContent(x_bin_ite,y_bin_ite,z_bin_ite)*1./(err*err)
			z_val_denum+= 1./(err*err)
		        entries+= histo.GetBinEntries(histo.GetBin(x_bin_ite,y_bin_ite,z_bin_ite))
		z_val = 0.
		if z_val_denum>0.:
		    z_val = z_val_num*1./z_val_denum
	        z_bin = int((z_val*1./(z_max-z_min))*z_bin_number)
	        op_histo.SetBinContent(x_bin_ite,y_bin_ite,z_bin,entries)
		
	return op_histo
		    
		  
	
	
	        
		 
	        
		
def GetFractionTH2Fs(num,den,kind):

        histo_output = ()
    
        for ass_ite in range(len(num)):
  
  	    x_bin_number = num[ass_ite].GetNbinsX()
  	    x_min = num[ass_ite].GetXaxis().GetBinLowEdge(1)
  	    x_max = num[ass_ite].GetXaxis().GetBinUpEdge(x_bin_number)
  
  	    y_bin_number = num[ass_ite].GetNbinsY()
  	    y_min = num[ass_ite].GetYaxis().GetBinLowEdge(1)
  	    y_max = num[ass_ite].GetYaxis().GetBinUpEdge(y_bin_number)
		
	    help = random.random()
	    op_name = "fraction" + str(help)
	    op_title = num[ass_ite].GetName()
  	   
	    op_histo = ROOT.TH2F(op_name,op_title,x_bin_number,x_min,x_max,y_bin_number,y_min,y_max)

  	    op_histo.GetXaxis().SetTitle(num[ass_ite].GetXaxis().GetTitle())
  	    op_histo.GetYaxis().SetTitle(num[ass_ite].GetYaxis().GetTitle())

            for x_bin_ite in range(1,x_bin_number+1): 
		    
                for y_bin_ite in range(1,y_bin_number+1):
 	    
	            val = 0.
	            err = 0.
	
	            if den[ass_ite].GetBinContent(x_bin_ite,y_bin_ite)!=0:
  	                if kind=="effic":
	                    val = num[ass_ite].GetBinContent(x_bin_ite,y_bin_ite)*1./den[ass_ite].GetBinContent(x_bin_ite,y_bin_ite)
		            err = math.sqrt( val*(1-val)/den[ass_ite].GetBinContent(x_bin_ite,y_bin_ite))
	                elif kind=="fakerate":
  	                    val = 1.- (num[ass_ite].GetBinContent(x_bin_ite,y_bin_ite)*1./den[ass_ite].GetBinContent(x_bin_ite,y_bin_ite))
		            err = math.sqrt( val*(1-val)/den[ass_ite].GetBinContent(x_bin_ite,y_bin_ite))

                    op_histo.SetBinContent(x_bin_ite,y_bin_ite, val)
                    op_histo.SetBinError(x_bin_ite,y_bin_ite, err)	
		 
            histo_output+=op_histo,
		    
		    		    
	return histo_output
	        
	        
		
def DrawTH1Fs(histos,savename,maximize=False,permutize=False,pointer=False,legend=False):
	
    dim = len(histos)
    ass_number = len(histos[0])
    x_axis_number = len(histos[0][0])
	
    ROOT.SetRootPalette(1)
		
    help = random.random()
    canv_name = "canv" + str(help)
  
    canv = ROOT.TCanvas(canv_name,"Guck mal",10,10,1200,850)
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetTitleX(0.2)
    ROOT.gStyle.SetTitleY(0.93)
    ROOT.gStyle.SetTitleW(1.)
    if legend:
        ROOT.gStyle.SetOptStat("mr")
    else:
        ROOT.gStyle.SetOptStat(0)
    canv.Divide(x_axis_number,dim)
   
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetBorderMode(0)
    ROOT.gPad.SetFillColor(0)
    
    for dim_ite in range(dim):
   
        for x_axis_ite in range(x_axis_number):
		
	    pad_num = dim_ite*x_axis_number + x_axis_ite + 1
	    
	    if maximize:
	         maximum = histos[dim_ite][0][x_axis_ite].GetMaximum()
		 for ass_ite in range(ass_number):
	             if histos[dim_ite][ass_ite][x_axis_ite].GetMaximum()>maximum:
			 maximum=histos[dim_ite][ass_ite][x_axis_ite].GetMaximum()
		 histos[dim_ite][0][x_axis_ite].SetMaximum(maximum + 0.1*maximum)
		    

            canv.cd(pad_num)
            ROOT.gPad.SetGrid()
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetLeftMargin(0.2)
            ROOT.gPad.SetTopMargin(0.15)
            ROOT.gPad.SetBottomMargin(0.15)
            ROOT.gPad.SetFillColor(0)

            histos[dim_ite][0][x_axis_ite].GetXaxis().CenterTitle()
            histos[dim_ite][0][x_axis_ite].GetYaxis().CenterTitle()
  
            histos[dim_ite][0][x_axis_ite].GetXaxis().SetTitleOffset(1.1)
            histos[dim_ite][0][x_axis_ite].GetYaxis().SetTitleOffset(1.1)
	    color = my_colors[0]
	    if permutize: 
		color = my_colorsP[0]
            histos[dim_ite][0][x_axis_ite].SetLineColor(color)
            histos[dim_ite][0][x_axis_ite].SetLineWidth(2)
	    if legend:
                histos[dim_ite][0][x_axis_ite].SetStats(1)
	    else:
                histos[dim_ite][0][x_axis_ite].SetStats(0)
		
	    if pointer:			
	        histos[dim_ite][0][x_axis_ite].SetMarkerStyle(20)
	        histos[dim_ite][0][x_axis_ite].SetMarkerSize(0.8)
	        histos[dim_ite][0][x_axis_ite].SetMarkerColor(color)
		
            histos[dim_ite][0][x_axis_ite].Draw()
	    
	    if legend:
		moveStatBox(canv, histos[dim_ite][0][x_axis_ite], color, 0)

            for ass_ite in range(1,ass_number):
	        color = my_colors[ass_ite]
	        if permutize: 
   		    color = my_colorsP[ass_ite]
    		histos[dim_ite][ass_ite][x_axis_ite].SetLineColor(color)
                histos[dim_ite][ass_ite][x_axis_ite].SetLineWidth(2)
                histos[dim_ite][ass_ite][x_axis_ite].SetStats(1)		
	        if pointer:			
	            histos[dim_ite][ass_ite][x_axis_ite].SetMarkerStyle(20)
	            histos[dim_ite][ass_ite][x_axis_ite].SetMarkerSize(0.8)
	            histos[dim_ite][ass_ite][x_axis_ite].SetMarkerColor(color)
		
                histos[dim_ite][ass_ite][x_axis_ite].Draw("sames")
	    
	        if legend: 
 		    moveStatBox(canv, histos[dim_ite][ass_ite][x_axis_ite], color, ass_ite)
		
    canv.SaveAs("../pictures/"+savename+".png")
    canv.SaveAs("../root/"+savename+".root")
    canv.SaveAs("../pdfs/"+savename+".pdf")
    
    return
	        
	        
		
def DrawTH2Fs(histos,savename,opt,legend=False):
	
    dim = len(histos)
    ass_number = len(histos[0])
	
    ROOT.SetRootPalette(1)
		
    help = random.random()
    canv_name = "canv2D" + str(help)
  
    canv2D = ROOT.TCanvas(canv_name,"Guck mal",10,10,1200,850)
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetOptStat(1000000000)
    canv2D.Divide(ass_number,dim)
   
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetBorderMode(0)
    ROOT.gPad.SetFillColor(0)
    
    for dim_ite in range(dim):
   
        for ass_ite in range(ass_number):
		
	    pad_num = dim_ite*ass_number + ass_ite + 1
	
            canv2D.cd(pad_num)
            ROOT.gPad.SetGrid()
            ROOT.gPad.SetRightMargin(0.15)
            ROOT.gPad.SetLeftMargin(0.1)
            ROOT.gPad.SetFillColor(0)

            histos[dim_ite][ass_ite].GetXaxis().CenterTitle()
            histos[dim_ite][ass_ite].GetYaxis().CenterTitle()
            histos[dim_ite][ass_ite].GetZaxis().CenterTitle()
  
            histos[dim_ite][ass_ite].GetXaxis().SetTitleOffset(1.1)
            histos[dim_ite][ass_ite].GetYaxis().SetTitleOffset(1.14)
            histos[dim_ite][ass_ite].GetZaxis().SetTitleOffset(1.1)
		
	    if legend:
                histos[dim_ite][ass_ite].SetStats(1)
		   	
            histos[dim_ite][ass_ite].Draw(opt)
	    
	    if legend:
		moveStatBox(canv2D, histos[dim_ite][ass_ite], 1, 0, 2)
				
    canv2D.SaveAs("../pictures/"+savename+".png")
    canv2D.SaveAs("../root/"+savename+".root") 
    canv2D.SaveAs("../pdfs/"+savename+".pdf")
    
    return
        

		
def DrawProfs(profs,help_dump,savename):
	
    x_axis_number = len(profs)
    dir_number = len(profs[0])
    ass_number = len(profs[0][0])
	
    ROOT.SetRootPalette(1)
		
    help = random.random()
    canv_name = "canv" + str(help)
  
    canv = ROOT.TCanvas(canv_name,"Guck mal",10,10,1200,850)
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetTitleX(0.1)
    ROOT.gStyle.SetTitleY(0.94)
    ROOT.gStyle.SetTitleW(0.8)
    canv.Divide(x_axis_number,dir_number)
   
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetBorderMode(0)
    ROOT.gPad.SetFillColor(0)
    
    for dir_ite in range(dir_number):
   
        for x_axis_ite in range(x_axis_number):
		
	    pad_num = dir_ite*x_axis_number + x_axis_ite + 1	 
	    
	    canv.cd(pad_num)
            ROOT.gPad.SetGrid()
            ROOT.gPad.SetRightMargin(0.1)
            ROOT.gPad.SetLeftMargin(0.2)
            ROOT.gPad.SetTopMargin(0.15)
            ROOT.gPad.SetBottomMargin(0.15)
            ROOT.gPad.SetFillColor(0)

            help_dump[x_axis_ite].GetXaxis().CenterTitle()
            help_dump[x_axis_ite].GetYaxis().CenterTitle()
            help_dump[x_axis_ite].GetZaxis().CenterTitle()
  
            help_dump[x_axis_ite].GetXaxis().SetTitleOffset(1.1)
            help_dump[x_axis_ite].GetYaxis().SetTitleOffset(1.14)
            help_dump[x_axis_ite].SetStats(0)
	
	    help_dump[x_axis_ite].Draw()
	
	    for ass_ite in range(ass_number):
	
	        profs[x_axis_ite][dir_ite][ass_ite].SetMarkerStyle(20)
	        profs[x_axis_ite][dir_ite][ass_ite].SetMarkerSize(0.8)
                profs[x_axis_ite][dir_ite][ass_ite].SetMarkerColor(my_colorsP[ass_ite])
	        profs[x_axis_ite][dir_ite][ass_ite].Draw("same")
	
				
    canv.SaveAs("../pictures/"+savename+".png")
    canv.SaveAs("../root/"+savename+".root")
    canv.SaveAs("../pdfs/"+savename+".pdf")
    
    return
		    