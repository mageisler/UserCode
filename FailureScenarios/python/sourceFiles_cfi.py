import FWCore.ParameterSet.Config as cms
import os, fcntl, fcntl, select, sys,subprocess
from numpy import loadtxt, savetxt, size
from YKuessel.TopCharge.RemoveDuplicates import RemoveDuplicates
       
def GetFileNames(abbrev, filename):
    print "Generate the file names to be used.\n Read in the ProdDetails..."
    dictOfLists = ReadProdDetails(filename)
    print dictOfLists.keys()
    listOfFilenames=cms.untracked.vstring()
    lines=""
    for sample, item in enumerate(dictOfLists["abbreviation"]):
        if item==abbrev:
            print "...sample found. now looking for the filenames in directory given by the productionfile summary. merged exists? "
            print dictOfLists["merged"][sample]
            print sample
            print item
            if (dictOfLists["merged"])[sample]=="-1":
                print " command: ls /pnfs/physik.rwth-aachen.de/cms"+(dictOfLists["output"])[sample]
                p = subprocess.Popen(["uberftp","grid-ftp","ls /pnfs/physik.rwth-aachen.de/cms"+(dictOfLists["output"])[sample]], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in p.stdout.readlines():
                    filename = line.strip().split()
                    #print filename
                    if len(filename)>8:
                        if (dictOfLists["outputFilename"])[sample].strip() in filename[8]:
                            #listOfFilenames.append((dictOfLists["output"])[sample]+"/"+filename[8])
  			    lines += (dictOfLists["output"])[sample]+filename[8] + "\n"
            else:
                print "Merged file exists and will be used."
                listOfFilenames.append("'"+(dictOfLists["merged"])[sample]+"'")
                lines += (dictOfLists["merged"])[sample]
		return listOfFilenames
    print "input file list generated. will be saved in tmp file "
    filenameTmp='listOfOutputFilesTmp'+abbrev
    print filenameTmp
    f = file(filenameTmp, 'w')
    f.write(lines)
    f.close()
    print "duplicates will be removed"
    RemoveDuplicates(filenameTmp, filenameTmp+'NoDuplicates')
    g=file(filenameTmp+'NoDuplicates')
    listOfFilenamesNoDuplicates=cms.untracked.vstring()
    for line in g.readlines():   
  	  listOfFilenamesNoDuplicates.append(line)
    print "GetFileNames done"
    return listOfFilenamesNoDuplicates

def ReadProdDetails(filename):
    print " \n  .. read in production details... \n "
    M=loadtxt(filename,  dtype='S')
    
    dictOfLists = {}
    for i in range(0,len(M)): #zeilen
        if M[i][0]=="Sample":
            for j in range(0,len(M[i])): #spalten                        
                list=[]
                for k in range(i+1,len(M)): #zeilen
                    list.append(M[k][j])
                dictOfLists[M[i][j]]=list

    return dictOfLists
