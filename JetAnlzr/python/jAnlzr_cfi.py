import FWCore.ParameterSet.Config as cms
import sys
				
# jet analyzer
    
dummy = cms.EDAnalyzer('JetAnlzr',
	genJets = cms.string("ak5GenJetsn"),
	recoJets = cms.vstring(),
	PileUpInfo = cms.string("addPileupInfo"),
	Outname = cms.string("")
)

def GetJAnlzrSequence(process,seq,corr):
	
    process.janlzr = cms.Sequence()
    
    janlzr_num = len(corr) + 1
    
    input_liste = seq.moduleNames()
    liste = [()]*janlzr_num
	
    for mod_name in input_liste:
	    
	appended = False
    
        for corr_ite in range(janlzr_num-1):
	
	    if corr[corr_ite] in str(mod_name):
              
	        liste[corr_ite]+= mod_name+"n",
		appended = True
		
	if not appended:
	    liste[janlzr_num-1]+= mod_name+"n",
	    
    corr+=["Uncorrected"]
	    
    for list_ite in range(janlzr_num):
        
	inst_name = "JetAnlzr"+corr[list_ite]
	
	reco_jets_help = cms.vstring()
	
	for jetC in liste[list_ite]:
	    reco_jets_help.append(jetC)
	    
	help = dummy.clone(
	           recoJets = reco_jets_help,
		   Outname = inst_name,
	       )
	       
	help.setLabel(inst_name)  
	setattr(process, inst_name, help)
        process.janlzr.__iadd__(help)
		
    return process.janlzr