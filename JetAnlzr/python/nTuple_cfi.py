import FWCore.ParameterSet.Config as cms
				
# nTuples production

dummy  = cms.EDProducer("CandViewNtpProducer",
	src = cms.InputTag("ak5GenJets"),
    	lazyParser = cms.untracked.bool(True),
    	prefix = cms.untracked.string(""),
    	eventInfo = cms.untracked.bool(True),
	variables = cms.VPSet(
		 cms.PSet(
		   tag = cms.untracked.string("pt"),
		   quantity = cms.untracked.string("pt")
		 ),
            	cms.PSet(
                   tag = cms.untracked.string("eta"),
                   quantity = cms.untracked.string("eta")
             	),
            	cms.PSet(
                   tag = cms.untracked.string("phi"),
                   quantity = cms.untracked.string("phi")
            	), 
	),     
)

def GetNTupsSequence(process,seq):
	
    process.ak5GenJetsn = dummy.clone()
	
    process.nTups = cms.Sequence()
    
    process.nTups.__iadd__(process.ak5GenJetsn)
    
    for mod_name in seq.moduleNames():
        
	inst_name = mod_name+"n"
	help = dummy.clone(
		           src = mod_name
	               )	
	help.setLabel(inst_name)  
	setattr(process, inst_name, help)
        process.nTups.__iadd__(help)
		
    return process.nTups
