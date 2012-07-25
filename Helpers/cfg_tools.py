import FWCore.ParameterSet.Config as cms

def RemoveModules(process,seq,lastMod):
	
    process.output = cms.Sequence()

    mod_list = str(seq).split("+")
    mods_length = len(mod_list)
    
    modsToRemove = ()
    
    for mod in mod_list:
	   
       if mod==lastMod:
	   process.output.__iadd__(getattr(process,mod))
           break
       else: 	   
           process.output.__iadd__(getattr(process,mod))
    
    return process.output
 