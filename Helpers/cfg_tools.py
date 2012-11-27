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



def changeString(oldInputParameter, oldSuffix, newSuffix, dict):
	
    for key in dict:
        if (key in str(oldInputParameter)):
	    new_input = cms.string(dict[key])	
            return new_input
	  
	
    if (oldSuffix in str(oldInputParameter)):	
			    
	help1 = str(oldInputParameter).split(oldSuffix)[0]
	help2 = help1.split("(")[1]
	help3 = "(" + help2[0]
	if len(help2)>1:
	    help3 = help2[1]				
	help4 = help3[1:len(help3)] + newSuffix
	new_input = cms.string(help4)
	
    else:	    
        new_input = oldInputParameter
	
    return new_input



def changeInputTag(oldInputParameter, oldSuffix, newSuffix, dict):
	
    for key in dict:
        if (key in str(oldInputParameter)):
	    new_input = cms.InputTag(dict[key])	
            return new_input
	
    if (oldSuffix in str(oldInputParameter)):	
			    
	help1 = str(oldInputParameter).split(oldSuffix)[0]
	help2 = help1.split("(")
	help3 = "(" + help2[0]
	if len(help2)>1:
	    help3 = help2[1]		
	help4 = help3[1:len(help3)] + newSuffix
	new_input = cms.InputTag(help4)
	
    else:	    
        new_input = oldInputParameter
	
    return new_input
    
	
def changeVInputTag(oldInputParameter, oldSuffix, newSuffix, dict):
			       
    new_input = cms.VInputTag() 
		    
    if "NoneType" in str(type(oldInputParameter)):		        
        return new_input
			    			
    for old_Input in oldInputParameter:	    	 
	    
        new_Input = changeInputTag(old_Input,oldSuffix,newSuffix,dict) 				    
        new_input.append(new_Input)
    
    return new_input
	    

	
def changePSet(oldInputParameter, oldSuffix, newSuffix, dict):
			       
    new_input = cms.PSet() 
		    
    for para in oldInputParameter.parameters_():
	    
        name = oldInputParameter.getParameter(para)	    
	input_type = type(name)
	
	if ('str' in str(input_type)):
	    
            new_Input = changeString(name,oldSuffix,newSuffix,dict) 				    
            setattr(new_input, str(para), new_Input)
	    continue
	
	elif ('.PSet' in str(input_type)):
	    
            new_Input = changePSet(name,oldSuffix,newSuffix,dict) 				    
            setattr(new_input, str(para), new_Input)
	    continue
	
	elif ('.VPSet' in str(input_type)):
	    
            new_Input = changeVPSet(name,oldSuffix,newSuffix,dict) 				    
            setattr(new_input, str(para), new_Input)
	    continue
	
	else:
	    
            new_Input = changeString(name,oldSuffix,newSuffix,dict) 				    
            setattr(new_input, str(para), name)
	    continue
    
    return new_input
	
	
def changeVPSet(oldInputParameter, oldSuffix, newSuffix, dict):
			       
    new_input = cms.VPSet() 		    
			
    for old_Input in oldInputParameter:	  
	      
	input_type = str(type(old_Input))
	
	if ('.PSet' in input_type):
	
            new_Input = changePSet(old_Input, oldSuffix, newSuffix,dict)
	    continue
	
	elif ('.VPSet' in input_type):
	
            new_Input = changeVPSet(old_Input, oldSuffix, newSuffix,dict)
	    continue
 		
        new_input.append(new_Input)
    
    return new_input
       
     
     
def clonePath(process,path,firstModule,oldSuffix,newSuffix,dict):
	
    output_path = cms.Path()

    module_list = str(path).split("+")
    
    alreadyInserted = ()
    
    modToChange = False
    
    for module in module_list:
	    
	if (module=="None") or (module in alreadyInserted):
	    continue
		
	temp = getattr( process, module ).clone()
	    
        if (module == firstModule) and (modToChange==False):
	    modToChange = True
	    
	if modToChange:	 
	       
	    inputParameters = temp.parameters_()
	    
	    for tag in inputParameters:
		
		oldInputParameter = inputParameters[tag]
		input_type = str(type(temp.getParameter(tag)))			
						
		if ( "Types.string" in input_type  ):
			       
	            new_input = changeString(oldInputParameter,oldSuffix,newSuffix,dict)
			   		
		    setattr(temp, tag, new_input)
		    continue			
						
		elif ( "Types.VInputTag" in input_type  ):
			       
	            new_input = changeVInputTag(oldInputParameter,oldSuffix,newSuffix,dict)
			   		
		    setattr(temp, tag, new_input)
		    continue	
			    
		elif ( "Types.InputTag" in input_type ):
		    
		    new_input = changeInputTag(oldInputParameter,oldSuffix,newSuffix,dict) 
			   		
		    setattr(temp, tag, new_input)
		    continue					
						
		elif ( "Types.VPSet" in input_type ):
			       
	            new_input = changeVPSet(oldInputParameter,oldSuffix,newSuffix,dict) 
			   		
		    setattr(temp, tag, new_input)
		    continue
			    
		elif ( "Types.PSet" in input_type ):
			       
	            new_input = changePSet(oldInputParameter,oldSuffix,newSuffix,dict) 
			   		
		    setattr(temp, tag, new_input)
		    continue		
		    	    
	new_module_label = module.replace(oldSuffix,newSuffix)
	     
	temp.setLabel(new_module_label)  
	setattr(process, new_module_label, temp)   
	    
        output_path.__iadd__( temp )
	alreadyInserted+= module,
	    
    return output_path