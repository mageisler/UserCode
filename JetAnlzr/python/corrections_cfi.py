import FWCore.ParameterSet.Config as cms
				
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *

corrDict = {"L123":cms.vstring('ak5chsChainL123'), "L23":cms.vstring('ak5chsChainL23'), "L123Residual":cms.vstring('ak5chsChainL123Residual'), "L23Residual":cms.vstring('ak5chsChainL23Residual')}
		    
# dummy
dummy = cms.EDProducer('PFJetCorrectionProducer',
    src        = cms.InputTag(''),
    correctors = cms.vstring('')
)

def GetCorrSequence(process,jetC,corr):

    process.ak5chsL1Producer = cms.ESProducer('L1OffsetCorrectionESProducer',
        level = cms.string('L1Offset'),
        algorithm = cms.string('AK5PFchs'),
        vertexCollection = cms.string('selectedPrimaryVertexQuality'),
        minVtxNdof = cms.int32(4),
    )
    process.ak5PFchsL2Relative = cms.ESProducer('LXXXCorrectionESProducer',
        level     = cms.string('L2Relative'),
        algorithm = cms.string('AK5PFchs')
    )
    process.ak5PFchsL3Absolute = cms.ESProducer('LXXXCorrectionESProducer',
        level     = cms.string('L3Absolute'),
        algorithm = cms.string('AK5PFchs')
    )
    process.ak5PFchsResidual = cms.ESProducer('LXXXCorrectionESProducer',
        level     = cms.string('L2L3Residual'),
        algorithm = cms.string('AK5PFchs')
    )	
    process.ak5chsChainL123  = cms.ESProducer('JetCorrectionESChain',
        correctors = cms.vstring('ak5chsL1Producer','ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )
    process.ak5chsChainL23 = cms.ESProducer('JetCorrectionESChain',
        correctors = cms.vstring('ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )
    process.ak5chsChainL123Residual  = cms.ESProducer('JetCorrectionESChain',
        correctors = cms.vstring('ak5chsL1Producer','ak5PFchsL2Relative','ak5PFchsL3Absolute','ak5PFchsResidual')
    )
    process.ak5chsChainL23Residual  = cms.ESProducer('JetCorrectionESChain',
        correctors = cms.vstring('ak5PFchsL2Relative','ak5PFchsL3Absolute','ak5PFchsResidual')
    )
    
    process.corr = cms.Sequence()
    
    for jet_name in jetC.moduleNames():
    
        for corr_name in corr:
		
            correction = corrDict[corr_name]
		
	    help = dummy.clone(
	               src = jet_name,
		       correctors = correction,
		   )
		             
	    help.setLabel(str(jet_name)+corr_name)  
	    setattr(process, str(jet_name)+corr_name, help)
            process.corr.__iadd__(help)
		
    return process.corr
