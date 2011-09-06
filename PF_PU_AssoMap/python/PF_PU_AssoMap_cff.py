import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQM_cfg import *
DQMStore.collateHistograms =cms.untracked.bool(True)
from MGeisler.PF_PU_AssoMap.pf_pu_assomap_cfi import *
