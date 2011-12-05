#ifndef PFCand_NoPU_WithAM_h
#define PFCand_NoPU_WithAM_h

// -*- C++ -*-
//
// Package:    PFCand_NoPU_WithAM
// Class:      PFCand_NoPU_WithAM
// 
/**\class PF_PU_AssoMap PFCand_NoPU_WithAM.cc MGeisler/PF_PU_AssoMap/plugins/PFCand_NoPU_WithAM.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Thu Dec  1 16:07:41 CET 2011
// $Id$
//
//

#include "MGeisler/PF_PU_AssoMap/interface/PFCand_NoPU_WithAM_Algos.h"

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

//
// class declaration
//

class PFCand_NoPU_WithAM : public edm::EDProducer {
   public:
      explicit PFCand_NoPU_WithAM(const edm::ParameterSet&);
      ~PFCand_NoPU_WithAM();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual bool TrackMatch(reco::TrackRef,reco::TrackRef);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::InputTag input_PFCandidates_;
      edm::InputTag input_VertexCollection_;
      edm::InputTag input_VertexTrackAssociationMap_;

      edm::InputTag ConversionsCollection_;

      edm::InputTag KshortCollection_;
      edm::InputTag LambdaCollection_;

      edm::InputTag NIVertexCollection_;
};


#endif
