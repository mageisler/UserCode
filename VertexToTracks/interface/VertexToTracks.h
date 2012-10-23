#ifndef VertexToTracks_h
#define VertexToTracks_h

/**\class VertexToTracks VertexToTracks.cc MGeisler/VertexToTracks/plugins/VertexToTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Mon Oct 22 16:12:39 CEST 2012
// $Id$
//
//

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MGeisler/VertexToTracks/interface/VertexToTracksAlgos.h"

//
// class declaration
//

class VertexToTracks : public edm::EDProducer, private VertexToTracksAlgos {
   public:
      explicit VertexToTracks(const edm::ParameterSet&);
      ~VertexToTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag input_TrackCollection_;

      edm::InputTag input_VertexCollection_;

      int input_MaxNumAssociations_;

};


#endif
