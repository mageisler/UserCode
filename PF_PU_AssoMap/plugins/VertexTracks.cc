// -*- C++ -*-
//
// Package:    VertexTracks
// Class:      VertexTracks
// 
/**\class VertexTracks VertexTracks.cc Validation/VertexTracks/src/VertexTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Tue Oct 18 16:24:26 CEST 2011
// $Id$
//
//
#include "MGeisler/PF_PU_AssoMap/interface/PF_PU_AssoMap.h"


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef vector<VertexRef > VertexRefV;


//
// class declaration
//

class VertexTracks : public edm::EDProducer {
   public:
      explicit VertexTracks(const edm::ParameterSet&);
      ~VertexTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      InputTag input_VertexCollection_;
      string input_generalTracksCollection_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
VertexTracks::VertexTracks(const edm::ParameterSet& iConfig)
{
   //register your products
   	produces<TrackCollection>();

   //now do what ever other initialization is needed

  	input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");

  	input_generalTracksCollection_= iConfig.getUntrackedParameter<string>("TrackCollection","default");
  
}


VertexTracks::~VertexTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
VertexTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	auto_ptr<TrackCollection> firstvertextracks(new TrackCollection() );

        Handle<BeamSpot> beamSpot;
     	iEvent.getByLabel("offlineBeamSpot",beamSpot); 
 
	//get the input track collection
  	Handle<TrackCollection> trkcollH;
  	iEvent.getByLabel(input_generalTracksCollection_,trkcollH);

	//get the input vertex collection
  	Handle<VertexCollection> vtxcollH;
  	iEvent.getByLabel(input_VertexCollection_,vtxcollH);
 	VertexRefV* qualvtxcoll = PF_PU_AssoMapAlgos::QualifiedVertices(vtxcollH,true,4.);
	
	//get first vertex from the qualified VertexCollection
        const VertexRef firstvertexref = qualvtxcoll->at(0);
	  	
	//loop over all tracks in the track collection	
  	for(unsigned int index_trck = 0;  index_trck < trkcollH->size(); ++index_trck) {

 	  TrackRef input_trackref(trkcollH,index_trck);
	  const TrackBaseRef& trackbaseRef = TrackBaseRef(input_trackref);

	  float w  = firstvertexref->trackWeight(trackbaseRef);
	      	  
      	  if ((w>0.) &&
	      (input_trackref->hitPattern().trackerLayersWithMeasurement() >= 3) &&
	      (input_trackref->hitPattern().pixelLayersWithMeasurement() +
	       input_trackref->hitPattern().numberOfValidStripLayersWithMonoAndStereo() >= 0) &&
	      (fabs(input_trackref->pt()) >= 0.1) &&
	      (input_trackref->eta() >= -5.0) && 
	      (input_trackref->eta() <= 5.0) &&
	      (fabs(input_trackref->dxy(beamSpot->position())) <= 120.0) &&
	      (fabs(input_trackref->dsz(beamSpot->position())) <= 300.0)  &&
	      (input_trackref->normalizedChi2()<=10000.0) &&
	      (input_trackref->qualityByName("highPurity"))) firstvertextracks->push_back((*input_trackref));
 
        }

	iEvent.put( firstvertextracks );
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
VertexTracks::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VertexTracks::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
VertexTracks::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
VertexTracks::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
VertexTracks::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
VertexTracks::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexTracks);
