//// -*- C++ -*-
//
// Package:    FirstVertexTracks
// Class:      FirstVertexTracks
// 
/**\class FirstVertexTracks FirstVertexTracks.cc RecoParticleFlow/FirstVertexTracks/src/FirstVertexTracks.cc

 Description: produces a TrackCollection of tracks associated to the first vertex

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Mon Apr 11 17:36:26 CEST 2011
// $Id: FirstVertexTracks.cc,v 1.1 2011/04/19 13:39:32 mgeisler Exp $
//
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


//
// class declaration
//
   
using namespace edm;
using namespace std;
using namespace reco;

class FirstVertexTracks : public edm::EDProducer {
   public:
      explicit FirstVertexTracks(const edm::ParameterSet&);
      ~FirstVertexTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;
     typedef vector<pair<TrackRef, float> > TrackQualityPairVector;
     typedef pair<TrackRef, float> TrackQualityPair;

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      InputTag input_VertexTrackAssociationMap_;
      InputTag input_VertexCollection_;
      bool input_generalTracksCollection_;
      bool input_GsfElectronCollection_;
      bool input_AssociationQuality_;
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
FirstVertexTracks::FirstVertexTracks(const edm::ParameterSet& iConfig)
{
   //now do what ever other initialization is needed

  	input_VertexTrackAssociationMap_ = iConfig.getParameter<InputTag>("VertexTrackAssociationMap");

  	input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");

  	input_generalTracksCollection_= iConfig.getUntrackedParameter<bool>("TrackCollection", true);

  	input_GsfElectronCollection_= iConfig.getUntrackedParameter<bool>("GsfElectronCollection", true);

  	input_AssociationQuality_= iConfig.getUntrackedParameter<bool>("AssociationQuality", true);


	produces<TrackCollection>();
  
}


FirstVertexTracks::~FirstVertexTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FirstVertexTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	auto_ptr<TrackCollection> firstvertextracks(new TrackCollection() );
  
	//get the input vertex<->track association map
  	Handle<TrackVertexAssMap> assomap;
  	iEvent.getByLabel(input_VertexTrackAssociationMap_,assomap);
 
	//get the input vertex collection
  	Handle<VertexCollection> vtxcoll;
  	iEvent.getByLabel(input_VertexCollection_,vtxcoll);

	//get the first vertex
        const VertexRef firstvertexref(vtxcoll,0);
	
	//look if an input for generalTrackCollection is set
 	if (input_generalTracksCollection_){

	  //loop over all vertices in the association map
          for (TrackVertexAssMap::const_iterator assomap_ite = assomap->begin(); assomap_ite != assomap->end(); assomap_ite++){

	    const VertexRef assomap_vertexref = assomap_ite->key;

  	    //take only the first vertex from the association map 
	    if ((assomap_vertexref)==(firstvertexref)){  
 
  	      const TrackQualityPairVector trckcoll = assomap_ite->val;

	      //get the tracks associated to the first vertex and store them in a track collection
	      for (unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){

		//since every general Track has a quality >= -1. no loop over all general Tracks necessary
	 	if (trckcoll[trckcoll_ite].second>=-1.)  
	          if (!(input_AssociationQuality_) || (trckcoll[trckcoll_ite].second>0.)) firstvertextracks->push_back(*(trckcoll[trckcoll_ite].first));

	      }

 	    }

	  }
	
	}
	
	//look if an input for GsfElectronCollection is set
 	if (input_GsfElectronCollection_!){	

	  //loop over all vertices in the association map
          for (TrackVertexAssMap::const_iterator assomap_ite = assomap->begin(); assomap_ite != assomap->end(); assomap_ite++){

 	    const VertexRef assomap_vertexref = assomap_ite->key;

	    //take only the first vertex from the association map 
	    if ((assomap_vertexref)==(firstvertexref)){  
 
	      const TrackQualityPairVector trckcoll = assomap_ite->val;

	      //get the tracks associated to the first vertex and store them in a track collection
	      for (unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){

		//since only a GsfElectron has a quality = -2. no loop over all GsfElectrons necessary
	        if (trckcoll[trckcoll_ite].second==-2.) firstvertextracks->push_back(*(trckcoll[trckcoll_ite].first));

	      }

            }

 	  }

	}

	iEvent.put( firstvertextracks );
 
}


// ------------ method called once each job just before starting event loop  ------------
void 
FirstVertexTracks::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FirstVertexTracks::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
FirstVertexTracks::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FirstVertexTracks::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FirstVertexTracks::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FirstVertexTracks::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FirstVertexTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FirstVertexTracks);