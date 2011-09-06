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
// $Id: FirstVertexTracks.cc,v 1.6 2011/06/01 13:15:47 mgeisler Exp $
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;
  typedef AssociationMap<OneToManyWithQuality< VertexCollection, GsfElectronCollection, float> > GsfVertexAssMap;

  typedef vector<pair<TrackRef, float> > TrackQualityPairVector;
  typedef pair<TrackRef, float> TrackQualityPair;

  typedef vector<pair<GsfElectronRef, float> > GsfQualityPairVector;
  typedef pair<GsfElectronRef, float> GsfQualityPair;


//
// class declaration
//

class FirstVertexTracks : public edm::EDProducer {
   public:
      explicit FirstVertexTracks(const edm::ParameterSet&);
      ~FirstVertexTracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual bool TrackMatch(reco::TrackRef,reco::TrackRef);
      virtual bool GsfTrackMatch(reco::GsfElectronRef,reco::GsfElectronRef);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      InputTag input_VertexTrackAssociationMap_;
      InputTag input_VertexCollection_;
      string input_generalTracksCollection_;
      string input_GsfElectronCollection_;
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

  	input_generalTracksCollection_= iConfig.getUntrackedParameter<string>("TrackCollection","default");

  	input_GsfElectronCollection_= iConfig.getUntrackedParameter<string>("GsfElectronCollection","default");


// 	if (input_generalTracksCollection_!="default") produces<TrackCollection>();
	if (input_GsfElectronCollection_!="default") produces<GsfElectronCollection>();
  	produces<TrackCollection>();
// 	produces<TrackCollection>("VCTracks");

}


FirstVertexTracks::~FirstVertexTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// ------------ method called once each job just before starting event loop  ------------
void 
FirstVertexTracks::beginJob()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FirstVertexTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	
	//look if an input for generalTrackCollection is set
 	if (input_generalTracksCollection_!="default"){

	  auto_ptr<TrackCollection> gT_firstvertextracks(new TrackCollection() );
  
	  //get the input vertex<->general track association map
  	  Handle<TrackVertexAssMap> GTassomap;
  	  iEvent.getByLabel(input_VertexTrackAssociationMap_,GTassomap);

	  const VertexRef assomap_vertexref = GTassomap->begin()->key;
	  const TrackQualityPairVector trckcoll = GTassomap->begin()->val;

	  TrackRef GTtrackref;

	  //get the tracks associated to the first vertex and store them in a track collection
	  for (unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){

	    bool isGtMatched = false;
	    GTtrackref = trckcoll[trckcoll_ite].first;
		
	    //if input Collection is not generalTracks loop over all tracks of the input collection to find the correct general Track
	    if (input_generalTracksCollection_!="generalTracks"){
 
	      //get the input track collection
  	      Handle<TrackCollection> input_trckcoll;
  	      iEvent.getByLabel(input_generalTracksCollection_,input_trckcoll);
	
	      TrackRef input_trackref;
	      unsigned index_input_trck = 0;

  	      for(TrackCollection::const_iterator input_trckcoll_ite = input_trckcoll->begin(); input_trckcoll_ite!=input_trckcoll->end(); input_trckcoll_ite++,index_input_trck++){

		input_trackref = TrackRef(input_trckcoll,index_input_trck);

   	  	if(TrackMatch(GTtrackref,input_trackref)){

		  isGtMatched = true;
		  break;
	 	   	      
		} 

	      }

	    } else isGtMatched = true;

	    if (isGtMatched) gT_firstvertextracks->push_back(*GTtrackref);

	  }

	  iEvent.put( gT_firstvertextracks );
	
	}
 
	//get the input vertex collection
  	Handle<VertexCollection> vtxcoll;
  	iEvent.getByLabel(input_VertexCollection_,vtxcoll);

	//get the first vertex
        const VertexRef firstvertexref(vtxcoll,0);

        //create TrackCollection for the VertexCollection
	auto_ptr<TrackCollection> VCtrks(new TrackCollection() );;
	
	//get first vertex from the VertexCollection
        const VertexRef firstvertexrefVC(vtxcoll,0);

     	typedef reco::Vertex::trackRef_iterator IT;

        for(IT iTrack=(*firstvertexrefVC).tracks_begin(); iTrack!=(*firstvertexrefVC).tracks_end(); ++iTrack) {
	 
      	  const reco::TrackBaseRef& vertexbaseRef = *iTrack;

      	  float w = (*firstvertexrefVC).trackWeight(vertexbaseRef);
      	  if (w>0.) VCtrks->push_back((**iTrack));
 
        }

// 	iEvent.put( VCtrks,"VCTracks" );

	//look if an input for GsfElectronCollection is set
 	if (input_GsfElectronCollection_!="default"){
  
	  //get the input vertex<->gsf electron track association map
  	  Handle<GsfVertexAssMap> Gsfassomap;
  	  iEvent.getByLabel(input_VertexTrackAssociationMap_,Gsfassomap);
	
	  auto_ptr<GsfElectronCollection> gsf_firstvertextracks(new GsfElectronCollection() );	

	  //loop over all vertices in the association map
          for (GsfVertexAssMap::const_iterator Gsfassomap_ite = Gsfassomap->begin(); Gsfassomap_ite != Gsfassomap->end(); Gsfassomap_ite++){

 	    const VertexRef assomap_vertexref = Gsfassomap_ite->key;

	    //take only the first vertex from the association map 
	    if ((assomap_vertexref)==(firstvertexref)){  
 
	      const GsfQualityPairVector gsfcoll = Gsfassomap_ite->val;

	      GsfElectronRef Gsfelectronref;

	      //get the tracks associated to the first vertex and store them in a track collection
	      for (unsigned int gsfcoll_ite = 0; gsfcoll_ite < gsfcoll.size(); gsfcoll_ite++){

		bool isGsfMatched = false;
		Gsfelectronref = gsfcoll[gsfcoll_ite].first;
		
		//if input Collection is not gsfElectrons loop over all tracks of the input collection to find the correct gsfElectron
	  	if (input_GsfElectronCollection_!="gsfElectrons"){	
 
	  	  //get the input track collection
  	  	  Handle<GsfElectronCollection> input_gsfcoll;
  	  	  iEvent.getByLabel(input_GsfElectronCollection_,input_gsfcoll);
	          
		  GsfElectronRef input_gsfref;
		  unsigned index_input_gsf = 0;

  	  	  for(GsfElectronCollection::const_iterator input_gsfcoll_ite = input_gsfcoll->begin(); input_gsfcoll_ite!=input_gsfcoll->end(); input_gsfcoll_ite++,index_input_gsf++){

		    input_gsfref = GsfElectronRef(input_gsfcoll,index_input_gsf);

   	  	    if(GsfTrackMatch(Gsfelectronref,input_gsfref)){

		      isGsfMatched = true;
		      break;
	 	   	      
		    } 

		  }

	  	} else isGsfMatched = true;

		if (isGsfMatched) gsf_firstvertextracks->push_back(*(gsfcoll[gsfcoll_ite].first));

	      }

            }

 	  }

	  iEvent.put( gsf_firstvertextracks );

	}
 
}

bool 
FirstVertexTracks::TrackMatch(reco::TrackRef trackref1,reco::TrackRef trackref2)
{

	return (
	  (*trackref1).eta() == (*trackref2).eta() &&
	  (*trackref1).phi() == (*trackref2).phi() &&
	  (*trackref1).chi2() == (*trackref2).chi2() &&
	  (*trackref1).ndof() == (*trackref2).ndof() &&
	  (*trackref1).p() == (*trackref2).p()
	);

}

bool 
FirstVertexTracks::GsfTrackMatch(reco::GsfElectronRef gsfref1,reco::GsfElectronRef gsfref2)
{

	return (
	  (*gsfref1).eta() == (*gsfref2).eta() &&
	  (*gsfref1).phi() == (*gsfref2).phi() &&
	  (*gsfref1).p() == (*gsfref2).p()
	);

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