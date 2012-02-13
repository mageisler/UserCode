// -*- C++ -*-
//
// Package:    TrackHitCount
// Class:      TrackHitCount
// 
/**\class TrackHitCount TrackHitCount.cc MGeisler/TrackHitCount/src/TrackHitCount.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Thu Jan 12 11:40:42 CET 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

//C++ include files
#include <vector>
#include <string>
#include <numeric>
#include <map>
#include <list>
#include <cctype>

//ROOT include files
#include "TCanvas.h"
#include "TChain.h"
#include "TProfile2D.h"
#include "TGaxis.h"
#include "TStyle.h"

using namespace edm;
using namespace std;
using namespace reco;
//
// class declaration
//

class TrackHitCount : public edm::EDAnalyzer {
   public:
      explicit TrackHitCount(const edm::ParameterSet&);
      ~TrackHitCount();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      string fileName_;
      TProfile2D* tkhitprof;
      char thpname[256];
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
TrackHitCount::TrackHitCount(const edm::ParameterSet& iConfig)

{	
   	//now do what ever initialization is needed

	fileName_ = iConfig.getUntrackedParameter<string>("FileName");
	sprintf(thpname,"%s_TrackerHitProfile",fileName_.c_str());

    	Service<TFileService> tfserv;

  	tkhitprof = tfserv->make<TProfile2D>(thpname,"Layers hit vs #phi vs #eta; #phi; #eta",30,-M_PI,M_PI,30,-3.,3.);

}


TrackHitCount::~TrackHitCount()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackHitCount::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    	Service<TFileService> tfserv;
    
	//get the input track collection     
  	Handle<TrackCollection> trkcollH;
  	iEvent.getByLabel("generalTracks",trkcollH);
	  	
	//loop over all tracks in the track collection	
  	for(TrackCollection::const_iterator trk_ite=trkcollH->begin(); trk_ite!=trkcollH->end(); ++trk_ite){

 	  if(trk_ite->pt()<0.3) continue;

	  unsigned validHit_counter = 0;

	  for(trackingRecHit_iterator hit_ite=trk_ite->recHitsBegin(); hit_ite!=trk_ite->recHitsEnd(); ++hit_ite){

	    if((*hit_ite)->isValid()) validHit_counter++;

	  }

	  tkhitprof->Fill(trk_ite->phi(),trk_ite->eta(),validHit_counter);

	}

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackHitCount::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackHitCount::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TrackHitCount::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackHitCount::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackHitCount::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackHitCount::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackHitCount::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackHitCount);
