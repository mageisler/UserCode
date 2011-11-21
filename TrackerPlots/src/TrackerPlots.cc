// -*- C++ -*-
//
// Package:    TrackerPlots
// Class:      TrackerPlots
// 
/**\class TrackerPlots TrackerPlots.cc MGeisler/TrackerPlots/src/TrackerPlots.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Fri Nov 18 16:10:56 CET 2011
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"


#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterfwd.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollectionfwd.h"
#include "CommonTools/TrackerMap/interface/TrackerMap.h"

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
#include "TH2F.h"
#include "TGaxis.h"
#include "TStyle.h"

using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

class TrackerPlots : public edm::EDAnalyzer {
   public:
      explicit TrackerPlots(const edm::ParameterSet&);
      ~TrackerPlots();

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
      TrackerMap* tkmap;
      char tmname[256];
      char savename[256];
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
TrackerPlots::TrackerPlots(const edm::ParameterSet& iConfig)

{	
   	//now do what ever initialization is needed

	fileName_ = iConfig.getUntrackedParameter<string>("FileName");
	sprintf(tmname,"%s_TrackerMap",fileName_.c_str());
	tkmap = new TrackerMap(tmname);

}


TrackerPlots::~TrackerPlots()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackerPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
	//get the input track collection     
  	Handle<TrackCollection> trkcollH;
  	iEvent.getByLabel("generalTracks",trkcollH);
	  	
	//loop over all tracks in the track collection	
  	for(TrackCollection::const_iterator trk_ite=trkcollH->begin(); trk_ite!=trkcollH->end(); ++trk_ite){

	  for(trackingRecHit_iterator hit_ite=trk_ite->recHitsBegin(); hit_ite!=trk_ite->recHitsEnd(); ++hit_ite){

	    if((*hit_ite)->geographicalId().det() == DetId::Tracker) tkmap->fill((*hit_ite)->geographicalId(),1.);

	  }

	}

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrackerPlots::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackerPlots::endJob() 
{
          sprintf(savename,"%s_TrackerMap.png",fileName_.c_str());
	  tkmap->save(true,0,0,savename); 
}

// ------------ method called when starting to processes a run  ------------
void 
TrackerPlots::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackerPlots::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackerPlots::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackerPlots::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackerPlots::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackerPlots);
