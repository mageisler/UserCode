// -*- C++ -*-
//
// Package:    VertexToTracks
// Class:      VertexToTracks
// 
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

#include "MGeisler/VertexToTracks/interface/VertexToTracks.h"

// system include files
#include <memory>
#include <vector>
#include <string>

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
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include <DataFormats/EgammaReco/interface/SuperCluster.h>
#include <DataFormats/EgammaReco/interface/SuperClusterFwd.h>
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

typedef edm::AssociationMap<edm::OneToManyWithQuality< reco::TrackCollection, reco::VertexCollection, float> > VertexToTrackAssMap;

typedef std::pair<reco::VertexRef, float> VertexQualityPair;
typedef std::pair<reco::TrackRef, VertexQualityPair> TrackVertexQuality;

//
// constructors and destructor
//
VertexToTracks::VertexToTracks(const edm::ParameterSet& iConfig):VertexToTracksAlgos(iConfig)
{
   //register your products

  	produces<VertexToTrackAssMap>();

   //now do what ever other initialization is needed

  	input_TrackCollection_ = iConfig.getParameter<InputTag>("TrackCollection");

  	input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");

  	input_MaxNumAssociations_= iConfig.getParameter<int>("MaxNumberOfAssociations");
  
}


VertexToTracks::~VertexToTracks()
{

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
VertexToTracks::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  	auto_ptr<VertexToTrackAssMap> vertextrackass(new VertexToTrackAssMap());
    
  	//get the input track collection     
  	Handle<TrackCollection> trkcollH;
  	iEvent.getByLabel(input_TrackCollection_, trkcollH);
    
  	//get the input vertex collection     
  	Handle<VertexCollection> vtxcollH;
  	iEvent.getByLabel(input_VertexCollection_, vtxcollH);

	int num_vertices = vtxcollH->size();

  	if ( (num_vertices==0) || !VertexToTracksAlgos::GetInputCollections(iEvent,iSetup) ){

  	  iEvent.put( vertextrackass );
          return;

        }

	if ( num_vertices < input_MaxNumAssociations_) num_vertices = input_MaxNumAssociations_;	
	    
  	//loop over all tracks of the track collection	
  	for ( size_t idxTrack = 0; idxTrack < trkcollH->size(); ++idxTrack ) {

    	  TrackRef trackref = TrackRef(trkcollH, idxTrack);

	  vector<VertexRef>* vtxColl_help = CreateVertexVector(vtxcollH);

	  int num_associations = 0;

	  // step 1: TrackWeight association check
    	  TrackVertexQuality tw_assoc = VertexToTracksAlgos::TrackWeightStep(trackref, vtxColl_help, iSetup);
    	  if ( tw_assoc.second.second != 0. ){

  	    vertextrackass->insert(tw_assoc.first, tw_assoc.second);

	    VertexRef associatedVertex = tw_assoc.second.first;
	    VertexToTracksAlgos::EraseVertex(vtxColl_help, associatedVertex);

  	    num_associations++;

          }

	  if ( num_associations == input_MaxNumAssociations_ ) continue;

	  // step 2: secondary vertices check
    	  TrackVertexQuality sv_assoc = VertexToTracksAlgos::SecondaryVertexCheck(trackref, vtxColl_help, iSetup);
    	  if ( sv_assoc.second.second != 0. ){

  	    vertextrackass->insert(sv_assoc.first, sv_assoc.second);

	    VertexRef associatedVertex = sv_assoc.second.first;
	    VertexToTracksAlgos::EraseVertex(vtxColl_help, associatedVertex);

  	    num_associations++;

          }

	  // step 3: add vertices until max number of associations per track is reached
	  for (int ite = num_associations; ite<input_MaxNumAssociations_; ite++){	    

    	    TrackVertexQuality fin_assoc = VertexToTracksAlgos::FinalAssociation(trackref, vtxColl_help, iSetup);

	    VertexRef associatedVertex = fin_assoc.second.first;
	    VertexToTracksAlgos::EraseVertex(vtxColl_help, associatedVertex);

	  }

	}

	iEvent.put( vertextrackass );
 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexToTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexToTracks);
