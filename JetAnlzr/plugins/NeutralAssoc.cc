// -*- C++ -*-
//
// Package:    NeutralAssoc
// Class:      NeutralAssoc
// 
/**\class NeutralAssoc NeutralAssoc.cc MGeisler/NeutralAssoc/src/NeutralAssoc.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Thu Jul 26 11:45:12 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
//
// class declaration
//

using namespace edm;
using namespace std;
using namespace reco;

class NeutralAssoc : public edm::EDAnalyzer {
   public:
      explicit NeutralAssoc(const edm::ParameterSet&);
      ~NeutralAssoc();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static VertexRef FindNeutralVertex(PFCandidateRef, Handle<VertexCollection>, Handle<BeamSpot>, double);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      InputTag pfLabel_;
      InputTag tpLabel_;

      InputTag vtxCollLabel_;
      InputTag beamSpotLabel_;

      InputTag puLabel_;
 

      TH1F* distance_histo;
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
NeutralAssoc::NeutralAssoc(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  pfLabel_ = iConfig.getParameter<InputTag>("PFCandidates");
  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  vtxCollLabel_ = iConfig.getParameter<InputTag>("VertexCollection");
  beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");

  puLabel_ = iConfig.getParameter<InputTag>("PileUp");

  //--------------


  Service<TFileService> tfs;

  distance_histo = tfs->make<TH1F>("distance_histo", "distance; distance / cm; # events", 1000, 0, 100);

}


NeutralAssoc::~NeutralAssoc()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
VertexRef
NeutralAssoc::FindNeutralVertex(PFCandidateRef candref, Handle<VertexCollection> vtxCollH, Handle<BeamSpot> bsH, double trackWeight)
{

 	math::XYZPoint xyz_in = candref->vertex();
        double theta_in = candref->theta();

	double x_new = xyz_in.x() - bsH->x(xyz_in.z());
	double y_new = xyz_in.y() - bsH->y(xyz_in.z());	

	double r_new = sqrt( x_new*x_new + y_new*y_new );

	double bs_dxdz2 = bsH->dxdz() * bsH->dxdz();
	double bs_dydz2 = bsH->dydz() * bsH->dydz();
	double bs_drdz = sqrt(bs_dxdz2 + bs_dydz2);

	double bs_theta = atan(bs_drdz);

	double theta_new = theta_in - bs_theta;

	double ztrack =  xyz_in.z() - (r_new/tan(theta_new));

	VertexRef foundVertexRef(vtxCollH,0);

	double dzmin = 5.;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxCollH->size(); ++index_vtx){

          VertexRef vertexref(vtxCollH,index_vtx);

	  int nTracks = sqrt(vertexref->tracksSize());
 
	  //find and store the closest vertex in z
          double distance = fabs(ztrack - vertexref->z());

	  double weightedDistance = max(0.0,distance-trackWeight*nTracks);	

          if(weightedDistance<dzmin) {
            dzmin = weightedDistance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;

}

// ------------ method called for each event  ------------
void
NeutralAssoc::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  	//get the pileup information  
  	Handle< vector<PileupSummaryInfo> > puinfoH;
  	iEvent.getByLabel(puLabel_,puinfoH);

  	//get the tracking vertices   
  	Handle<TrackingVertexCollection> theTVsH ;
  	iEvent.getByLabel(tpLabel_,theTVsH);

  	vector< vector<TrackingVertexRef> > realTrackingV;
  	vector< vector<bool> > realTrackingVAvail;
  
  	for (unsigned int puinfo_ite=0;puinfo_ite<puinfoH->size();++puinfo_ite){ 
  
    	  vector<TrackingVertexRef> help;
	
    	  for(int coll_ite=0; coll_ite<=puinfoH->at(puinfo_ite).getPU_NumInteractions(); coll_ite++){

      	    TrackingVertexRef tv_help;
      	    help.push_back(tv_help);

    	  }

    	  realTrackingV.push_back(help);

  	}

  	int minBC = puinfoH->at(0).getBunchCrossing();
  
  	for(TrackingVertexCollection::size_type tv_ite=0; tv_ite<theTVsH->size(); tv_ite++){

    	  TrackingVertexRef tvr(theTVsH, tv_ite);

    	  if( (tvr->eventId().event()<(signed int)realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).size()) &&
              (realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()).isNull()) &&
              (tvr->nSourceTracks()==0) ){

      	    realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()) = tvr;

    	  }

  	}

	/****************/
  
	//get the input pfCandidateCollection
  	Handle<PFCandidateCollection> theCandsH;
  	iEvent.getByLabel(pfLabel_,theCandsH);

  	//get the tracking particles   
   	Handle<TrackingParticleCollection>  theTPsH;
  	iEvent.getByLabel(tpLabel_,theTPsH);

  	//get the vertex collection  
   	Handle<VertexCollection>  vtxCollH;
  	iEvent.getByLabel(vtxCollLabel_,vtxCollH);

  	//get the offline beam spot
   	Handle<BeamSpot>  bsH;
  	iEvent.getByLabel(beamSpotLabel_,bsH);

	for (unsigned pfc_idx=0; pfc_idx<theCandsH->size(); pfc_idx++){
     
          PFCandidateRef candref(theCandsH,pfc_idx);

	  if ( (candref->trackRef().isNull()) && (candref->charge()==0) ){
          
	    double z_reco = (NeutralAssoc::FindNeutralVertex(candref,vtxCollH,bsH,0.01))->position().z();	   
  
  	    for(TrackingParticleCollection::size_type tp_ite=0; tp_ite<theTPsH->size(); tp_ite++){

      	      TrackingParticleRef tpr(theTPsH, tp_ite);

	      if( (tpr->charge()==0) &&
	          (deltaR(tpr->eta(),tpr->phi(),candref->eta(),candref->phi())<0.1) &&
                  (fabs(1- (tpr->p()*1./candref->p()))<=0.1) ){

	        double z_gen =  (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().z();

	        distance_histo->Fill(fabs(z_reco-z_gen));	       

	        break;

	      }

	    }

	  }

	}

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
NeutralAssoc::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NeutralAssoc);
