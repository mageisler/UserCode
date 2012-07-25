// -*- C++ -*-
//
// Package:    PerfectAssociation
// Class:      PerfectAssociation
// 
/**\class PerfectAssociation PerfectAssociation.cc MGeisler/PerfectAssociation/src/PerfectAssociation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Mon Jul 16 09:26:58 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
 
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;
using namespace reco;


//
// class declaration
//

class PerfectAssociation : public edm::EDProducer {
   public:
      explicit PerfectAssociation(const edm::ParameterSet&);
      ~PerfectAssociation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static bool FindNeutralVertex(PFCandidateRef, Handle<VertexCollection>, Handle<BeamSpot>, double);

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      InputTag tcLabel_;
      InputTag pfLabel_;
      InputTag gpLabel_;
      InputTag tpLabel_;

      //**************

      InputTag vtxCollLabel_;
      InputTag beamSpotLabel_;

      int assoType_;
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
PerfectAssociation::PerfectAssociation(const edm::ParameterSet& iConfig)
{
   //register your products

  	produces<TrackCollection>();
  	produces<PFCandidateCollection>();


   //now do what ever initialization is needed

  	tcLabel_ = iConfig.getParameter<InputTag>("trackCollection");
  	pfLabel_ = iConfig.getParameter<InputTag>("PFCandidateCollection");
  	gpLabel_ = iConfig.getParameter<InputTag>("genParticles");
  	tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");


	//**************

  	beamSpotLabel_ = iConfig.getParameter<InputTag>("BeamSpot");
  	vtxCollLabel_ = iConfig.getParameter<InputTag>("VertexCollection");

  	assoType_ = iConfig.getParameter<int>("AssociationType");
  
}


PerfectAssociation::~PerfectAssociation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool
PerfectAssociation::FindNeutralVertex(PFCandidateRef candref, Handle<VertexCollection> vtxCollH, Handle<BeamSpot> bsH, double trackWeight)
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


	VertexRef firstVertexRef(vtxCollH,0);

	VertexRef foundVertexRef;

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

	if(dzmin<5.) return (foundVertexRef == firstVertexRef);
	  else return true;

}

// ------------ method called to produce the data  ------------
void
PerfectAssociation::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	auto_ptr<PFCandidateCollection> theOutputPfc(new PFCandidateCollection() );
	auto_ptr<TrackCollection> theOutputGen(new TrackCollection());
	auto_ptr<TrackCollection> theOutputTP(new TrackCollection());

  	//associate reco tracks to tracking particles
  	ESHandle<TrackAssociatorBase> theAssociator;
  	iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
  	TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product(); 
  
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

  	//associate reco tracks to tracking particles
  	Handle<TrackCollection> theTracksH;
  	iEvent.getByLabel(tcLabel_,theTracksH);
    	Handle<View<Track> >  theTracksV;
	iEvent.getByLabel(tcLabel_,theTracksV);

  	//get the tracking particles   
  	Handle<GenParticleCollection> theGensH;
  	iEvent.getByLabel(gpLabel_,theGensH);

    	RecoToSimCollection recSimColl;
    	recSimColl=theTrackAssociator_->associateRecoToSim(theTracksV,theTPsH,&iEvent,&iSetup);

	for (unsigned pfc_idx=0; pfc_idx<theCandsH->size(); pfc_idx++){
     
          PFCandidateRef candref(theCandsH,pfc_idx);

	  TrackRef trkref = candref->trackRef();

	  if ( trkref.isNull() ){

  	    switch(assoType_){
	      case 1:{
                theOutputPfc->push_back(*candref);
                break;
              }
	      case 2:{
                if(candref->charge()==0){
                  if(PerfectAssociation::FindNeutralVertex(candref,vtxCollH,bsH,0.01)){   
	            theOutputPfc->push_back(*candref);
                  }   
                }else{
                  theOutputPfc->push_back(*candref);
                }
                break;
              }
	      default:{
                break;
              }
            }
                     

          } else {

	    const TrackBaseRef& trkBref = TrackBaseRef(trkref);

	    for (unsigned gen_idx=0; gen_idx<theGensH->size(); gen_idx++){

	      GenParticleRef genref(theGensH,gen_idx);

	      if((trkref->charge() == genref->charge()) &&
	         (fabs(deltaR(trkref->eta(),trkref->phi(),genref->eta(),genref->phi()))<=0.1) &&
	         (fabs(1.- trkref->p()/genref->p())<=0.2))
	   
	        theOutputGen->push_back(*trkref);

	    }	

      	    vector<pair<TrackingParticleRef, double> > assTPs;
      	    if(recSimColl.find(trkBref) != recSimColl.end()) assTPs = recSimColl[trkBref]; 

      	    if(assTPs.size()!=0){

              for (unsigned int tp_ite=0;tp_ite<assTPs.size();++tp_ite){

                TrackingParticle trackpart = *(assTPs[tp_ite].first);

                if ((trackpart.eventId().event() == 0) && (trackpart.eventId().bunchCrossing() == 0)){
                  theOutputTP->push_back(*trkref);
                  theOutputPfc->push_back(*candref);
                  break;
                }

              }

       	    }

	  }

	}

// 	cout << "GENs\t\t\t" << theGensH->size() << endl;
//         cout << "Reco-Tracks\t\t" << theTracksH->size() << endl;
//         cout << "Assoc-Tracks I\t\t" << theOutputGen->size() << endl;
//         cout << "Assoc-Tracks II\t\t" << theOutputTP->size() << endl;

	iEvent.put( theOutputTP );
	iEvent.put( theOutputPfc );

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PerfectAssociation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PerfectAssociation);
