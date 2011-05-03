// -*- C++ -*-
//
// Package:    PF_PU_AssoMap
// Class:      PF_PU_AssoMap
// 
/**\class PF_PU_AssoMap PF_PU_AssoMap.cc RecoParticleFlow/PF_PU_AssoMap/plugins/PF_PU_AssoMap.cc

 Description: Produces a map with association between tracks and their particular most probable vertex with a quality of this association

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Tue Apr  5 18:19:28 CEST 2011
// $Id$
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
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "TMath.h"


//
// class declaration
//
   
using namespace edm;
using namespace std;
using namespace reco;


class PF_PU_AssoMap : public edm::EDProducer {
   public:
      explicit PF_PU_AssoMap(const edm::ParameterSet&);
      ~PF_PU_AssoMap();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual float CalculateWeight(const reco::Vertex vertex, 
			             const reco::TrackRef trackref);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      edm::InputTag input_TrackCollection_;
      edm::InputTag input_VertexCollection_;
      edm::InputTag input_GsfElectronCollection_;
      bool input_VertexQuality_;
      double input_VertexMinNdof_;
      bool input_ClosestVertex_;
};

//
// constructors and destructor
//
PF_PU_AssoMap::PF_PU_AssoMap(const edm::ParameterSet& iConfig)
{ 
   //now do what ever initialization is needed

  	input_TrackCollection_ = iConfig.getParameter<InputTag>("TrackCollection");

  	input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");

  	input_GsfElectronCollection_= iConfig.getParameter<InputTag>("GsfElectronCollection");

  	input_VertexQuality_= iConfig.getUntrackedParameter<bool>("VertexQuality", true);

  	input_VertexMinNdof_= iConfig.getUntrackedParameter<double>("VertexMinNdof", 4.);

  	input_ClosestVertex_= iConfig.getUntrackedParameter<bool>("ClosestVertex", true);


	produces<TrackVertexAssMap>();
  
}


PF_PU_AssoMap::~PF_PU_AssoMap()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PF_PU_AssoMap::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  	auto_ptr<TrackVertexAssMap> trackvertexass(new TrackVertexAssMap() );
 
	//get the input vertex collection
  	Handle<VertexCollection> vtxcoll;
  	iEvent.getByLabel(input_VertexCollection_,vtxcoll);

	//######################### part of the track collection association ###############################
  
	//get the input track collection
  	Handle<TrackCollection> trkcoll;
  	iEvent.getByLabel(input_TrackCollection_,trkcoll);
	  
	//variables for the best vertex for the track
	float bestweight;
   	VertexRef bestvertexref;
	TrackRef trackref;
	
	unsigned index_trck = 0;
	  	
	//loop over all tracks in the track collection
  	for(TrackCollection::const_iterator trck_ite= trkcoll->begin(); trck_ite!=trkcoll->end(); ++trck_ite,++index_trck) {	

	  trackref = TrackRef(trkcoll,index_trck);
		
	  bestweight = 0.;

	  unsigned index_vtx = 0;

	  //loop over all vertices in the vertex collection
  	  for(VertexCollection::const_iterator vtx_ite= vtxcoll->begin(); vtx_ite!=vtxcoll->end(); ++vtx_ite,++index_vtx) {

            const VertexRef vertexref(vtxcoll,index_vtx);

	    //check only those vertices with a goor quality
	    if (!(input_VertexQuality_) ||  (((*vertexref).ndof()>=input_VertexMinNdof_) && !((*vertexref).isFake()))){
 
	        //get the most probable vertex for the track
	 	float weight = CalculateWeight(*vertexref,trackref);
	     	if (weight>bestweight){
		  bestweight = weight;
	   	  bestvertexref = vertexref;
		}

	    }

	  }

	  //if no vertex is found with a probability greater 0.00001
	  if (bestweight<0.00001){
	    //choose the closest vertex in z & bestweight set to -1
	    if (input_ClosestVertex_){

  	      double dzmin = 10000;
              double ztrack = (*trackref).referencePoint().z();
	      unsigned iVertex = 0;

	      index_vtx = 0;
          
	      //loop over all vertices with a good quality in the vertex collection
              for(VertexCollection::const_iterator vtx_ite= vtxcoll->begin();vtx_ite!=vtxcoll->end();++vtx_ite,++index_vtx) {

                const VertexRef vertexref(vtxcoll,index_vtx);

	        if (!(input_VertexQuality_) ||  (((*vertexref).ndof()>=input_VertexMinNdof_) && !((*vertexref).isFake()))){
 
	          //find and store the closest vertex in z
                  double dz = fabs(ztrack - (*vertexref).z());
                  if(dz<dzmin) {
                    dzmin = dz; 
                    iVertex = index_vtx;
                  }

                }
	
	      }

	      bestvertexref = VertexRef(vtxcoll,iVertex);	
              bestweight = -1.;      

	    } 
            else {
	      //choose always the first vertex & bestweight set to -1
	      bestvertexref = VertexRef(vtxcoll,0);	
              bestweight = -1.;      

	    }
	  }

	  //insert the best vertex, track and the quality of this association in the map
          trackvertexass->insert(bestvertexref,make_pair(trackref,bestweight));

	}


	//######################### part of the gsf electron collection association ###############################
 
	//get the input particle flow collection
  	Handle<GsfElectronCollection> gsfcoll;
  	iEvent.getByLabel(input_GsfElectronCollection_,gsfcoll);

  	//choose the closest vertex in z & bestweight set to -2 for all gsfelectrons
	
	unsigned index_gsf = 0;
	  	
	//loop over all tracks in the track collection
  	for(GsfElectronCollection::const_iterator gsf_ite= gsfcoll->begin(); gsf_ite!=gsfcoll->end(); ++gsf_ite,++index_gsf) {

	  const GsfElectronRef gsfelecref(gsfcoll,index_gsf);
	  trackref = (*gsfelecref).closestCtfTrackRef();

	  //if gsfelectron's closestCtfTrack is a null reference 
   	  if (trackref.isNull()){

     	    unsigned int index_trck=0;
     	    int ibest=-1;
     	    unsigned int sharedhits_max=0;
     	    float dr_min=1000;
     	    
	    //search the general track that shares the most hits with the electron seed
     	    for(TrackCollection::const_iterator trck_ite= trkcoll->begin(); trck_ite!=trkcoll->end(); ++trck_ite,++index_trck){
       
	      unsigned int sharedhits=0;
       
              float dph= fabs(trck_ite->phi()-gsf_ite->phi()); 
       	      if (dph>TMath::Pi()) dph-= TMath::TwoPi();
       	      float det=fabs(trck_ite->eta()-gsf_ite->eta());
       	      float dr =sqrt(dph*dph+det*det);  
              
	      //loop over all valid hits of the chosen general track
       	      for(trackingRecHit_iterator trackhit_ite= trck_ite->recHitsBegin();trackhit_ite!=trck_ite->recHitsEnd();++trackhit_ite){

	        if (!(*trackhit_ite)->isValid()) continue;

	         //loop over all valid hits of the electron seed 
             	for(TrajectorySeed::const_iterator gsfhit_ite= gsf_ite->gsfTrack()->extra()->seedRef()->recHits().first;gsfhit_ite!=gsf_ite->gsfTrack()->extra()->seedRef()->recHits().second;++gsfhit_ite){
           	  
		  if (!(gsfhit_ite->isValid())) continue;
           	  if((*trackhit_ite)->sharesInput(&*(gsfhit_ite),TrackingRecHit::all))  sharedhits++; 
         	
		}       
         
       	      }
       
 
       	      if ((sharedhits>sharedhits_max) || ((sharedhits==sharedhits_max)&&(dr<dr_min))){

                sharedhits_max=sharedhits;
                dr_min=dr;
                ibest=index_trck;

       	      }

    	    }

      	    trackref = TrackRef(trkcoll,ibest);

   	  }
		  
	  double dzmin = 10000;
	  double gsfelectron_z = (*gsfelecref).trackPositionAtVtx().z();
	  unsigned iVertex = 0;

          unsigned index_vtx = 0;
          
	  //loop over all vertices with a good quality in the vertex collection
          for(VertexCollection::const_iterator vtx_ite= vtxcoll->begin();vtx_ite!=vtxcoll->end();++vtx_ite,++index_vtx) {

            const VertexRef vertexref(vtxcoll,index_vtx);

	    if (!(input_VertexQuality_) ||  (((*vertexref).ndof()>=input_VertexMinNdof_) && !((*vertexref).isFake()))){
 
	      //find and store the closest vertex in z
              double dz = fabs(gsfelectron_z - (*vertexref).z());
              if(dz<dzmin) {
                dzmin = dz; 
                iVertex = index_vtx;
              }

            }
	
          }

	  bestvertexref = VertexRef(vtxcoll,iVertex);	
          bestweight = -2.;

	  //insert the best vertex, gsfelectron and the quality of this association in the map
          trackvertexass->insert(bestvertexref,make_pair(trackref,bestweight));
	
	}

	iEvent.put( trackvertexass );
}

float 
PF_PU_AssoMap::CalculateWeight(const reco::Vertex vertex, const reco::TrackRef trackref) 
{
 
	const reco::TrackBaseRef& trackbaseRef = TrackBaseRef(trackref);

  	float weight=0.;

    	typedef reco::Vertex::trackRef_iterator IT;
    
        //loop on tracks in vertices
    	for(IT iTrack=vertex.tracks_begin(); iTrack!=vertex.tracks_end(); ++iTrack) {
	 
      	  const reco::TrackBaseRef& vertexbaseRef = *iTrack;

       	  // one of the tracks in the vertex is the same as 
      	  // the track considered in the function
      	  float w = vertex.trackWeight(vertexbaseRef);
      	  if(vertexbaseRef == trackbaseRef) {
	    if (w > weight){
	      weight=w;
	    }	 	
          }
        }

  return weight;

}

// ------------ method called once each job just before starting event loop  ------------
void 
PF_PU_AssoMap::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PF_PU_AssoMap::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
PF_PU_AssoMap::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PF_PU_AssoMap::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PF_PU_AssoMap::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PF_PU_AssoMap::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PF_PU_AssoMap::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PF_PU_AssoMap);
