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
// $Id: PF_PU_AssoMap.cc,v 1.10 2011/09/06 14:02:31 mgeisler Exp $
//
//

#include "MGeisler/PF_PU_AssoMap/interface/PF_PU_AssoMap.h"

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
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversion.h"
#include "DataFormats/ParticleFlowReco/interface/PFConversionFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;
  typedef AssociationMap<OneToManyWithQuality< VertexCollection, GsfElectronCollection, float> > GsfVertexAssMap;

  typedef pair<TrackRef, float> TrackQualityPair;
  typedef vector<TrackQualityPair > TrackQualityPairVector;

  typedef pair<VertexRef, TrackQualityPair> VertexTrackQuality;

  typedef vector<VertexRef > VertexRefV;

//
// constructors and destructor
//
PF_PU_AssoMap::PF_PU_AssoMap(const edm::ParameterSet& iConfig)
{ 
   //now do what ever initialization is needed

  	input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");
  	input_TrackCollection_ = iConfig.getUntrackedParameter<string>("TrackCollection","default");
  	input_GsfElectronCollection_= iConfig.getUntrackedParameter<string>("GsfElectronCollection","default");


  	input_VertexQuality_= iConfig.getUntrackedParameter<bool>("VertexQuality", true);
  	input_VertexMinNdof_= iConfig.getUntrackedParameter<double>("VertexMinNdof", 4.);


  	input_VertexAssOneDim_= iConfig.getUntrackedParameter<bool>("VertexAssOneDim", true);
  	input_VertexAssClosest_= iConfig.getUntrackedParameter<bool>("VertexAssClosest", true);
  	input_VertexAssUseAbsDistance_= iConfig.getUntrackedParameter<bool>("VertexAssUseAbsDistance", true);


  	input_UseGsfElectronVertex_= iConfig.getUntrackedParameter<bool>("UseGsfElectronVertex", true);
  	input_UseCtfAssVertexForGsf_= iConfig.getUntrackedParameter<bool>("UseCtfAssVertexForGsf", false);

  	ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  	KshortCollection_= iConfig.getParameter<InputTag>("V0KshortCollection");
  	LambdaCollection_= iConfig.getParameter<InputTag>("V0LambdaCollection");

  	NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

	produces<TrackVertexAssMap>();
	if (input_GsfElectronCollection_!="default") produces<GsfVertexAssMap>();
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

	//variables for the best vertex for the track
	float bestweight;
   	VertexRef bestvertexref;
	TrackRef trackref;

     	auto_ptr<TrackVertexAssMap> trackvertexass(new TrackVertexAssMap() );

	//get the conversion collection for the gamma conversions
	Handle<ConversionCollection> convCollH;
	iEvent.getByLabel(ConversionsCollection_, convCollH);
	const ConversionCollection& convColl = *(convCollH.product());

	//get the vertex composite candidate collection for the Kshort's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
	iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);

	//get the vertex composite candidate collection for the Lambda's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
	iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);

	//get the displaced vertex collection for nuclear interactions
	Handle<PFDisplacedVertexCollection> displVertexCollH;
	iEvent.getByLabel(NIVertexCollection_, displVertexCollH);

	//get the input vertex collection
  	Handle<VertexCollection> vtxcoll;
  	iEvent.getByLabel(input_VertexCollection_,vtxcoll);
 	VertexRefV* qualvtxcoll = PF_PU_AssoMapAlgos::QualifiedVertices(vtxcoll,input_VertexQuality_,input_VertexMinNdof_,displVertexCollH,
									vertCompCandCollKshortH,vertCompCandCollLambdaH);

	if(qualvtxcoll->size()==0){
	  if(vtxcoll->size()>0){
            const VertexRef vertexref(vtxcoll,0);
	    qualvtxcoll->push_back(vertexref);
	  } else return;
	} 

	//######################### part of the general track collection association ###############################
	
	//look if an input for trackcollection is set
 	if (input_TrackCollection_!="default"){
    
	  //get the input track collection     
  	  Handle<TrackCollection> trkcoll;
  	  iEvent.getByLabel(input_TrackCollection_,trkcoll);
	  	
	  //loop over all tracks in the track collectio	
  	  for(unsigned int index_trck = 0;  index_trck < trkcoll->size(); ++index_trck) {

	    VertexTrackQuality VtxTrkQualAss;

 	    trackref = TrackRef(trkcoll,index_trck);
	    const TrackBaseRef& trackbaseRef = TrackBaseRef(trackref);

	    //First round of association:
	    //Looking for a qualified vertex with a track weight > 0
	    VtxTrkQualAss = PF_PU_AssoMapAlgos::TrackWeightAssociation(trackbaseRef,qualvtxcoll);

	    //Test if the track comes from a gamma conversion
	    //If so reassociate the track to the vertex of the photon
	    if(PF_PU_AssoMapAlgos::ComesFromConversion(trackref,convColl)){

// 	      VertexRef help = VtxTrkQualAss.first; 
	      VtxTrkQualAss = PF_PU_AssoMapAlgos::FindConversionVertex(trackref,qualvtxcoll);  	
// 	      if (help!=VtxTrkQualAss.first) cout << "New gamma Association" << endl;

	    }

	    //Test if the track comes from a Kshort or Lambda decay
	    //If so reassociate the track to the vertex of the V0
	    VertexCompositeCandidate V0;
	    if(PF_PU_AssoMapAlgos::ComesFromV0Decay(trackref,vertCompCandCollKshortH,vertCompCandCollLambdaH,&V0)){

// 	      VertexRef help = VtxTrkQualAss.first; 
	      VtxTrkQualAss = PF_PU_AssoMapAlgos::FindV0Vertex(trackref,V0,qualvtxcoll);	
// 	      if (help!=VtxTrkQualAss.first) cout << "New Kshort Association" << endl;

	    }

	    //Test if the track comes from a nuclear interaction
	    //If so reassociate the track to the vertex of the incomming particle
	    PFDisplacedVertex displVtx;
	    if(PF_PU_AssoMapAlgos::ComesFromNI(trackref,displVertexCollH,&displVtx)){

//            VertexRef help = VtxTrkQualAss.first; 
	      VtxTrkQualAss = PF_PU_AssoMapAlgos::FindNIVertex(trackref,displVtx,qualvtxcoll); 	
// 	      if (help!=VtxTrkQualAss.first) cout << "New nuclear interaction Association" << endl;

	    }

	    //If no vertex is found with a probability greater 0.00001
	    //and no reassociation was done look for the closest vertex in z in 3d
	    //or associate the track allways to the first vertex
	    if ((VtxTrkQualAss.second.second<=0.00001) && (VtxTrkQualAss.second.second!=-1.)){

	      // if input_VertexAssOneDim_ == true association done by closest in z or allways first vertex
	      if (input_VertexAssOneDim_){ 
	        if (input_VertexAssClosest_)
  		  //associate to closest vertex in z 
                  VtxTrkQualAss = PF_PU_AssoMapAlgos::AssociateClosestInZ(trackref,qualvtxcoll);
		else{
	          //choose always the first vertex from the vertex collection & bestweight set to -1
		  VtxTrkQualAss = make_pair(qualvtxcoll->at(0),make_pair(trackref,-1.));
		}
	      }
	      else
		VtxTrkQualAss = PF_PU_AssoMapAlgos::AssociateClosest3D(trackref,qualvtxcoll,iSetup,input_VertexAssUseAbsDistance_);

	    } 
	    
	    //insert the best vertex and the pair of track and the quality of this association in the map
            trackvertexass->insert(VtxTrkQualAss.first,VtxTrkQualAss.second);

	  }
	  
	  iEvent.put( PF_PU_AssoMapAlgos::SortAssociationMap(&(*trackvertexass)) );

	}
	
	//leave if no input for gsfelectroncollection is set
 	if (input_GsfElectronCollection_=="default") return;

	//######################### part of the gsf electron collection association ###############################

  	auto_ptr<GsfVertexAssMap> gsfvertexass(new GsfVertexAssMap() );
 
	//get the input gsfelectron collection
  	Handle<GsfElectronCollection> gsfcoll;
  	iEvent.getByLabel(input_GsfElectronCollection_,gsfcoll);

	GsfElectronRef gsfelectronref;
	
	unsigned index_gsf = 0;
	  	
	//loop over all tracks in the gsfelectron collection
  	for(GsfElectronCollection::const_iterator gsf_ite= gsfcoll->begin(); gsf_ite!=gsfcoll->end(); ++gsf_ite,++index_gsf) {

	  gsfelectronref = GsfElectronRef(gsfcoll,index_gsf);
	  trackref = (*gsfelectronref).closestCtfTrackRef();

	  //if gsfelectron's closestCtfTrack is a null reference 
   	  if (trackref.isNull()){
    
	    //get the CTFtrack collection     
  	    Handle<TrackCollection> CTFtrkcoll;
  	    iEvent.getByLabel("generalTracks",CTFtrkcoll);

     	    unsigned int index_trck=0;
     	    int ibest=-1;
     	    unsigned int sharedhits_max=0;
     	    float dr_min=1000;
     	    
	    //search the general track that shares the most hits with the electron seed
     	    for(TrackCollection::const_iterator trck_ite= CTFtrkcoll->begin(); trck_ite!=CTFtrkcoll->end(); ++trck_ite,++index_trck){
       
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

      	    trackref = TrackRef(CTFtrkcoll,ibest);

   	  }//END OF if gsfelectron's closestCtfTrack is a null reference
		  
	  if ((input_UseCtfAssVertexForGsf_) && (input_TrackCollection_!="default")){

	    //loop over all vertices in the general tracks association map
            for (TrackVertexAssMap::const_iterator assomap_ite = trackvertexass->begin(); assomap_ite != trackvertexass->end(); assomap_ite++){

	      const VertexRef assomap_vertexref = assomap_ite->key; 
  	      const TrackQualityPairVector trckcoll = assomap_ite->val;

	      TrackRef GTtrackref;

	      for (unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){
		
		GTtrackref = trckcoll[trckcoll_ite].first;

  	        if (GTtrackref==trackref){
                  bestvertexref = assomap_vertexref;
		  break;
		}

	      }

	    }

	  } else {

	    double dzmin = 10000;
	    unsigned iVertex = 0;
            unsigned index_vtx = 0;

  	    double gsfelectron_z;
	    if (input_UseGsfElectronVertex_) gsfelectron_z = (*gsfelectronref).trackPositionAtVtx().z();
       	      else gsfelectron_z = (*trackref).referencePoint().z();
          
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

	  }	
          bestweight = -2.;

	  //insert the best vertex, gsfelectron and the quality of this association in the map
          gsfvertexass->insert(bestvertexref,make_pair(gsfelectronref,bestweight));
	
	}

	iEvent.put( gsfvertexass );
}

// ------------ method called once each job just before starting event loop  ------------
void 
PF_PU_AssoMap::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PF_PU_AssoMap::endJob() 
{
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