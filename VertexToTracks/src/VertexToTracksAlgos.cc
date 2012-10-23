#include "MGeisler/VertexToTracks/interface/VertexToTracksAlgos.h"

#include <vector>
#include <string>
#include <sstream> 
#include <algorithm>

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"


#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"

#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

const double eMass = 0.000511;
const double kMass = 0.49765;
const double lamMass = 1.11568;
const double piMass = 0.1396;


/*************************************************************************************/
/* dedicated constructor for the algorithms                                          */ 
/*************************************************************************************/

VertexToTracksAlgos::VertexToTracksAlgos(const edm::ParameterSet& iConfig)
  : maxNumWarnings_(3),
    numWarnings_(0)
{

  	cleanedColls_ = iConfig.getParameter<bool>("GetCleanedCollections");
  
  	ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  	KshortCollection_= iConfig.getParameter<InputTag>("V0KshortCollection");
  	LambdaCollection_= iConfig.getParameter<InputTag>("V0LambdaCollection");

  	NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  	input_BeamSpot_= iConfig.getParameter<InputTag>("BeamSpot");

  	ignoremissingpfcollection_ = iConfig.getParameter<bool>("ignoreMissingCollection");

}

/*************************************************************************************/
/* get all needed collections at the beginning                                       */ 
/*************************************************************************************/

bool 
VertexToTracksAlgos::GetInputCollections(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  	//get the offline beam spot
  	iEvent.getByLabel(input_BeamSpot_, beamspotH);

  	//get the conversion collection for the gamma conversions
  	iEvent.getByLabel(ConversionsCollection_, convCollH);
	cleanedConvCollP = VertexToTracksAlgos::GetCleanedConversions(convCollH,beamspotH,cleanedColls_);

  	//get the vertex composite candidate collection for the Kshort's
  	iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
	cleanedKshortCollP = VertexToTracksAlgos::GetCleanedKshort(vertCompCandCollKshortH,beamspotH,cleanedColls_);
  
  	//get the vertex composite candidate collection for the Lambda's
  	iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
	cleanedLambdaCollP = VertexToTracksAlgos::GetCleanedLambda(vertCompCandCollLambdaH,beamspotH,cleanedColls_);
  
  	//get the displaced vertex collection for nuclear interactions
  	missingColls = false;
  	if(!iEvent.getByLabel(NIVertexCollection_,displVertexCollH)){
          if (ignoremissingpfcollection_){

    	    missingColls = true; 

            if ( numWarnings_ < maxNumWarnings_ ) {
	      edm::LogWarning("VertexToTracksAlgos::GetInputCollections")
	        << "No Extra objects available in input file --> skipping reconstruction of photon conversions && displaced vertices !!" << std::endl;
	      ++numWarnings_;
            }

  	  } else {

	    return false;

	  }
	} else {

	  cleanedNICollP = VertexToTracksAlgos::GetCleanedNI(displVertexCollH,beamspotH,true);

	}

     	iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);

	return true;

}

/*************************************************************************************/
/* create helping vertex vector to remove associated vertices                        */ 
/*************************************************************************************/

std::vector<reco::VertexRef>* 
VertexToTracksAlgos::CreateVertexVector(edm::Handle<reco::VertexCollection> vtxcollH)
{

	vector<VertexRef>* output = new vector<VertexRef>();
	
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollH->size(); ++index_vtx){

          VertexRef vertexref(vtxcollH,index_vtx);

	  output->push_back(vertexref);

	}

	return output;

}

/****************************************************************************/
/* erase one vertex from the vertex vector                                  */ 
/****************************************************************************/

void 
VertexToTracksAlgos::EraseVertex(std::vector<reco::VertexRef>* vtxcollV, reco::VertexRef toErase)
{
	
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollV->size(); ++index_vtx){

          VertexRef vertexref = vtxcollV->at(index_vtx);

	  if ( vertexref == toErase ){
            vtxcollV->erase(vtxcollV->begin()+index_vtx);
	    break;
	  }

	}

}

/****************************************************************************************************************/
/* function to find the vertex with the highest TrackWeight for a certain track with the according quality      */ 
/****************************************************************************************************************/

TrackVertexQuality 
VertexToTracksAlgos::TrackWeightStep(const TrackRef&  trackRef, std::vector<reco::VertexRef>* vtxcollV, const EventSetup& iSetup) 
{

	TrackVertexQuality tw_assoc = VertexToTracksAlgos::TrackWeightAssociation(trackRef,vtxcollV);

	if ( tw_assoc.second.second == -1. ) return make_pair(trackRef, make_pair(tw_assoc.second.first, 0.));

	TransientTrack transtrk(trackRef, &(*bFieldH) );
    	transtrk.setBeamSpot(*beamspotH);
    	transtrk.setES(iSetup);

	float distance = (IPTools::absoluteImpactParameter3D(transtrk, *(tw_assoc.second.first))).second.value();

	return make_pair(trackRef, make_pair(tw_assoc.second.first, -1.*distance));

}

/*************************************************************************************/
/* function to check for the compatibility to a secondary vertex                     */ 
/*************************************************************************************/

TrackVertexQuality 
VertexToTracksAlgos::SecondaryVertexCheck(const TrackRef& trackRef, std::vector<reco::VertexRef>* vtxcollV, const EventSetup& iSetup) 
{
	
	if ( !missingColls ){

      	  // Test if the track comes from a photon conversion:
      	  // If so, try to find the vertex of the mother particle
	  Conversion gamma;
          if ( VertexToTracksAlgos::ComesFromConversion(trackRef, *cleanedConvCollP, &gamma) ){
  	    return VertexToTracksAlgos::FindConversionVertex(trackRef, gamma, bFieldH, iSetup, beamspotH, vtxcollV);
          }

      	  // Test if the track comes from a Kshort or Lambda decay:
      	  // If so, reassociate the track to the vertex of the V0
	  VertexCompositeCandidate V0;
	  if ( VertexToTracksAlgos::ComesFromV0Decay(trackRef, *cleanedKshortCollP, *cleanedLambdaCollP, &V0) ) {
            return VertexToTracksAlgos::FindV0Vertex(trackRef, V0, bFieldH, iSetup, beamspotH, vtxcollV);	
	  }

      	  // Test if the track comes from a nuclear interaction:
      	  // If so, reassociate the track to the vertex of the incoming particle 
	  PFDisplacedVertex displVtx;
	  if ( VertexToTracksAlgos::ComesFromNI(trackRef, *cleanedNICollP, &displVtx) ){
	    return VertexToTracksAlgos::FindNIVertex(trackRef, displVtx, bFieldH, iSetup, beamspotH, vtxcollV);
	  }

	}

	return make_pair(trackRef, make_pair(vtxcollV->at(0), 0.));

}

/*******************************************************************************/
/* function to fill up the remaining associations                              */ 
/*******************************************************************************/

TrackVertexQuality 
VertexToTracksAlgos::FinalAssociation(const TrackRef& trackRef, std::vector<reco::VertexRef>* vtxcollV, const EventSetup& iSetup) 
{

	TransientTrack transtrk(trackRef, &(*bFieldH) );
    	transtrk.setBeamSpot(*beamspotH);
    	transtrk.setES(iSetup);

	VertexRef closestVertex = VertexToTracksAlgos::FindClosest3D(transtrk, vtxcollV);
 
	float distance = (IPTools::absoluteImpactParameter3D(transtrk, *closestVertex)).second.value();
	if ( distance > 1000. ) distance = 1000.;

     	return make_pair(trackRef, make_pair(closestVertex, distance));

}

/****member functions ****/

/*************************************************************************************/
/* function to find the vertex with the highest TrackWeight for a certain track      */ 
/*************************************************************************************/

TrackVertexQuality 
VertexToTracksAlgos::TrackWeightAssociation(const TrackRef& trackRef, std::vector<reco::VertexRef>* vtxcollV) 
{

	VertexRef bestvertexref = vtxcollV->at(0);		
 	float bestweight = 0.;

	const TrackBaseRef& trackbaseRef = TrackBaseRef(trackRef);

	//loop over all vertices in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollV->size(); ++index_vtx){

          VertexRef vertexref = vtxcollV->at(index_vtx);

     	  //get the most probable vertex for the track
	  float weight = vertexref->trackWeight(trackbaseRef);
	  if(weight>bestweight){
  	    bestweight = weight;
	    bestvertexref = vertexref;
 	  } 

	}

	if ( bestweight>1.e-5 ){ 
	  //found a vertex with a track weight
	  //return weight == 0., meaning that a vertex could be found
  	  return make_pair(trackRef, make_pair(bestvertexref, 0.));
	} else { 
	  //found no vertex with a track weight
	  //return weight == -1., meaning that no vertex could be found
  	  return make_pair(trackRef, make_pair(bestvertexref, -1.));
	}

}

/****************************************************************************************/
/* function to calculate the deltaR between a vector and a vector connecting two points */ 
/****************************************************************************************/

double
VertexToTracksAlgos::dR(math::XYZPoint vtx_pos, math::XYZVector vtx_mom, edm::Handle<reco::BeamSpot> bsH)
{

	double bs_x = bsH->x0();
	double bs_y = bsH->y0();
	double bs_z = bsH->z0();

     	double connVec_x = vtx_pos.x() - bs_x;
	double connVec_y = vtx_pos.y() - bs_y;
	double connVec_z = vtx_pos.z() - bs_z;

     	double connVec_r = sqrt(connVec_x*connVec_x + connVec_y*connVec_y + connVec_z*connVec_z);
	double connVec_theta = acos(connVec_z*1./connVec_r);

	double connVec_eta = -1.*log(tan(connVec_theta*1./2.));
	double connVec_phi = atan2(connVec_y,connVec_x);

	return deltaR(vtx_mom.eta(),vtx_mom.phi(),connVec_eta,connVec_phi);
    
}


/*************************************************************************************/
/* function to find the closest vertex in 3D for a certain track                     */ 
/*************************************************************************************/

VertexRef 
VertexToTracksAlgos::FindClosest3D(TransientTrack transtrk, std::vector<reco::VertexRef>* vtxcollV)
{

	VertexRef foundVertexRef = vtxcollV->at(0);

	double d3min = 1e5;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<vtxcollV->size(); ++index_vtx){

          VertexRef vertexref = vtxcollV->at(index_vtx);

          double distance = 1e5;	        
          pair<bool,Measurement1D> IpPair = IPTools::absoluteImpactParameter3D(transtrk, *vertexref);
 
	  if(IpPair.first)
            distance = IpPair.second.value();	

          if(distance<d3min) {
            d3min = distance; 
            foundVertexRef = vertexref;
          }
	
	}

	return foundVertexRef;
}


/*************************************************************************************/
/* function to filter the conversion collection                                      */ 
/*************************************************************************************/

auto_ptr<ConversionCollection> 
VertexToTracksAlgos::GetCleanedConversions(edm::Handle<reco::ConversionCollection> convCollH, Handle<BeamSpot> bsH, bool cleanedColl)
{
     	auto_ptr<ConversionCollection> cleanedConvColl(new ConversionCollection() );

	for (unsigned int convcoll_idx=0; convcoll_idx<convCollH->size(); convcoll_idx++){

	  ConversionRef convref(convCollH,convcoll_idx);

 	  if(!cleanedColl){   
            cleanedConvColl->push_back(*convref);
	    continue;
          }

	  if( (convref->nTracks()==2) &&
              (fabs(convref->pairInvariantMass())<=0.1) ){
    
            cleanedConvColl->push_back(*convref);

	  }

	}

  	return cleanedConvColl;

}


/*************************************************************************************/
/* function to find out if the track comes from a gamma conversion                   */ 
/*************************************************************************************/

bool
VertexToTracksAlgos::ComesFromConversion(const TrackRef trackref, ConversionCollection cleanedConvColl, Conversion* gamma)
{

	for(unsigned int convcoll_ite=0; convcoll_ite<cleanedConvColl.size(); convcoll_ite++){

	  if(ConversionTools::matchesConversion(trackref,cleanedConvColl.at(convcoll_ite))){
	
	    *gamma = cleanedConvColl.at(convcoll_ite);
	    return true;

  	  }

  	}

	return false;
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a conversion            */ 
/*************************************************************************************/

TrackVertexQuality
VertexToTracksAlgos::FindConversionVertex(const reco::TrackRef trackref, reco::Conversion gamma, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, edm::Handle<reco::BeamSpot> bsH, std::vector<reco::VertexRef>* vtxcollV)
{ 

	math::XYZPoint conv_pos = gamma.conversionVertex().position();

	math::XYZVector conv_mom(gamma.refittedPair4Momentum().x(),
	                         gamma.refittedPair4Momentum().y(),
	                         gamma.refittedPair4Momentum().z());

	Track photon(trackref->chi2(), trackref->ndof(), conv_pos, conv_mom, 0, trackref->covariance());

    	TransientTrack transpho(photon, &(*bFieldH) );
    	transpho.setBeamSpot(*bsH);
    	transpho.setES(iSetup);

	VertexRef foundVertexRef = FindClosest3D(transpho, vtxcollV);

	TransientTrack transtrk(trackref, &(*bFieldH) );
    	transtrk.setBeamSpot(*bsH);
    	transtrk.setES(iSetup);

	float distance = (IPTools::absoluteImpactParameter3D(transtrk, *foundVertexRef)).second.value();
	if ( distance > 1000. ) distance = 1000.;
	float quality = distance + 1000.;

	return make_pair(trackref, make_pair(foundVertexRef, quality));	

}


/*************************************************************************************/
/* function to filter the Kshort collection                                          */ 
/*************************************************************************************/

auto_ptr<VertexCompositeCandidateCollection>
VertexToTracksAlgos::GetCleanedKshort(Handle<VertexCompositeCandidateCollection> KshortsH, Handle<BeamSpot> bsH, bool cleanedColl)
{

     	auto_ptr<VertexCompositeCandidateCollection> cleanedKaonColl(new VertexCompositeCandidateCollection() );

	for (unsigned int kscoll_idx=0; kscoll_idx<KshortsH->size(); kscoll_idx++){

	  VertexCompositeCandidateRef ksref(KshortsH,kscoll_idx);

 	  if(!cleanedColl){   
            cleanedKaonColl->push_back(*ksref);
	    continue;
	  }

  	  VertexDistance3D distanceComputer;

          GlobalPoint dec_pos = RecoVertex::convertPos(ksref->vertex());    

       	  GlobalError decayVertexError = GlobalError(ksref->vertexCovariance(0,0), ksref->vertexCovariance(0,1), ksref->vertexCovariance(1,1), ksref->vertexCovariance(0,2), ksref->vertexCovariance(1,2), ksref->vertexCovariance(2,2));
	
      	  math::XYZVector dec_mom(ksref->momentum().x(),
	                          ksref->momentum().y(),
	                          ksref->momentum().z());    

      	  GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
      	  GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());
      
   	  double kaon_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(dec_pos, decayVertexError))).significance();

	  if ((ksref->vertex().rho()>=3.) &&
              (ksref->vertexNormalizedChi2()<=3.) &&
              (fabs(ksref->mass() - kMass)<=0.01) &&
              (kaon_significance>15.) &&
              (VertexToTracksAlgos::dR(ksref->vertex(),dec_mom,bsH)<=0.3) ){
  
            cleanedKaonColl->push_back(*ksref);

       	  }

	}

	return cleanedKaonColl;

}


/*************************************************************************************/
/* function to filter the Lambda collection                                          */ 
/*************************************************************************************/

auto_ptr<VertexCompositeCandidateCollection>
VertexToTracksAlgos::GetCleanedLambda(Handle<VertexCompositeCandidateCollection> LambdasH, Handle<BeamSpot> bsH, bool cleanedColl)
{

     	auto_ptr<VertexCompositeCandidateCollection> cleanedLambdaColl(new VertexCompositeCandidateCollection() );

	for (unsigned int lambdacoll_idx=0; lambdacoll_idx<LambdasH->size(); lambdacoll_idx++){

	  VertexCompositeCandidateRef lambdaref(LambdasH,lambdacoll_idx);

 	  if(!cleanedColl){   
            cleanedLambdaColl->push_back(*lambdaref);
	    continue;
          }

  	  VertexDistance3D distanceComputer;

          GlobalPoint dec_pos = RecoVertex::convertPos(lambdaref->vertex());    

       	  GlobalError decayVertexError = GlobalError(lambdaref->vertexCovariance(0,0), lambdaref->vertexCovariance(0,1), lambdaref->vertexCovariance(1,1), lambdaref->vertexCovariance(0,2), lambdaref->vertexCovariance(1,2), lambdaref->vertexCovariance(2,2));
	
      	  math::XYZVector dec_mom(lambdaref->momentum().x(),
	                          lambdaref->momentum().y(),
	                          lambdaref->momentum().z());    

      	  GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
      	  GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());
      
   	  double lambda_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(dec_pos, decayVertexError))).significance();

	  if ((lambdaref->vertex().rho()>=3.) &&
              (lambdaref->vertexNormalizedChi2()<=3.) &&
              (fabs(lambdaref->mass() - lamMass)<=0.005) &&
              (lambda_significance>15.) &&
              (VertexToTracksAlgos::dR(lambdaref->vertex(),dec_mom,bsH)<=0.3) ){
  
            cleanedLambdaColl->push_back(*lambdaref);

	  }

	}

	return cleanedLambdaColl;

}

/*************************************************************************************/
/* function to find out if the track comes from a V0 decay                           */ 
/*************************************************************************************/

bool
VertexToTracksAlgos::ComesFromV0Decay(const TrackRef trackref, VertexCompositeCandidateCollection cleanedKshort, VertexCompositeCandidateCollection cleanedLambda, VertexCompositeCandidate* V0)
{

	//the part for the reassociation of particles from Kshort decays
	for(VertexCompositeCandidateCollection::const_iterator iKS=cleanedKshort.begin(); iKS!=cleanedKshort.end(); iKS++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iKS->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iKS; 
	    return true;

	  }

	}

	//the part for the reassociation of particles from Lambda decays
	for(VertexCompositeCandidateCollection::const_iterator iLambda=cleanedLambda.begin(); iLambda!=cleanedLambda.end(); iLambda++){

	  const RecoChargedCandidate *dauCand1 = dynamic_cast<const RecoChargedCandidate*>(iLambda->daughter(0));
 	  TrackRef dauTk1 = dauCand1->track();
	  const RecoChargedCandidate *dauCand2 = dynamic_cast<const RecoChargedCandidate*>(iLambda->daughter(1));
 	  TrackRef dauTk2 = dauCand2->track();

   	  if((trackref==dauTk1) || (trackref==dauTk2)){
	  
	    *V0 = *iLambda; 
	    return true;

	  }

	}

	return false;
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a V0                    */ 
/*************************************************************************************/

TrackVertexQuality
VertexToTracksAlgos::FindV0Vertex(const TrackRef trackref, VertexCompositeCandidate V0_vtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, Handle<BeamSpot> bsH, std::vector<reco::VertexRef>* vtxcollV)
{ 

	math::XYZPoint dec_pos = V0_vtx.vertex();

	math::XYZVector dec_mom(V0_vtx.momentum().x(),
	                        V0_vtx.momentum().y(),
	                        V0_vtx.momentum().z());

	Track V0(trackref->chi2(), trackref->ndof(), dec_pos, dec_mom, 0, trackref->covariance());

    	TransientTrack transV0(V0, &(*bFieldH) );
    	transV0.setBeamSpot(*bsH);
    	transV0.setES(iSetup);

	VertexRef foundVertexRef = FindClosest3D(transV0, vtxcollV); 

	TransientTrack transtrk(trackref, &(*bFieldH) );
    	transtrk.setBeamSpot(*bsH);
    	transtrk.setES(iSetup);
 
	float distance = (IPTools::absoluteImpactParameter3D(transtrk, *foundVertexRef)).second.value();
	if ( distance > 1000. ) distance = 1000.;
	float quality = distance + 1000.;

	return make_pair(trackref, make_pair(foundVertexRef, quality));		

}


/*************************************************************************************/
/* function to filter the nuclear interaction collection                             */ 
/*************************************************************************************/

auto_ptr<PFDisplacedVertexCollection>
VertexToTracksAlgos::GetCleanedNI(Handle<PFDisplacedVertexCollection> NuclIntH, Handle<BeamSpot> bsH, bool cleanedColl)
{

     	auto_ptr<PFDisplacedVertexCollection> cleanedNIColl(new PFDisplacedVertexCollection() );

	for (PFDisplacedVertexCollection::const_iterator niref=NuclIntH->begin(); niref!=NuclIntH->end(); niref++){


	  if( (niref->isFake()) || !(niref->isNucl()) ) continue;

	  if(!cleanedColl){
	    cleanedNIColl->push_back(*niref);
	    continue;
          }

  	  VertexDistance3D distanceComputer;

      	  GlobalPoint ni_pos = RecoVertex::convertPos(niref->position());    
      	  GlobalError interactionVertexError = RecoVertex::convertError(niref->error());

      	  math::XYZVector ni_mom(niref->primaryMomentum().x(),
	                         niref->primaryMomentum().y(),
	                         niref->primaryMomentum().z());

      	  GlobalPoint bsPosition = RecoVertex::convertPos(bsH->position());
      	  GlobalError bsError = RecoVertex::convertError(bsH->covariance3D());
      
   	  double nuclint_significance = (distanceComputer.distance(VertexState(bsPosition,bsError), VertexState(ni_pos, interactionVertexError))).significance();

	  if ((niref->position().rho()>=3.) &&
              (nuclint_significance>15.) &&
              (VertexToTracksAlgos::dR(niref->position(),ni_mom,bsH)<=0.3) ){
  
            cleanedNIColl->push_back(*niref);

	  }

	}

	return cleanedNIColl;

}


/*************************************************************************************/
/* function to find out if the track comes from a nuclear interaction                */ 
/*************************************************************************************/

bool
VertexToTracksAlgos::ComesFromNI(const TrackRef trackref, PFDisplacedVertexCollection cleanedNI, PFDisplacedVertex* displVtx)
{

	//the part for the reassociation of particles from nuclear interactions
	for(PFDisplacedVertexCollection::const_iterator iDisplV=cleanedNI.begin(); iDisplV!=cleanedNI.end(); iDisplV++){

	  if(iDisplV->trackWeight(trackref)>1.e-5){
	  
	    *displVtx = *iDisplV; 
	    return true;

	  }

	}

	return false;
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a nuclear interaction   */ 
/*************************************************************************************/

TrackVertexQuality
VertexToTracksAlgos::FindNIVertex(const TrackRef trackref, PFDisplacedVertex displVtx, ESHandle<MagneticField> bFieldH, const EventSetup& iSetup, Handle<BeamSpot> bsH, std::vector<reco::VertexRef>* vtxcollV)
{

	TransientTrack transtrk(trackref, &(*bFieldH) );
    	transtrk.setBeamSpot(*bsH);
    	transtrk.setES(iSetup);

	TrackCollection refittedTracks = displVtx.refittedTracks();

	if((displVtx.isTherePrimaryTracks()) || (displVtx.isThereMergedTracks())){

	  for(TrackCollection::const_iterator trkcoll_ite=refittedTracks.begin(); trkcoll_ite!=refittedTracks.end(); trkcoll_ite++){
	
	    const TrackBaseRef retrackbaseref = displVtx.originalTrack(*trkcoll_ite); 

	    if(displVtx.isIncomingTrack(retrackbaseref)){

              TrackVertexQuality VOAssociation = VertexToTracksAlgos::TrackWeightAssociation(retrackbaseref.castTo<TrackRef>(), vtxcollV);

	      if(VOAssociation.second.second==0.){
 
	        float distance = (IPTools::absoluteImpactParameter3D(transtrk, *(VOAssociation.second.first))).second.value();
	        if ( distance > 1000. ) distance = 1000.;
	        float quality = distance + 1000.;

                return make_pair(trackref, make_pair(VOAssociation.second.first, quality));

	      }

    	      TransientTrack transIncom(*retrackbaseref, &(*bFieldH) );
    	      transIncom.setBeamSpot(*bsH);
    	      transIncom.setES(iSetup);

	      VertexRef foundVertexRef = FindClosest3D(transIncom, vtxcollV); 
 
	      float distance = (IPTools::absoluteImpactParameter3D(transtrk, *foundVertexRef)).second.value();	
	      if ( distance > 1000. ) distance = 1000.;
	      float quality = distance + 1000.;

	      return make_pair(trackref, make_pair(foundVertexRef, quality));

	    }

	  }

	}

	math::XYZPoint ni_pos = displVtx.position();

	math::XYZVector ni_mom(displVtx.primaryMomentum().x(),
	                       displVtx.primaryMomentum().y(),
	                       displVtx.primaryMomentum().z());

	Track incom(trackref->chi2(), trackref->ndof(), ni_pos, ni_mom, 0, trackref->covariance());

    	TransientTrack transIncom(incom, &(*bFieldH) );
    	transIncom.setBeamSpot(*bsH);
    	transIncom.setES(iSetup);

	VertexRef foundVertexRef = FindClosest3D(transIncom, vtxcollV); 
 
	float distance = (IPTools::absoluteImpactParameter3D(transtrk, *foundVertexRef)).second.value();
	if ( distance > 1000. ) distance = 1000.;
	float quality = distance + 1000.;

     	return make_pair(trackref, make_pair(foundVertexRef, quality));

}