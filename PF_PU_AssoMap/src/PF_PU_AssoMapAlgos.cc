#include "MGeisler/PF_PU_AssoMap/interface/PF_PU_AssoMapAlgos.h"

#include <vector>
#include <string>

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
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

  typedef pair<TrackRef, float> TrackQualityPair;
  typedef vector< TrackQualityPair > TrackQualityPairVector;
  typedef pair<VertexRef, TrackQualityPair> VertexTrackQuality;

  typedef pair <VertexRef, float>  VertexPtsumPair;
  typedef vector< VertexPtsumPair > VertexPtsumVector;

  typedef vector<VertexRef > VertexRefV;


/*******************************************************************************************/
/* function to create a vertex collection containing only the vertices that pass the cuts  */ 
/*******************************************************************************************/

VertexRefV* 
PF_PU_AssoMapAlgos::QualifiedVertices(Handle<VertexCollection> input_vtxcoll, bool input_VertexQuality_, double input_VertexMinNdof_) 
{

	VertexRefV* qualvtxcoll(new VertexRefV() );

  	for(unsigned int index_vtx = 0;  index_vtx < input_vtxcoll->size(); ++index_vtx) {

          const VertexRef vertexref(input_vtxcoll,index_vtx);

	  if(!(input_VertexQuality_) ||  (((*vertexref).ndof()>=input_VertexMinNdof_) && !((*vertexref).isFake())))
	    qualvtxcoll->push_back(vertexref);

	}

	return qualvtxcoll;

}

/*************************************************************************************/
/* function to achieve the weight of the association track<->vertex from the vertex  */ 
/*************************************************************************************/

float 
PF_PU_AssoMapAlgos::CalculateWeight(const Vertex vertex, const TrackRef trackref) 
{
 
	const TrackBaseRef& trackbaseRef = TrackBaseRef(trackref);

  	float weight=0.;

    	typedef Vertex::trackRef_iterator IT;
    
        //loop on tracks in vertices
    	for(IT iTrack=vertex.tracks_begin(); iTrack!=vertex.tracks_end(); ++iTrack) {
	 
      	  const TrackBaseRef& vertexbaseRef = *iTrack;

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


/*************************************************************************************/
/* function to associate the track to the closest vertex in z                        */ 
/*************************************************************************************/

VertexTrackQuality
PF_PU_AssoMapAlgos::AssociateClosestInZ(const TrackRef trackref, VertexRefV* qualvtxcoll)
{

	VertexRef bestvertexref;

	double dzmin = 10000;
        double ztrack = (*trackref).referencePoint().z();
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx = 0;  index_vtx < qualvtxcoll->size(); ++index_vtx) {

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
 
	  //find and store the closest vertex in z
          double dz = fabs(ztrack - (*vertexref).z());
          if(dz<dzmin) {
            dzmin = dz; 
            bestvertexref = vertexref;
          }
	
	}

	return make_pair(bestvertexref,make_pair(trackref,-1.));
}


/*******************************************************************************************/
/* function to associate the track to the closest vertex in 3D, absolue distance or sigma  */ 
/*******************************************************************************************/

VertexTrackQuality
PF_PU_AssoMapAlgos::AssociateClosest3D(const TrackRef trackref, VertexRefV* qualvtxcoll, 
				       const edm::EventSetup& iSetup, bool input_VertexAssUseAbsDistance_)
{

	VertexRef bestvertexref;

	ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	TransientTrack genTrkTT(trackref, &(*bFieldHandle) );

  	double IpMin = 10000;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx = 0;  index_vtx < qualvtxcoll->size(); ++index_vtx) {

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
	        
	  double genTrk3DIp = -1.;
	  double genTrk3DIpSig = -1.;
	  pair<bool,Measurement1D> genTrk3DIpPair = IPTools::absoluteImpactParameter3D(genTrkTT, *vertexref);

	  if (genTrk3DIpPair.first){
    	    genTrk3DIp = genTrk3DIpPair.second.value();
	    genTrk3DIpSig = genTrk3DIpPair.second.significance();
	  }
 
	  //find and store the closest vertex
	  if(input_VertexAssUseAbsDistance_){
            if(genTrk3DIp<IpMin){
              IpMin = genTrk3DIp; 
              bestvertexref = vertexref;
            }
	  }else{
            if(genTrk3DIpSig<IpMin){
              IpMin = genTrk3DIp; 
              bestvertexref = vertexref;
            }
	  }

        }	

	return make_pair(bestvertexref,make_pair(trackref,-1.));
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a conversion            */ 
/*************************************************************************************/

VertexRef
PF_PU_AssoMapAlgos::FindConversionVertex(const TrackRef trackref, VertexRefV* qualvtxcoll)
{

	double dzmin = 10000;

	double ztrackfirst = trackref->innerPosition().z();
	double tracktheta = trackref->innerMomentum().theta();
	double radius = trackref->innerPosition().rho();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	unsigned iVertex = 0;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx = 0;  index_vtx < qualvtxcoll->size(); ++index_vtx) {

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
 
	  //find and store the closest vertex in z
          double dz = fabs(ztrack - (*vertexref).z());
          if(dz<dzmin) {
            dzmin = dz; 
            iVertex = index_vtx;
          }
	
	}

	return qualvtxcoll->at(iVertex);
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a V0                    */ 
/*************************************************************************************/

VertexRef
PF_PU_AssoMapAlgos::FindV0Vertex(const TrackRef trackref, VertexCompositeCandidate V0, VertexRefV* qualvtxcoll)
{

	double dzmin = 10000;

        double ztrackfirst = trackref->innerPosition().z();

	double tracktheta = V0.p4().theta();
	double radius = trackref->innerPosition().rho();

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	unsigned iVertex = 0;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx = 0;  index_vtx < qualvtxcoll->size(); ++index_vtx) {

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
 
	  //find and store the closest vertex in z
          double dz = fabs(ztrack - (*vertexref).z());
          if(dz<dzmin) {
            dzmin = dz; 
            iVertex = index_vtx;
          }
	
	}

	return qualvtxcoll->at(iVertex);
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a nuclear interaction   */ 
/*************************************************************************************/

VertexRef
PF_PU_AssoMapAlgos::FindNIVertex(const TrackRef trackref, PFDisplacedVertex displVtx, VertexRefV* qualvtxcoll)
{

	double dzmin = 10000;

	double ztrackfirst = 0.;
	double radius = 0.;
	double tracktheta = 0.;

	if ((displVtx.isTherePrimaryTracks()) || (displVtx.isThereMergedTracks())){

	  TrackCollection refittedTracks = displVtx.refittedTracks();

	  for(TrackCollection::const_iterator trkcoll_ite = refittedTracks.begin(); trkcoll_ite != refittedTracks.end(); trkcoll_ite++){
	
	    TrackBaseRef retrackbaseref = displVtx.originalTrack(*trkcoll_ite);

	    if(displVtx.isIncomingTrack(retrackbaseref)){

              ztrackfirst = retrackbaseref->innerPosition().z();
	      radius = retrackbaseref->innerPosition().rho();
	      tracktheta = retrackbaseref->innerMomentum().theta();

	      break; 

	    }

	  }

	}else{

	  math::XYZTLorentzVector mom_sec = displVtx.secondaryMomentum((string) "PI", true);

          ztrackfirst = trackref->innerPosition().z();
	  radius = trackref->innerPosition().rho();     
	  tracktheta = mom_sec.theta();

	}	

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	unsigned iVertex = 0;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx = 0;  index_vtx < qualvtxcoll->size(); ++index_vtx) {

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
 
	  //find and store the closest vertex in z
          double dz = fabs(ztrack - (*vertexref).z());
          if(dz<dzmin) {
            dzmin = dz; 
            iVertex = index_vtx;
          }
	
	}

	return qualvtxcoll->at(iVertex);
}


/*****************************************************************************************/
/* function to create a collection of tracks from secondary vertices                     */ 
/*****************************************************************************************/

auto_ptr<TrackCollection>  
PF_PU_AssoMapAlgos::SecondaryTracks(TrackVertexAssMap* trackvertexassInput) 
{

	auto_ptr<TrackCollection> secondarytrks(new TrackCollection() );

	//loop over all vertices in the association map
        for(TrackVertexAssMap::const_iterator assomap_ite = trackvertexassInput->begin(); assomap_ite != trackvertexassInput->end(); assomap_ite++){

  	  const TrackQualityPairVector trckcoll = assomap_ite->val;
 
	  TrackRef trackref;

	  for(unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){

	    trackref = trckcoll[trckcoll_ite].first;
 	    if (trackref->trackerExpectedHitsInner().numberOfLostHits()>0) secondarytrks->push_back(*trackref);

	  }

	}

	return secondarytrks;

}


/*****************************************************************************************/
/* function to sort the vertices in the AssociationMap by the sum of (pT - pT_Error)**2  */ 
/*****************************************************************************************/

auto_ptr<TrackVertexAssMap>  
PF_PU_AssoMapAlgos::SortAssociationMap(TrackVertexAssMap* trackvertexassInput) 
{
	//create a new TrackVertexAssMap for the Output which will be sorted
     	auto_ptr<TrackVertexAssMap> trackvertexassOutput(new TrackVertexAssMap() );

	//Create and fill a vector of pairs of vertex and the summed (pT-pT_Error)**2 of the tracks associated to the vertex 
	VertexPtsumVector vertexptsumvector;

	//loop over all vertices in the association map
        for(TrackVertexAssMap::const_iterator assomap_ite = trackvertexassInput->begin(); assomap_ite != trackvertexassInput->end(); assomap_ite++){

	  const VertexRef assomap_vertexref = assomap_ite->key;
  	  const TrackQualityPairVector trckcoll = assomap_ite->val;

	  float ptsum = 0;
 
	  TrackRef trackref;

	  //get the tracks associated to the vertex and calculate the manipulated pT**2
	  for(unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++){

	    trackref = trckcoll[trckcoll_ite].first;
	    ptsum+=((*trackref).pt() - (*trackref).ptError())*((*trackref).pt() - (*trackref).ptError());

	  }

	  vertexptsumvector.push_back(make_pair(assomap_vertexref,ptsum));

	}

	while (vertexptsumvector.size()!=0){

	  VertexRef vertexref_highestpT;
	  float highestpT = 0.;
	  int highestpT_index = 0;

	  for(unsigned int vtxptsumvec_ite = 0; vtxptsumvec_ite < vertexptsumvector.size(); vtxptsumvec_ite++){
 
 	    if(vertexptsumvector[vtxptsumvec_ite].second>highestpT){

	      vertexref_highestpT = vertexptsumvector[vtxptsumvec_ite].first;
	      highestpT = vertexptsumvector[vtxptsumvec_ite].second;
	      highestpT_index = vtxptsumvec_ite;
	
	    }

	  }
	  
	  //loop over all vertices in the association map
          for(TrackVertexAssMap::const_iterator assomap_ite = trackvertexassInput->begin(); assomap_ite != trackvertexassInput->end(); assomap_ite++){

	    const VertexRef assomap_vertexref = assomap_ite->key;
  	    const TrackQualityPairVector trckcoll = assomap_ite->val;

	    //if the vertex from the association map the vertex with the highest manipulated pT 
	    //insert all associated tracks in the output Association Map
	    if(assomap_vertexref==vertexref_highestpT) 
	      for(unsigned int trckcoll_ite = 0; trckcoll_ite < trckcoll.size(); trckcoll_ite++) 
	        trackvertexassOutput->insert(assomap_vertexref,trckcoll[trckcoll_ite]);
 
	  }

	  vertexptsumvector.erase(vertexptsumvector.begin()+highestpT_index);	

	}

  	return trackvertexassOutput;

}