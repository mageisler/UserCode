#ifndef VertexToTracksAlgos_h
#define VertexToTracksAlgos_h

/**\class VertexToTracks VertexToTracksAlgos.cc MGeisler/VertexToTracks/src/VertexToTracksAlgos.cc

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

#include <string>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"

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
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

  typedef edm::AssociationMap<edm::OneToManyWithQuality< reco::TrackCollection, reco::VertexCollection, float> > VertexToTrackAssMap;

  typedef std::pair<reco::VertexRef, float> VertexQualityPair;
  typedef std::pair<reco::TrackRef, VertexQualityPair> TrackVertexQuality;

class VertexToTracksAlgos{

 public:

   //dedicated constructor for the algorithms
   VertexToTracksAlgos(const edm::ParameterSet&);

   //get all needed collections at the beginning
   bool GetInputCollections(edm::Event&, const edm::EventSetup&);

   //create helping vertex vector to remove associated vertices
   std::vector<reco::VertexRef>* CreateVertexVector(edm::Handle<reco::VertexCollection>);

   //erase one vertex from the vertex vector
   void EraseVertex(std::vector<reco::VertexRef>*, reco::VertexRef);
  
   //function to find the vertex with the highest TrackWeight for a certain track with the according quality
   TrackVertexQuality TrackWeightStep(const reco::TrackRef&, std::vector<reco::VertexRef>*, const edm::EventSetup&);
  
   //function to check for the compatibility to a secondary vertex
   TrackVertexQuality SecondaryVertexCheck(const reco::TrackRef&, std::vector<reco::VertexRef>*, const edm::EventSetup&);
  
   //function to fill up the remaining associations
   TrackVertexQuality FinalAssociation(const reco::TrackRef&, std::vector<reco::VertexRef>*, const edm::EventSetup&);

 protected:
  //protected functions 

 private: 

  // private methods for internal usage
  
   //function to find the vertex with the highest TrackWeight for a certain track
   static TrackVertexQuality TrackWeightAssociation(const reco::TrackRef&, std::vector<reco::VertexRef>*);

   //function to calculate the deltaR between a vector and a vector connecting two points
   static double dR(math::XYZPoint, math::XYZVector, edm::Handle<reco::BeamSpot>);

   //function to find the closest vertex in 3D for a certain track
   static reco::VertexRef FindClosest3D(reco::TransientTrack, std::vector<reco::VertexRef>*);

   //function to filter the conversion collection
   static std::auto_ptr<reco::ConversionCollection> GetCleanedConversions(edm::Handle<reco::ConversionCollection>, 
                                                                          edm::Handle<reco::BeamSpot>, bool);

   //function to find out if the track comes from a gamma conversion
   static bool ComesFromConversion(const reco::TrackRef, reco::ConversionCollection, reco::Conversion*);
        
   static TrackVertexQuality FindConversionVertex(const reco::TrackRef, reco::Conversion, 	
                                                  edm::ESHandle<MagneticField>, const edm::EventSetup&, 
				                  edm::Handle<reco::BeamSpot>, std::vector<reco::VertexRef>*); 


   //function to filter the Kshort collection
   static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedKshort(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);

   //function to filter the Lambda collection
   static std::auto_ptr<reco::VertexCompositeCandidateCollection> GetCleanedLambda(edm::Handle<reco::VertexCompositeCandidateCollection>, edm::Handle<reco::BeamSpot>, bool);  
    
   //function to find out if the track comes from a V0 decay
   static bool ComesFromV0Decay(const reco::TrackRef, reco::VertexCompositeCandidateCollection, 
	 	 	  	reco::VertexCompositeCandidateCollection, reco::VertexCompositeCandidate*);
   
   static TrackVertexQuality FindV0Vertex(const reco::TrackRef, reco::VertexCompositeCandidate, 
                                          edm::ESHandle<MagneticField>, const edm::EventSetup&, 
					  edm::Handle<reco::BeamSpot>, std::vector<reco::VertexRef>*);


   //function to filter the nuclear interaction collection
   static std::auto_ptr<reco::PFDisplacedVertexCollection> GetCleanedNI(edm::Handle<reco::PFDisplacedVertexCollection>, edm::Handle<reco::BeamSpot>, bool); 

   //function to find out if the track comes from a nuclear interaction
   static bool ComesFromNI(const reco::TrackRef, reco::PFDisplacedVertexCollection, reco::PFDisplacedVertex*);
   
   static TrackVertexQuality FindNIVertex(const reco::TrackRef, reco::PFDisplacedVertex, 
                                          edm::ESHandle<MagneticField>, const edm::EventSetup&, 
 	 	                          edm::Handle<reco::BeamSpot>, std::vector<reco::VertexRef>*);

  // ----------member data ---------------------------

   bool cleanedColls_;

   edm::InputTag ConversionsCollection_;
   edm::Handle<reco::ConversionCollection> convCollH;
   std::auto_ptr<reco::ConversionCollection> cleanedConvCollP;

   edm::InputTag KshortCollection_;
   edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollKshortH;
   std::auto_ptr<reco::VertexCompositeCandidateCollection> cleanedKshortCollP;

   edm::InputTag LambdaCollection_;
   edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
   std::auto_ptr<reco::VertexCompositeCandidateCollection> cleanedLambdaCollP;

   edm::InputTag NIVertexCollection_;
   edm::Handle<reco::PFDisplacedVertexCollection> displVertexCollH;
   std::auto_ptr<reco::PFDisplacedVertexCollection> cleanedNICollP;

   edm::InputTag input_BeamSpot_;
   edm::Handle<reco::BeamSpot> beamspotH;

   edm::ESHandle<MagneticField> bFieldH;

   bool ignoremissingpfcollection_;
   bool missingColls;	 

   int maxNumWarnings_;	    //  print Warning if TrackExtra objects don't exist in input file,
   int numWarnings_;        //  but only a few times
    
};

#endif
