#ifndef PFCand_NoPU_WithAM_Algos_h
#define PFCand_NoPU_WithAM_Algos_h

// -*- C++ -*-
//
// Package:    PFCand_NoPU_WithAM
// Class:      PFCand_NoPU_WithAM
// 
/**\class PF_PU_AssoMap PFCand_NoPU_WithAM.cc MGeisler/PF_PU_AssoMap/plugins/PFCand_NoPU_WithAM.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Thu Dec  1 16:07:41 CET 2011
// $Id$
//
//

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

  typedef vector<VertexRef > VertexRefV;

class PFCand_NoPU_WithAM_Algos{
 public:

   //function to create a vertex collection containing only the vertices that pass the cuts
   static VertexRefV* QualifiedVertices(Handle<VertexCollection>, bool, double);

   //function to find the closest vertex in z for a certain point
   static VertexRef FindClosestInZ(double, VertexRefV*);

   //function to compare two pfcandidates
   static bool Match(const PFCandidatePtr, const RecoCandidate*);

   //function to find out if the track comes from a gamma conversion
   static bool ComesFromConversion(const PFCandidatePtr, Handle<ConversionCollection>, VertexRefV*, VertexRef*);  

   //function to find the best vertex for a muon
   static VertexRef FindPFCandVertex(const PFCandidatePtr, VertexRefV*, const edm::EventSetup&);   

   //function to find out if the track comes from a V0 decay
   static bool ComesFromV0Decay(const PFCandidatePtr, Handle<VertexCompositeCandidateCollection>, 
	 	 	  	Handle<VertexCompositeCandidateCollection>, VertexRefV*, VertexRef*); 

   //function to find out if the track comes from a nuclear interaction
   static bool ComesFromNI(const PFCandidatePtr, Handle<PFDisplacedVertexCollection>, PFDisplacedVertex*, const edm::EventSetup&);
   
   static VertexRef FindNIVertex(const PFCandidatePtr, PFDisplacedVertex, VertexRefV*, bool, const edm::EventSetup&);  

 protected:
  //protected functions 

 private: 


};

#endif
