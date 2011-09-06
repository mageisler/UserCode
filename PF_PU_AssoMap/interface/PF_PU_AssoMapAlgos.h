#ifndef PF_PU_AssoMapAlgos_h
#define PF_PU_AssoMapAlgos_h


/**\class PF_PU_AssoMap PF_PU_AssoMap.cc RecoParticleFlow/PF_PU_AssoMap/plugins/PF_PU_AssoMap.cc

 Description: Produces a map with association between tracks and their particular most probable vertex with a quality of this association

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

  typedef pair<TrackRef, float> TrackQualityPair;
  typedef pair<VertexRef, TrackQualityPair> VertexTrackQuality;

  typedef vector<VertexRef > VertexRefV;

class PF_PU_AssoMapAlgos{
 public:

   static float CalculateWeight(const Vertex, const TrackRef);

   static VertexRefV* QualifiedVertices(Handle<VertexCollection>, bool, double);

   static VertexTrackQuality AssociateClosestInZ(const TrackRef, VertexRefV*);
   
   static VertexTrackQuality AssociateClosest3D(const TrackRef, VertexRefV*, const edm::EventSetup&, bool);

   static auto_ptr<TrackCollection> SecondaryTracks(TrackVertexAssMap*);
   
   static VertexRef FindConversionVertex(const TrackRef, VertexRefV*);
   
   static VertexRef FindV0Vertex(const TrackRef, VertexCompositeCandidate, VertexRefV*);
   
   static VertexRef FindNIVertex(const TrackRef, PFDisplacedVertex, VertexRefV*);

   static auto_ptr<TrackVertexAssMap> SortAssociationMap(TrackVertexAssMap*); 

 protected:
  //protected functions 

 private: 


};

#endif
