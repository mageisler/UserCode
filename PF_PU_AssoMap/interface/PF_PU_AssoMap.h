#ifndef PF_PU_AssoMap_h
#define PF_PU_AssoMap_h


/**\class PF_PU_AssoMap PF_PU_AssoMap.cc RecoParticleFlow/PF_PU_AssoMap/plugins/PF_PU_AssoMap.cc

 Description: Produces a map with association between tracks and their particular most probable vertex with a quality of this association

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MGeisler/PF_PU_AssoMap/interface/PF_PU_AssoMapAlgos.h"
   
using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

class PF_PU_AssoMap : public edm::EDProducer {
   public:
      explicit PF_PU_AssoMap(const edm::ParameterSet&);
      ~PF_PU_AssoMap();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      InputTag input_VertexCollection_;
      string input_TrackCollection_;
      string input_GsfElectronCollection_;

      bool input_VertexQuality_;
      double input_VertexMinNdof_;

      bool input_VertexAssOneDim_;
      bool input_VertexAssClosest_;
      bool input_VertexAssUseAbsDistance_;

      bool input_UseGsfElectronVertex_;
      bool input_UseCtfAssVertexForGsf_;

      InputTag ConversionsCollection_;

      InputTag KshortCollection_;
      InputTag LambdaCollection_;

      InputTag NIVertexCollection_;

};


#endif
