// -*- C++ -*-
//
// Package:    ReassocTest
// Class:      ReassocTest
// 
/**\class ReassocTest ReassocTest.cc MGeisler/ReassocTest/src/ReassocTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Mon Jul 23 09:31:25 CEST 2012
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CommonTools/RecoUtils/interface/PF_PU_AssoMapAlgos.h"

#include "DataFormats/Math/interface/Point3D.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// ROOT include files
#include <TH1F.h>
#include <TH2F.h>

//
// class declaration
//

class ReassocTest : public edm::EDAnalyzer, protected PF_PU_AssoMapAlgos {
   public:
      explicit ReassocTest(const edm::ParameterSet&);
      ~ReassocTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      static double GetZwrtBeamSpot(edm::Handle<reco::BeamSpot>, const math::XYZPoint, double);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      // ----------member data ---------------------------

      edm::InputTag trackLabel_;

      edm::InputTag tpLabel_;

      edm::InputTag puLabel_;

      TrackingParticleSelector* TpSelector;


      edm::InputTag input_VertexCollection_;
      edm::Handle<reco::VertexCollection> vtxcollH;

      bool input_VertexAssOneDim_;
      bool input_VertexAssClosest_;
      bool input_VertexAssUseAbsDistance_;

      bool UseBeamSpotCompatibility_;

      bool ignoremissingpfcollection_;

      double input_PtCut_;
      double input_BSCut_;
      double input_nTrack_;

      bool cleanedColls_;

      edm::InputTag input_BeamSpot_;
      edm::Handle<reco::BeamSpot> beamspotH;

      edm::InputTag ConversionsCollection_;
      edm::Handle<reco::ConversionCollection> convCollH;
      std::auto_ptr<reco::ConversionCollection> cleanedConvCollP;

      edm::InputTag KshortCollection_;
      edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollKshortH;
      std::auto_ptr<reco::VertexCompositeCandidateCollection> cleanedKshortP;

      edm::InputTag LambdaCollection_;
      edm::Handle<reco::VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
      std::auto_ptr<reco::VertexCompositeCandidateCollection> cleanedLambdaP;

      edm::InputTag NIVertexCollection_;
      edm::Handle<reco::PFDisplacedVertexCollection> displVertexCollH;
      std::auto_ptr<reco::PFDisplacedVertexCollection> cleanedNIP;
 

      TH1F* conv_distance_histo;
      TH1F* dec_distance_histo;
      TH1F* ni_distance_histo;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
ReassocTest::ReassocTest(const edm::ParameterSet& iConfig): PF_PU_AssoMapAlgos(iConfig)
{
   //now do what ever initialization is needed

  trackLabel_ = iConfig.getParameter<InputTag>("TrackCollection");

  tpLabel_ = iConfig.getParameter<InputTag>("TrackingParticles");

  puLabel_ = iConfig.getParameter<InputTag>("PULabel");

  //configure TP selectors

  using namespace reco::modules;

  ParameterSet TpSelectorPSet = iConfig.getParameter<ParameterSet>("TpSelector");

  TpSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(TpSelectorPSet));

  input_VertexCollection_= iConfig.getParameter<InputTag>("VertexCollection");

  input_VertexAssOneDim_= iConfig.getUntrackedParameter<bool>("VertexAssOneDim", true);
  input_VertexAssClosest_= iConfig.getUntrackedParameter<bool>("VertexAssClosest", true);
  input_VertexAssUseAbsDistance_= iConfig.getUntrackedParameter<bool>("VertexAssUseAbsDistance", true);
  
  ConversionsCollection_= iConfig.getParameter<InputTag>("ConversionsCollection");

  KshortCollection_= iConfig.getParameter<InputTag>("V0KshortCollection");
  LambdaCollection_= iConfig.getParameter<InputTag>("V0LambdaCollection");

  NIVertexCollection_= iConfig.getParameter<InputTag>("NIVertexCollection");

  UseBeamSpotCompatibility_= iConfig.getUntrackedParameter<bool>("UseBeamSpotCompatibility", false);
  input_BeamSpot_= iConfig.getParameter<InputTag>("BeamSpot");

  ignoremissingpfcollection_ = iConfig.getParameter<bool>("ignoreMissingCollection");

  //--------------

  input_PtCut_ = iConfig.getParameter<double>("TrackPtCut");
  input_BSCut_ = iConfig.getParameter<double>("BeamSpotCompatibilityCut");
  input_nTrack_ = iConfig.getParameter<double>("nTrackWeight");
  cleanedColls_ = iConfig.getUntrackedParameter<bool>("GetCleanedCollections", true);

  //--------------


  Service<TFileService> tfs;

  conv_distance_histo = tfs->make<TH1F>("conv_distance_histo", "distance; distance / cm; # events", 1000, 0, 100);
  dec_distance_histo = tfs->make<TH1F>("dec_distance_histo", "distance; distance / cm; # events", 1000, 0, 100);
  ni_distance_histo = tfs->make<TH1F>("ni_distance_histo", "distance; distance / cm; # events", 1000, 0, 100);

}


ReassocTest::~ReassocTest()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ReassocTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product(); 

  //get the general tracks
  Handle<TrackCollection>  theTracksH ;
  Handle<View<Track> >  theTracksV ;
  iEvent.getByLabel(trackLabel_, theTracksH);
  iEvent.getByLabel(trackLabel_, theTracksV);

  //get the tracking particles   
  Handle<TrackingParticleCollection>  theTrackingPsH ;
  iEvent.getByLabel(tpLabel_,theTrackingPsH);

  //get the tracking vertices   
  Handle<TrackingVertexCollection>  theTrackingVsH ;
  iEvent.getByLabel(tpLabel_,theTrackingVsH);

  SimToRecoCollection simRecColl;
  simRecColl=theTrackAssociator_->associateSimToReco(theTracksV,theTrackingPsH,&iEvent,&iSetup);

  RecoToSimCollection recSimColl;
  recSimColl=theTrackAssociator_->associateRecoToSim(theTracksV,theTrackingPsH,&iEvent,&iSetup);
  

  //get the offline beam spot
  iEvent.getByLabel(input_BeamSpot_, beamspotH);

  //get the conversion collection for the gamma conversions
  iEvent.getByLabel(ConversionsCollection_, convCollH);
  cleanedConvCollP = PF_PU_AssoMapAlgos::GetCleanedConversions(convCollH,beamspotH,cleanedColls_);

  //get the vertex composite candidate collection for the Kshort's
  iEvent.getByLabel(KshortCollection_, vertCompCandCollKshortH);
  cleanedKshortP = PF_PU_AssoMapAlgos::GetCleanedKshort(vertCompCandCollKshortH,beamspotH,cleanedColls_);
  
  //get the vertex composite candidate collection for the Lambda's
  iEvent.getByLabel(LambdaCollection_, vertCompCandCollLambdaH);
  cleanedLambdaP = PF_PU_AssoMapAlgos::GetCleanedLambda(vertCompCandCollLambdaH,beamspotH,cleanedColls_);

  iEvent.getByLabel(NIVertexCollection_,displVertexCollH);
  cleanedNIP = PF_PU_AssoMapAlgos::GetCleanedNI(displVertexCollH,beamspotH,cleanedColls_);
	 
  //get the input vertex collection
  iEvent.getByLabel(input_VertexCollection_, vtxcollH);

  //get the pileup information  
  Handle< vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel(puLabel_,puinfoH);

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
  
  for(TrackingVertexCollection::size_type tv_ite=0; tv_ite<theTrackingVsH->size(); tv_ite++){

    TrackingVertexRef tvr(theTrackingVsH, tv_ite);

    if( (tvr->eventId().event()<(signed int)realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).size()) &&
        (realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()).isNull()) &&
        (tvr->nSourceTracks()==0) ){

      realTrackingV.at(tvr->eventId().bunchCrossing()-minBC).at(tvr->eventId().event()) = tvr;

    }

  }

  for(size_t idxTrack=0; idxTrack<theTracksH->size(); ++idxTrack){

    TrackRef trackref = TrackRef(theTracksH, idxTrack);
    TrackBaseRef trackbaseref = TrackBaseRef(trackref);

    Conversion gamma;
    if(PF_PU_AssoMapAlgos::ComesFromConversion(trackref,*cleanedConvCollP,&gamma)){

      double z_reco = (PF_PU_AssoMapAlgos::FindConversionVertex(trackref,gamma,beamspotH,vtxcollH,0.01)).first->position().z(); 

      vector<pair<TrackingParticleRef, double> > tp;
      if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];

      if (tp.size()!=0) {

        TrackingParticleRef tpr = tp.at(0).first;

        double z_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().z();
      
        conv_distance_histo->Fill(fabs(z_reco-z_gen));     
    
      }

    }

    VertexCompositeCandidate V0;
    if(PF_PU_AssoMapAlgos::ComesFromV0Decay(trackref,*cleanedKshortP,*cleanedLambdaP,&V0)){

      double z_reco = (PF_PU_AssoMapAlgos::FindV0Vertex(trackref,V0,beamspotH,vtxcollH,0.01)).first->position().z(); 

      vector<pair<TrackingParticleRef, double> > tp;
      if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];

      if (tp.size()!=0) {

        TrackingParticleRef tpr = tp.at(0).first;

        double z_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().z();
      
        dec_distance_histo->Fill(fabs(z_reco-z_gen));     
    
      }

    }

    PFDisplacedVertex displVtx;
    if(PF_PU_AssoMapAlgos::ComesFromNI(trackref,*cleanedNIP,&displVtx)){

      double z_reco = (PF_PU_AssoMapAlgos::FindNIVertex(trackref,displVtx,beamspotH,vtxcollH,0.01)).first->position().z(); 

      vector<pair<TrackingParticleRef, double> > tp;
      if(recSimColl.find(trackbaseref) != recSimColl.end()) tp = recSimColl[trackbaseref];

      if (tp.size()!=0) {

        TrackingParticleRef tpr = tp.at(0).first;

        double z_gen = (realTrackingV.at(tpr->eventId().bunchCrossing()-minBC).at(tpr->eventId().event()))->position().z();
      
        ni_distance_histo->Fill(fabs(z_reco-z_gen));     
    
      }

    }
  
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ReassocTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ReassocTest);
