#include "MGeisler/TrackValidator/interface/TrackValidator.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackValidator::TrackValidator(const edm::ParameterSet& iConfig):TrackValidatorAlgos(iConfig)
{
   //now do what ever initialization is needed

  tcRefLabel_ = iConfig.getParameter<InputTag>("tcRefLabel");
  tcLabels_ = iConfig.getParameter<vector<InputTag> >("tcLabel");

  pfLabels_ = iConfig.getParameter<vector<InputTag> >("pfLabel");

  tpLabel_ = iConfig.getParameter<InputTag>("TPLabel");

  puLabel_ = iConfig.getParameter<InputTag>("PULabel");

  Service<TFileService> tfs;
  vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  for(unsigned tcl=0; tcl<tcLabels_.size(); tcl++){

    TrackValidatorAlgos::initialize();

    InputTag tColl = tcLabels_[tcl];
    string dirName = "";
    dirName += tColl.label();

    subDir->push_back(tfs->mkdir(dirName));
    TrackValidatorAlgos::BookHistos(subDir->at(tcl));

  }

  for(unsigned pfl=0; pfl<pfLabels_.size(); pfl++){

    TrackValidatorAlgos::initializePF();

    InputTag pfColl = pfLabels_[pfl];
    string dirName = "";
    dirName += pfColl.label();

    subDir->push_back(tfs->mkdir(dirName));
    TrackValidatorAlgos::BookHistosPF(subDir->at(pfl+tcLabels_.size()));

  }


  ignoremissingtkcollection_ = iConfig.getParameter<bool>("ignoremissingtrackcollection");

  photonPtMin_ = iConfig.getParameter<double>("photonPtMin");
  photonEtaMin_ = iConfig.getParameter<double>("photonEtaMin");
  photonEtaMax_ = iConfig.getParameter<double>("photonEtaMax");
  photonLip_ = iConfig.getParameter<double>("photonLip");
  photonTip_ = iConfig.getParameter<double>("photonTip");

}

TrackValidator::~TrackValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called when starting to processes a run  ------------
void 
TrackValidator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called for each event  ------------
void
TrackValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{       

  //associate reco tracks to tracking particles
  ESHandle<TrackAssociatorBase> theAssociator;
  iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
  TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product(); 

  //get the tracking particles   
  Handle<TrackingParticleCollection>  TPCollectionH ;
  iEvent.getByLabel(tpLabel_,TPCollectionH);
  const TrackingParticleCollection tPC = *(TPCollectionH.product());

  //get reference track collection from the event
  Handle<View<Track> >  RefTrackCollectionH;
  iEvent.getByLabel(tcRefLabel_,RefTrackCollectionH);

  RecoToSimCollection recSimCollRef;
  recSimCollRef=theTrackAssociator_->associateRecoToSim(RefTrackCollectionH,TPCollectionH,&iEvent,&iSetup);

  //get the pileup information  
  Handle< vector<PileupSummaryInfo> > puinfoH;
  iEvent.getByLabel(puLabel_,puinfoH);
  PileupSummaryInfo puinfo;      
  
  for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      puinfo=(*puinfoH)[puinfo_ite];
      break;
    }
  }

  int npu = puinfo.getPU_NumInteractions();

  // ########################################################
  // part of the charged particle analysis
  // ########################################################

  //loop over input collections
  for(unsigned tcl=0; tcl<tcLabels_.size(); tcl++){ 

    //get track collection from the event
    Handle<View<Track> >  trackCollectionH;
    if(!iEvent.getByLabel(tcLabels_[tcl],trackCollectionH)&&ignoremissingtkcollection_) continue;

    RecoToSimCollection recSimColl;
    recSimColl=theTrackAssociator_->associateRecoToSim(trackCollectionH,TPCollectionH,&iEvent,&iSetup);

    SimToRecoCollection simRecColl;
    simRecColl=theTrackAssociator_->associateSimToReco(trackCollectionH,TPCollectionH,&iEvent,&iSetup);

    // ########################################################
    // fill simulation histograms (LOOP OVER TRACKINGPARTICLES)
    // ########################################################

    for (TrackingParticleCollection::size_type tp_ite=0; tp_ite<TPCollectionH->size(); tp_ite++){

      TrackingParticleRef tpr(TPCollectionH, tp_ite);
      TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

      TrackValidatorAlgos::fill_generic_simTrack_histos(tcl,tp);

      const Track* matchedTrackPointer=0;
      vector<pair<RefToBase<Track>, double> > rt;

      if(simRecColl.find(tpr) != simRecColl.end()){
	rt = (vector<pair<RefToBase<Track>, double> >) simRecColl[tpr];
	if (rt.size()!=0) {
	  matchedTrackPointer = rt.begin()->first.get();
        }
      }

      TrackValidatorAlgos::fill_recoAssociated_simTrack_histos(tcl,tp,matchedTrackPointer,npu);


    }

    // #####################################################
    // fill reconstruction histograms (LOOP OVER RECOTRACKS)
    // #####################################################

    for(View<Track>::size_type rt_ite=0; rt_ite<trackCollectionH->size(); ++rt_ite){

      RefToBase<Track> track(trackCollectionH, rt_ite);

      vector<pair<TrackingParticleRef, double> > tp;
      if(recSimColl.find(track) != recSimColl.end()) tp = recSimColl[track];

      bool isMatched = false;
      bool isSigMatched = false;

      if (tp.size()!=0) {

        isMatched = true;
        for (unsigned int tp_ite=0;tp_ite<tp.size();++tp_ite){ 

          TrackingParticle trackpart = *(tp[tp_ite].first);

          if ((trackpart.eventId().event() == 0) && (trackpart.eventId().bunchCrossing() == 0)){
            isSigMatched = true;
            break;
          }

        }

      }

      TrackValidatorAlgos::fill_simAssociated_recoTrack_histos(tcl,*track,isMatched,isSigMatched,npu);

    }

    // #####################################################
    // fill pileup related histograms (LOOP OVER REF TRACKS)
    // #####################################################

    for(View<Track>::size_type ref_ite=0; ref_ite<RefTrackCollectionH->size(); ++ref_ite){

      RefToBase<Track> refTrack(RefTrackCollectionH, ref_ite);

      vector<pair<TrackingParticleRef, double> > refTp;
      if(recSimCollRef.find(refTrack) != recSimCollRef.end()) refTp = recSimCollRef[refTrack];

      bool removedTrack = true;

      for(View<Track>::size_type rt_ite=0; rt_ite<trackCollectionH->size(); ++rt_ite){

        RefToBase<Track> track(trackCollectionH, rt_ite);
	if(TrackValidatorAlgos::findRefTrack(*refTrack,*track)){
          removedTrack=false;
          break;
        }

      } 

      bool isMatched = false;
      bool isSigMatched = false;

      if (refTp.size()!=0) {

        isMatched = true;
        for (unsigned int tp_ite=0;tp_ite<refTp.size();++tp_ite){ 

          TrackingParticle trackpart = *(refTp[tp_ite].first);

          if ((trackpart.eventId().event() == 0) && (trackpart.eventId().bunchCrossing() == 0)){
            isSigMatched = true;
            break;
          }

        }

      }    

      if(isMatched) TrackValidatorAlgos::fill_removedRecoTrack_histos(tcl,*refTrack,isSigMatched,removedTrack,npu);

    }

    // ###########################
    // fill independent histograms
    // ###########################

    TrackValidatorAlgos::fill_independent_histos(tcl,npu,trackCollectionH->size()); 

  }

  // ########################################################
  // part of the uncharged particle analysis
  // ########################################################

  //loop over input collections
  for(unsigned pfl=0; pfl<pfLabels_.size(); pfl++){

    //get particle flow collection from the event
    Handle<PFCandidateCollection>  pfcCollectionH;
    if(!iEvent.getByLabel(pfLabels_[pfl],pfcCollectionH)&&ignoremissingtkcollection_) continue;

    auto_ptr<PFCandidateCollection> photons(new PFCandidateCollection() );
   
    for(unsigned pfc_ite=0;pfc_ite<pfcCollectionH->size();pfc_ite++) {
     
      PFCandidatePtr candptr(pfcCollectionH,pfc_ite);
      if((candptr->particleId()==PFCandidate::gamma) &&
         (candptr->pt()>=photonPtMin_) &&
         (candptr->momentum().eta()>=photonEtaMin_) &&
         (candptr->momentum().eta()<=photonEtaMax_) &&
         (fabs(candptr->vertex().z())<=photonTip_) &&
         (sqrt(candptr->vertex().perp2())<=photonLip_)) photons->push_back(*candptr);

    }

    TrackValidatorAlgos::fill_photon_related_histos(pfl,TPCollectionH,photons,npu);

  }

}

// ------------ method called when ending the processing of a run  ------------
void 
TrackValidator::endRun(edm::Run const&, edm::EventSetup const&)
{

  //loop over input track collections
  for(unsigned tcl=0; tcl<tcLabels_.size(); tcl++){

    TrackValidatorAlgos::fillFractionHistosFromVectors(tcl);
    TrackValidatorAlgos::fillHistosFromVectors(tcl); 
 
  }

  //loop over input pf collections
  for(unsigned pfl=0; pfl<pfLabels_.size(); pfl++){

    TrackValidatorAlgos::fillFractionHistosFromVectorsPF(pfl);
    TrackValidatorAlgos::fillHistosFromVectorsPF(pfl); 
 
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackValidator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackValidator);
