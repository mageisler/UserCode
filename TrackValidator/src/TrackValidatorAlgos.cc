#include "MGeisler/TrackValidator/interface/TrackValidatorAlgos.h"

// system include files
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"

#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"

// ROOT include files
#include <TH1F.h>
#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

TrackValidatorAlgos::TrackValidatorAlgos(const edm::ParameterSet& iConfig)
{
  //parameters for vs_eta plots
  minEta  = -2.5;  
  maxEta  = 2.5;
  nintEta = 50;

  //parameters for vs_pt plots
  minpt  = 0.1;
  maxpt  = 100.;
  nintpt = 40;
  
  //parameters for Pileup plots
  minVertcount  = -0.5;
  maxVertcount  = 59.5;
  nintVertcount = 60;
  
  //parameters for track number plots
  minTrackcount  = -0.5;
  maxTrackcount  = 999.5;
  nintTrackcount = 1000;
  
  //parameters for p contribution for jets plots
  minContribution  = 0.;
  maxContribution  = 1.;
  nintContribution = 200;

  //configure TP selectors

  using namespace reco::modules;

  ParameterSet generalTpSignalSelectorPSet = iConfig.getParameter<ParameterSet>("generalTpSelector");

  generalTpSignalSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpSignalSelectorPSet));

  ParameterSet generalTpPUSelectorPSet = generalTpSignalSelectorPSet;
  Entry sOname("signalOnly",false,true);
  generalTpPUSelectorPSet.insert(true,"signalOnly",sOname);

  generalTpPUSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpPUSelectorPSet));

  // fix for the LogScale by Ryan
  useLogpt_ = iConfig.getParameter<bool>("UseLogPt");

  useJetWeighting_ = iConfig.getUntrackedParameter<bool>("useJetWeighting",false);
  genJetCollLabel_ = iConfig.getUntrackedParameter<string>("genJetCollLabel");
  jetCollLabel_ = iConfig.getUntrackedParameter<string>("jetCollLabel");

  if(useLogpt_){
    maxpt=log10(maxpt);
    if(minpt > 0){
      minpt=log10(minpt);
    }else{
      minpt=log10(0.1);
    }
  }
  
  //parameters for photon selection

  etaPhiDistance=0.01;
  BARL = 1.4442; 
  END_LO = 1.566;
  END_HI = 2.5;

  MCphotonPtMin_ = iConfig.getParameter<double>("photonPtMin");
  MCphotonLip_ = iConfig.getParameter<double>("photonLip");
  MCphotonTip_ = iConfig.getParameter<double>("photonTip");

}


void 
TrackValidatorAlgos::CreateIntervalVectors()
{

  // eta vectors

  double eta_step=(maxEta-minEta)/nintEta;
  etaintervals.push_back(minEta);

  for (int k=1;k<nintEta+1;k++) {

    double d=minEta+k*eta_step;
    etaintervals.push_back(d);

  }

  // pt vectors

  double pt_step=(maxpt-minpt)/nintpt;
  ptintervals.push_back(minpt);

  for (int k=1;k<nintpt+1;k++) {

    double d;
    if(useLogpt_){
      d=pow(10,minpt+k*pt_step);
    }else{
      d=minpt+k*pt_step;
    }
    ptintervals.push_back(d);

  }

  // npu vectors

  double stepVertcount=(maxVertcount-minVertcount)/nintVertcount;
  vertcountintervals.push_back(minVertcount);

  for (int k=1;k<nintVertcount+1;k++) {

    double d=minVertcount+k*stepVertcount;
    vertcountintervals.push_back(d);

  }   

}

void 
TrackValidatorAlgos::GetInputCollections(const Event& iEvent){

  //get jet collection from the event
  if(useJetWeighting_){
    iEvent.getByLabel(genJetCollLabel_,genJetCollH);
    iEvent.getByLabel(jetCollLabel_,jetCollH);
  }

}

void 
TrackValidatorAlgos::setUpVectors()
{

  // eta vectors
  vector<double> etaintervalsh;

  for (int k=1;k<nintEta+1;k++) {
    etaintervalsh.push_back(0.);
  }

  allSignalTP_eta.push_back(etaintervalsh);
  allRT_eta.push_back(etaintervalsh);
  assSignalTP_eta.push_back(etaintervalsh);
  assSignalRT_eta.push_back(etaintervalsh);
  allSigRT_eta.push_back(etaintervalsh);
  allPURT_eta.push_back(etaintervalsh);
  allAssPURT_eta.push_back(etaintervalsh);
  allRemovedRT_eta.push_back(etaintervalsh);
  removedSigRT_eta.push_back(etaintervalsh);
  removedPURT_eta.push_back(etaintervalsh);

  // pt vectors
  vector<double> ptintervalsh;

  for (int k=1;k<nintpt+1;k++) {
    ptintervalsh.push_back(0.);
  }

  allSignalTP_pt.push_back(ptintervalsh);
  allRT_pt.push_back(ptintervalsh);
  assSignalTP_pt.push_back(ptintervalsh);
  assSignalRT_pt.push_back(ptintervalsh);
  allSigRT_pt.push_back(ptintervalsh);
  allPURT_pt.push_back(ptintervalsh);
  allRemovedRT_pt.push_back(ptintervalsh);
  removedSigRT_pt.push_back(ptintervalsh);
  removedPURT_pt.push_back(ptintervalsh);

  // npu vectors
  vector<double> vertcountintervalsh;

  for (int k=1;k<nintVertcount+1;k++) {
    vertcountintervalsh.push_back(0.);
  }   

  allSignalTP_npu.push_back(vertcountintervalsh);
  allRT_npu.push_back(vertcountintervalsh);
  assSignalTP_npu.push_back(vertcountintervalsh);
  assSignalRT_npu.push_back(vertcountintervalsh);
  allPURT_npu.push_back(vertcountintervalsh);
  allSigRT_npu.push_back(vertcountintervalsh);
  allRemovedRT_npu.push_back(vertcountintervalsh);
  removedSigRT_npu.push_back(vertcountintervalsh);
  removedPURT_npu.push_back(vertcountintervalsh);


  //counting vectors

  sim_tracks.push_back(0);

}


void 
TrackValidatorAlgos::setUpVectorsPF()
{

  // eta vectors
  vector<double> etaintervalsh;

  for (int k=1;k<nintEta+1;k++) {
    etaintervalsh.push_back(0.);
  }

  allSignalPhoton_eta.push_back(etaintervalsh);
  allRecoPhoton_eta.push_back(etaintervalsh);
  assSignalPhoton_eta.push_back(etaintervalsh);
  signalRecoPhoton_eta.push_back(etaintervalsh);
  allRemovedRecoPhoton_eta.push_back(etaintervalsh);
  removedPURecoPhoton_eta.push_back(etaintervalsh);
  removedSignalRecoPhoton_eta.push_back(etaintervalsh);
  allSignalRecoPhoton_eta.push_back(etaintervalsh);
  allPURecoPhoton_eta.push_back(etaintervalsh);

  // pt vectors
  vector<double> ptintervalsh;

  for (int k=1;k<nintpt+1;k++) {
    ptintervalsh.push_back(0.);
  }

  allSignalPhoton_pt.push_back(ptintervalsh);
  allRecoPhoton_pt.push_back(ptintervalsh);
  assSignalPhoton_pt.push_back(ptintervalsh);
  signalRecoPhoton_pt.push_back(ptintervalsh);
  allRemovedRecoPhoton_pt.push_back(etaintervalsh);
  removedSignalRecoPhoton_pt.push_back(etaintervalsh);
  removedPURecoPhoton_pt.push_back(etaintervalsh);
  allSignalRecoPhoton_pt.push_back(etaintervalsh);
  allPURecoPhoton_pt.push_back(etaintervalsh);

  // npu vectors
  vector<double> vertcountintervalsh;

  for (int k=1;k<nintVertcount+1;k++) {
    vertcountintervalsh.push_back(0.);
  }   

  allSignalPhoton_npu.push_back(vertcountintervalsh);
  allRecoPhoton_npu.push_back(vertcountintervalsh);
  assSignalPhoton_npu.push_back(vertcountintervalsh);
  signalRecoPhoton_npu.push_back(vertcountintervalsh);
  allRemovedRecoPhoton_npu.push_back(etaintervalsh);
  removedSignalRecoPhoton_npu.push_back(etaintervalsh);
  removedPURecoPhoton_npu.push_back(etaintervalsh);
  allSignalRecoPhoton_npu.push_back(etaintervalsh);
  allPURecoPhoton_npu.push_back(etaintervalsh);

}

void 
TrackValidatorAlgos::BookHistos(TFileDirectory subDir) 
{

  effic_npu_Contr.push_back(subDir.make<TH2F>("effic_npu_Contr", "efficiency vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  num_simul_tracks_npu_Contr.push_back(subDir.make<TH2F>("num_simul_tracks_npu_Contr", "Number of simulated tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  num_assoc_npu_Contr.push_back(subDir.make<TH2F>("num_assoc(simToReco)_npu_Contr", "Number of associated simulated tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  fakerate_npu_Contr.push_back(subDir.make<TH2F>("fakerate_npu_Contr", "fakerate vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  fakerate_npu_Contr_help.push_back(subDir.make<TH2F>("fakerate_npu_Contr_help", "HELP vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  num_reco_tracks_npu_Contr.push_back(subDir.make<TH2F>("num_reco_tracks_npu_Contr", "Number of reconstructed tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));
  num_assoc2_npu_Contr.push_back(subDir.make<TH2F>("num_assoc(recoToSim)_npu_Contr", "Number of associated  reco tracks vs npu vs contribution", nintVertcount, minVertcount, maxVertcount, nintContribution, minContribution, maxContribution));

  weights.push_back(subDir.make<TH1F>("weights", "weights", 1000, 0., 2.));

  //Book PileUp related histograms

  PU_effic_eta.push_back(subDir.make<TH1F>("PU_effic_eta", "PU_effic vs eta", nintEta, minEta, maxEta));
  PU_effic_pt.push_back(subDir.make<TH1F>("PU_effic_pt", "PU_effic vs pt", nintpt, minpt, maxpt));
  PU_effic_npu.push_back(subDir.make<TH1F>("PU_effic_npu", "PU_effic vs npu", nintVertcount, minVertcount, maxVertcount));

  PU_fakerate_1_eta.push_back(subDir.make<TH1F>("PU_fakerate_1_eta", "PU_fakerate 1 vs eta", nintEta, minEta, maxEta));
  PU_fakerate_1_pt.push_back(subDir.make<TH1F>("PU_fakerate_1_pt", "PU_fakerate 1 vs pt", nintpt, minpt, maxpt));
  PU_fakerate_1_npu.push_back(subDir.make<TH1F>("PU_fakerate_1_npu", "PU_fakerate 1 vs npu", nintVertcount, minVertcount, maxVertcount));

  PU_fakerate_2_eta.push_back(subDir.make<TH1F>("PU_fakerate_2_eta", "PU_fakerate 2 vs eta", nintEta, minEta, maxEta));
  PU_fakerate_2_pt.push_back(subDir.make<TH1F>("PU_fakerate_2_pt", "PU_fakerate 2 vs pt", nintpt, minpt, maxpt));
  PU_fakerate_2_npu.push_back(subDir.make<TH1F>("PU_fakerate_2_npu", "PU_fakerate 2 vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book efficiency and fakerate histograms

  effic_eta.push_back(subDir.make<TH1F>("effic_eta", "effic vs eta", nintEta, minEta, maxEta));
  effic_pt.push_back(subDir.make<TH1F>("effic_pt", "effic vs pt", nintpt, minpt, maxpt));
  effic_npu.push_back(subDir.make<TH1F>("effic_npu", "effic vs npu", nintVertcount, minVertcount, maxVertcount));

  fakerate_eta.push_back(subDir.make<TH1F>("fakerate_eta", "fakerate vs eta", nintEta, minEta, maxEta));
  fakerate_pt.push_back(subDir.make<TH1F>("fakerate_pt", "fakerate vs pt", nintpt, minpt, maxpt));
  fakerate_npu.push_back(subDir.make<TH1F>("fakerate_npu", "fakerate vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book simulation related histograms

  num_simul_tracks.push_back(subDir.make<TH1F>("num_simul_tracks", "Number of simulated tracks", nintTrackcount, minTrackcount, maxTrackcount));

  num_track_simul_eta.push_back(subDir.make<TH1F>("num_track_simul_eta", "Number of simulated tracks vs eta", nintEta, minEta, maxEta));
  num_track_simul_pt.push_back(subDir.make<TH1F>("num_track_simul_pt", "Number of simulated tracks vs pt", nintpt, minpt, maxpt));
  num_track_simul_npu.push_back(subDir.make<TH1F>("num_track_simul_npu", "Number of simulated tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_simul_vertex.push_back(subDir.make<TH1F>("num_simul_vertex", "Number of simulated vertices", nintVertcount, minVertcount, maxVertcount));

  //Book reconstruction related histograms

  num_reco_tracks.push_back(subDir.make<TH1F>("num_reco_tracks", "Number of reconstructed tracks", nintTrackcount, minTrackcount, maxTrackcount));

  num_track_reco_eta.push_back(subDir.make<TH1F>("num_track_reco_eta", "Number of reconstructed tracks vs eta", nintEta, minEta, maxEta));
  num_track_reco_pt.push_back(subDir.make<TH1F>("num_track_reco_pt", "Number of reconstructed tracks vs pt", nintpt, minpt, maxpt));
  num_track_reco_npu.push_back(subDir.make<TH1F>("num_track_reco_npu", "Number of reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_removed_reco_signal_eta.push_back(subDir.make<TH1F>("num_removed_reco_signal_eta", "Number of removed reconstructed signal tracks vs eta", nintEta, minEta, maxEta));
  num_removed_reco_eta.push_back(subDir.make<TH1F>("num_removed_reco_eta", "Number of removed reconstructed tracks vs eta", nintEta, minEta, maxEta));

  num_removed_reco_PU_eta.push_back(subDir.make<TH1F>("num_removed_reco_PU_eta", "Number of removed reconstructed pileup tracks vs eta", nintEta, minEta, maxEta));
  num_reco_PU_eta.push_back(subDir.make<TH1F>("num_reco_PU_eta", "Number of reconstructed pileup tracks vs eta", nintEta, minEta, maxEta));

  num_removed_reco_signal_pt.push_back(subDir.make<TH1F>("num_removed_reco_signal_pt", "Number of removed reconstructed signal tracks vs pt", nintpt, minpt, maxpt));
  num_removed_reco_pt.push_back(subDir.make<TH1F>("num_removed_reco_pt", "Number of removed reconstructed tracks vs pt", nintpt, minpt, maxpt));
  num_removed_reco_PU_pt.push_back(subDir.make<TH1F>("num_removed_reco_PU_pt", "Number of removed reconstructed pileup tracks vs pt", nintpt, minpt, maxpt));

  num_removed_reco_signal_npu.push_back(subDir.make<TH1F>("num_removed_reco_signal_npu", "Number of removed reconstructed signal tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_removed_reco_npu.push_back(subDir.make<TH1F>("num_removed_reco_npu", "Number of removed reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_removed_reco_PU_npu.push_back(subDir.make<TH1F>("num_removed_reco_PU_npu", "Number of removed reconstructed pileup tracks vs npu", nintVertcount, minVertcount, maxVertcount));

  num_track_reco_PU_eta.push_back(subDir.make<TH1F>("num_track_reco_PU_eta", "Number of reco tracks from PileUp vs eta", nintEta, minEta, maxEta));
  num_track_reco_signal_eta.push_back(subDir.make<TH1F>("num_track_reco_signal_eta", "Number of reco tracks from Signal vs eta", nintEta, minEta, maxEta));

  num_track_reco_PU_pt.push_back(subDir.make<TH1F>("num_track_reco_PU_pt", "Number of reco tracks from PileUp vs pt", nintpt, minpt, maxpt));
  num_track_reco_signal_pt.push_back(subDir.make<TH1F>("num_track_reco_signal_pt", "Number of reco tracks from Signal vs pt", nintpt, minpt, maxpt));

  num_track_reco_PU_npu.push_back(subDir.make<TH1F>("num_track_reco_PU_npu", "Number of reco tracks from PileUp vs npu", nintVertcount, minVertcount, maxVertcount));
  num_track_reco_signal_npu.push_back(subDir.make<TH1F>("num_track_reco_signal_npu", "Number of reco tracks from Signal vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book association related histograms

  num_assoc_eta.push_back(subDir.make<TH1F>("num_assoc(simToReco)_eta", "Number of associated simulated tracks vs eta", nintEta, minEta, maxEta));
  num_assoc2_eta.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_eta", "Number of associated reconstructed tracks vs eta", nintEta, minEta, maxEta));

  num_assoc_pt.push_back(subDir.make<TH1F>("num_assoc(simToReco)_pt", "Number of associated simulated tracks vs pt", nintpt, minpt, maxpt));
  num_assoc2_pt.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_pt", "Number of associated reconstructed tracks vs pt", nintpt, minpt, maxpt));

  num_assoc_npu.push_back(subDir.make<TH1F>("num_assoc(simToReco)_npu", "Number of associated simulated tracks vs npu", nintVertcount, minVertcount, maxVertcount));
  num_assoc2_npu.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_npu", "Number of associated reconstructed tracks vs npu", nintVertcount, minVertcount, maxVertcount));

}

void 
TrackValidatorAlgos::BookHistosPF(TFileDirectory subDir) 
{

  //Book PileUp related histograms

  photon_PU_effic_eta.push_back(subDir.make<TH1F>("photon_PU_effic_eta", "Photon PU_effic vs eta", nintEta, minEta, maxEta));
  photon_PU_effic_pt.push_back(subDir.make<TH1F>("photon_PU_effic_pt", "Photon PU_effic vs pt", nintpt, minpt, maxpt));
  photon_PU_effic_npu.push_back(subDir.make<TH1F>("photon_PU_effic_npu", "Photon PU_effic vs npu", nintVertcount, minVertcount, maxVertcount));

  num_removedPURecoPhoton_eta.push_back(subDir.make<TH1F>("num_removedPURecoPhoton_eta", "Number of removed reconstructed pileup photons vs eta", nintEta, minEta, maxEta));
  num_removedPURecoPhoton_pt.push_back(subDir.make<TH1F>("num_removedPURecoPhoton_pt", "Number of removed reconstructed pileup photons vs pt", nintpt, minpt, maxpt));
  num_removedPURecoPhoton_npu.push_back(subDir.make<TH1F>("num_removedPURecoPhoton_npu", "Number of removed reconstructed pileup photons vs npu", nintVertcount, minVertcount, maxVertcount));

  num_allPURecoPhoton_eta.push_back(subDir.make<TH1F>("num_allPURecoPhoton_eta", "Number of reconstructed pileup photons vs eta", nintEta, minEta, maxEta));
  num_allPURecoPhoton_pt.push_back(subDir.make<TH1F>("num_allPURecoPhoton_pt", "Number of reconstructed pileup photons vs pt", nintpt, minpt, maxpt));
  num_allPURecoPhoton_npu.push_back(subDir.make<TH1F>("num_allPURecoPhoton_npu", "Number of reconstructed pileup photons vs npu", nintVertcount, minVertcount, maxVertcount));


  photon_PU_fakerate_1_eta.push_back(subDir.make<TH1F>("photon_PU_fakerate_1_eta", "Photon PU_fakerate 1 vs eta", nintEta, minEta, maxEta));
  photon_PU_fakerate_1_pt.push_back(subDir.make<TH1F>("photon_PU_fakerate_1_pt", "Photon PU_fakerate 1 vs pt", nintpt, minpt, maxpt));
  photon_PU_fakerate_1_npu.push_back(subDir.make<TH1F>("photon_PU_fakerate_1_npu", "Photon PU_fakerate 1 vs npu", nintVertcount, minVertcount, maxVertcount));

  photon_PU_fakerate_2_eta.push_back(subDir.make<TH1F>("photon_PU_fakerate_2_eta", "Photon PU_fakerate 2 vs eta", nintEta, minEta, maxEta));
  photon_PU_fakerate_2_pt.push_back(subDir.make<TH1F>("photon_PU_fakerate_2_pt", "Photon PU_fakerate 2 vs pt", nintpt, minpt, maxpt));
  photon_PU_fakerate_2_npu.push_back(subDir.make<TH1F>("photon_PU_fakerate_2_npu", "Photon PU_fakerate 2 vs npu", nintVertcount, minVertcount, maxVertcount));

  num_removedSignalRecoPhoton_eta.push_back(subDir.make<TH1F>("num_removedSignalRecoPhoton_eta", "Number of removed reconstructed signal photons vs eta", nintEta, minEta, maxEta));
  num_removedSignalRecoPhoton_pt.push_back(subDir.make<TH1F>("num_removedSignalRecoPhoton_pt", "Number of removed reconstructed signal photons vs pt", nintpt, minpt, maxpt));
  num_removedSignalRecoPhoton_npu.push_back(subDir.make<TH1F>("num_removedSignalRecoPhoton_npu", "Number of removed reconstructed signal photons vs npu", nintVertcount, minVertcount, maxVertcount));

  num_allRemovedRecoPhoton_eta.push_back(subDir.make<TH1F>("num_allRemovedRecoPhoton_eta", "Number of removed reconstructed photons vs eta", nintEta, minEta, maxEta));
  num_allRemovedRecoPhoton_pt.push_back(subDir.make<TH1F>("num_allRemovedRecoPhoton_pt", "Number of removed reconstructed pileup vs pt", nintpt, minpt, maxpt));
  num_allRemovedRecoPhoton_npu.push_back(subDir.make<TH1F>("num_allRemovedRecoPhoton_npu", "Number of removed reconstructed photons vs npu", nintVertcount, minVertcount, maxVertcount));

  num_allSignalRecoPhoton_eta.push_back(subDir.make<TH1F>("num_allSignalRecoPhoton_eta", "Number of reconstructed signal photons vs eta", nintEta, minEta, maxEta));
  num_allSignalRecoPhoton_pt.push_back(subDir.make<TH1F>("num_allSignalRecoPhoton_pt", "Number of reconstructed signal photons vs pt", nintpt, minpt, maxpt));
  num_allSignalRecoPhoton_npu.push_back(subDir.make<TH1F>("num_allSignalRecoPhoton_npu", "Number of reconstructed signal photons vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book simulation related histograms

  num_photon_simul.push_back(subDir.make<TH1F>("num_photon_simul", "Number of simulated photons", nintTrackcount, minTrackcount, maxTrackcount));
  num_photon_simul_eta.push_back(subDir.make<TH1F>("num_photon_simul_eta", "Number of simulated photons vs eta", nintEta, minEta, maxEta));

  num_photon_simul_pt.push_back(subDir.make<TH1F>("num_photon_simul_pt", "Number of simulated photons vs pt", nintpt, minpt, maxpt));
  num_photon_simul_npu.push_back(subDir.make<TH1F>("num_photon_simul_npu", "Number of simulated photons vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book reconstruction related histograms

  num_photon_reco.push_back(subDir.make<TH1F>("num_photon_reco", "Number of reconstructed photons", nintTrackcount, minTrackcount, maxTrackcount));
  num_photon_reco_eta.push_back(subDir.make<TH1F>("num_photon_reco_eta", "Number of reconstructed photons vs eta", nintEta, minEta, maxEta));

  num_photon_reco_pt.push_back(subDir.make<TH1F>("num_photon_reco_pt", "Number of reconstructed photons vs pt", nintpt, minpt, maxpt));
  num_photon_reco_npu.push_back(subDir.make<TH1F>("num_photon_reco_npu", "Number of reconstructed photons vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book association related histograms

  num_assoc_photon_eta.push_back(subDir.make<TH1F>("num_assoc(simToReco)_photon_eta", "Number of associated simulated photons vs eta", nintEta, minEta, maxEta));
  num_assoc2_photon_eta.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_photon_eta", "Number of associated reconstructed photons vs eta", nintEta, minEta, maxEta));

  num_assoc_photon_pt.push_back(subDir.make<TH1F>("num_assoc(simToReco)_photon_pt", "Number of associated simulated photons vs pt", nintpt, minpt, maxpt));
  num_assoc2_photon_pt.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_photon_pt", "Number of associated reconstructed photons vs pt", nintpt, minpt, maxpt));

  num_assoc_photon_npu.push_back(subDir.make<TH1F>("num_assoc(simToReco)_photon_npu", "Number of associated simulated photons vs npu", nintVertcount, minVertcount, maxVertcount));
  num_assoc2_photon_npu.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_photon_npu", "Number of associated reconstructed photons vs npu", nintVertcount, minVertcount, maxVertcount));

  //Book efficiency and fakerate histograms

  photon_effic_eta.push_back(subDir.make<TH1F>("photon_effic_eta", "photon effic vs eta", nintEta, minEta, maxEta));
  photon_effic_pt.push_back(subDir.make<TH1F>("photon_effic_pt", "photon effic vs pt", nintpt, minpt, maxpt));
  photon_effic_npu.push_back(subDir.make<TH1F>("photon_effic_npu", "photon effic vs npu", nintVertcount, minVertcount, maxVertcount));

  photon_fakerate_eta.push_back(subDir.make<TH1F>("photon_fakerate_eta", "photon fakerate vs eta", nintEta, minEta, maxEta));
  photon_fakerate_pt.push_back(subDir.make<TH1F>("photon_fakerate_pt", "photon fakerate vs pt", nintpt, minpt, maxpt));
  photon_fakerate_npu.push_back(subDir.make<TH1F>("photon_fakerate_npu", "photon fakerate vs npu", nintVertcount, minVertcount, maxVertcount));

}

void 
TrackValidatorAlgos::fill_independent_histos(int counter, int npu, int rt, unsigned st)
{

  num_simul_vertex.at(counter)->Fill(npu);
  num_simul_tracks.at(counter)->Fill(st);
  num_reco_tracks.at(counter)->Fill(rt);

}

void
TrackValidatorAlgos::fill_recoAssociated_simTrack_histos(int counter, TrackingParticle* tp, const Track* track, int npu, unsigned* st)
{

  bool isMatched = track;
  double tp_eta = tp->momentum().eta();

  if((*generalTpSignalSelector)(*tp)){

    double weight = 1.;
    if(useJetWeighting_) weight = getTrackWeight(tp,&(*genJetCollH));

    (*st)++;

    num_simul_tracks_npu_Contr[counter]->Fill(npu,weight);
    if(isMatched) num_assoc_npu_Contr[counter]->Fill(npu,weight);

    //effic vs eta
    for(unsigned int f=0; f<etaintervals.size()-1; f++){
      if(tp_eta>etaintervals[f]&&
	 tp_eta<etaintervals[f+1]){
	allSignalTP_eta[counter][f]+=weight;
	if(isMatched){
	  assSignalTP_eta[counter][f]+=weight;
	}
        break;
      }
    } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

    //effic vs pt
    for(unsigned int f=0; f<ptintervals.size()-1; f++){
      if(sqrt(tp->momentum().perp2())>ptintervals[f]&&
	 sqrt(tp->momentum().perp2())<ptintervals[f+1]){
        allSignalTP_pt[counter][f]+=weight; 
        if(isMatched){
	  assSignalTP_pt[counter][f]+=weight;
        }	
        break;      
      }
    } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){

    //effic vs num pileup vertices
    for(unsigned int f=0; f<vertcountintervals.size()-1; f++){
      if(npu == vertcountintervals[f]){
        allSignalTP_npu[counter][f]+=weight;
        if(isMatched){
          assSignalTP_npu[counter][f]+=weight;
        }
        break;
      }    
    }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){

  }


}

void 
TrackValidatorAlgos::fill_simAssociated_recoTrack_histos(int counter, const Track& track, bool isMatched, bool isSigMatched, int npu, TpDoubV tp)
{

  double weight = 1.;
  if(useJetWeighting_) weight = getTrackWeightReco(track,&(*jetCollH));

  num_reco_tracks_npu_Contr[counter]->Fill(npu,weight);
  if(isSigMatched) num_assoc2_npu_Contr[counter]->Fill(npu,weight);

  //fake rate vs eta
  for (unsigned int f=0; f<etaintervals.size()-1; f++){
    if (track.eta()>etaintervals[f]&&
        track.eta()<etaintervals[f+1]) {
      allRT_eta[counter][f]+=weight;
      if (isSigMatched){
	assSignalRT_eta[counter][f]+=weight;
      }else{
        if(isMatched) allAssPURT_eta[counter][f]+=weight;
      }
      break;
    }
  } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

  //fake rate vs pt
  for (unsigned int f=0; f<ptintervals.size()-1; f++){
    if (sqrt(track.momentum().perp2())>ptintervals[f]&&
        sqrt(track.momentum().perp2())<ptintervals[f+1]) {
      allRT_pt[counter][f]+=weight;
      if (isSigMatched){
	assSignalRT_pt[counter][f]+=weight;
      }	  
      break;    
    }
  } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){

  //fake rate vs num pileup vertices
  for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
    if (npu == vertcountintervals[f]) {
      allRT_npu[counter][f]+=weight;
      if (isSigMatched){
        assSignalRT_npu[counter][f]+=weight;
      }
      break;
    }    
  }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){


}

void 
TrackValidatorAlgos::fill_removedRecoTrack_histos(int counter, const Track& refTrack, bool isSigMatched,  bool isRemoved, int npu, TpDoubV tp)
{

  double weight = 1.;
  if(useJetWeighting_){
    weight = getTrackWeightReco(refTrack,&(*jetCollH));
    weights[counter]->Fill(weight);
  }

  // vs eta
  for (unsigned int f=0; f<etaintervals.size()-1; f++){
    if (refTrack.eta()>etaintervals[f]&&
        refTrack.eta()<etaintervals[f+1]) {
      if(isSigMatched){
        allSigRT_eta[counter][f]+=weight; 
       }else{
        allPURT_eta[counter][f]+=weight; 
      }   
     if(isRemoved){
        allRemovedRT_eta[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_eta[counter][f]+=weight;
        }else{
          removedPURT_eta[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(refTrack.eta()>etaintervals[counter][f]&&refTrack.eta()<etaintervals[counter][f+1){
  } // END for(unsigned int f=0; f<etaintervals.size()-1; f++){

  // vs pt
  for (unsigned int f=0; f<ptintervals.size()-1; f++){
    if (sqrt(refTrack.momentum().perp2())>ptintervals[f]&&
        sqrt(refTrack.momentum().perp2())<ptintervals[f+1]){
      if(isSigMatched){
        allSigRT_pt[counter][f]+=weight; 
      }else{
        allPURT_pt[counter][f]+=weight; 
      }   
      if(isRemoved){
        allRemovedRT_pt[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_pt[counter][f]+=weight; 
        }else{
	  removedPURT_pt[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(sqrt(refTrack.momentum().perp2())>ptintervals[counter][f]&&sqrt(refTrack.momentum().perp2())<ptintervals[counter][f+1]){
  } // END for(unsigned int f=0; f<ptintervals.size()-1; f++){

  // vs npu
  for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
    if (npu == vertcountintervals[f]) {
      if(isSigMatched){
        allSigRT_npu[counter][f]+=weight; 
      }else{
        allPURT_npu[counter][f]+=weight; 
      }   
      if(isRemoved){
        allRemovedRT_npu[counter][f]+=weight;
        if(isSigMatched){
          removedSigRT_npu[counter][f]+=weight;
        }else{
  	  removedPURT_npu[counter][f]+=weight; 
        }
      } // END if(isRemoved){
    } // END if(npu == vertcountintervals[counter][f]){
  } // END for(unsigned int f=0; f<vertcountintervals.size()-1; f++){
   
}

void 
TrackValidatorAlgos::fill_photon_related_histos(int counter, vector<PhotonMCTruth>  mcPhotons, auto_ptr<PFCandidateCollection> photons, auto_ptr<PFCandidateCollection> RefPhotons, SimVertex MainInt, int npu)
{

  //fill simulation related histos

  unsigned num_photons_simul = 0;

  for(vector<PhotonMCTruth>::const_iterator mcPho=mcPhotons.begin(); mcPho !=mcPhotons.end(); mcPho++){

    float mcPhi= (*mcPho).fourMomentum().phi();
    float mcEta= (*mcPho).fourMomentum().pseudoRapidity();
    mcEta = etaTransformation(mcEta, (*mcPho).primaryVertex().z() );
    float PX = (*mcPho).fourMomentum().px();
    float PY = (*mcPho).fourMomentum().py();
    float PT = sqrt(PX*PX + PY*PY);

    if((photonSelector(*mcPho,MainInt)) && (fabs(mcEta) <= BARL || (fabs(mcEta) >= END_LO && fabs(mcEta) <=END_HI))){

      bool isMatched = false;

      num_photons_simul++;

      // Loop over recontructed photons
      for(PFCandidateConstIterator iPho = photons->begin(); iPho != photons->end(); iPho++){

        if(photonMatching(mcPhi,mcEta,PT,*iPho)){

          isMatched = true;
          break;

        }

      }

      //effic vs eta
      for(unsigned int f=0; f<etaintervals.size()-1; f++){
        if(mcEta>etaintervals[f]&&
	   mcEta<etaintervals[f+1]){
	  allSignalPhoton_eta[counter][f]++;
	  if(isMatched){
	    assSignalPhoton_eta[counter][f]++;
	  }
          break;
        }
      } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

      //effic vs pt
      for(unsigned int f=0; f<ptintervals.size()-1; f++){
        if(PT>ptintervals[f]&&
	   PT<ptintervals[f+1]){
          allSignalPhoton_pt[counter][f]++; 
          if(isMatched){
	    assSignalPhoton_pt[counter][f]++;
          }	
          break;      
        }
      } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){

      //effic vs num pileup vertices
      for(unsigned int f=0; f<vertcountintervals.size()-1; f++){
        if(npu == vertcountintervals[f]){
          allSignalPhoton_npu[counter][f]++;
          if(isMatched){
            assSignalPhoton_npu[counter][f]++;
          }
          break;
        }    
      }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){

    }

  }

  num_photon_simul.at(counter)->Fill(num_photons_simul);

  //fill reconstruction related histos

  for(PFCandidateConstIterator iPho = photons->begin(); iPho != photons->end(); iPho++){

    bool isSigMatched = false;

    for(vector<PhotonMCTruth>::const_iterator mcPho=mcPhotons.begin(); mcPho !=mcPhotons.end(); mcPho++){

      float mcPhi= (*mcPho).fourMomentum().phi();
      float mcEta= (*mcPho).fourMomentum().pseudoRapidity();
      mcEta = etaTransformation(mcEta, (*mcPho).primaryVertex().z() );
      float PX = (*mcPho).fourMomentum().px();
      float PY = (*mcPho).fourMomentum().py();
      float PT = sqrt(PX*PX + PY*PY);

      if(photonMatching(mcPhi,mcEta,PT,*iPho)){

        if((mcPho->primaryVertex().x()==MainInt.position().x()) &&
           (mcPho->primaryVertex().y()==MainInt.position().y()) &&
           (mcPho->primaryVertex().z()==MainInt.position().z())) isSigMatched=true;
        break;

      }

    }

    //fake rate vs eta
    for (unsigned int f=0; f<etaintervals.size()-1; f++){
      if (iPho->eta()>etaintervals[f]&&
          iPho->eta()<etaintervals[f+1]) {
        allRecoPhoton_eta[counter][f]++;
        if(isSigMatched){
  	  signalRecoPhoton_eta[counter][f]++;
        }
        break;
      }
    } // END for (unsigned int f=0; f<etaintervals.size()-1; f++){

    //fake rate vs pt
    for (unsigned int f=0; f<ptintervals.size()-1; f++){
      if (sqrt(iPho->momentum().perp2())>ptintervals[f]&&
          sqrt(iPho->momentum().perp2())<ptintervals[f+1]) {
        allRecoPhoton_pt[counter][f]++; 
        if(isSigMatched){
  	  signalRecoPhoton_pt[counter][f]++;
        }	  
        break;    
      }
    } // End for (unsigned int f=0; f<ptintervals.size()-1; f++){
 
    //fake rate vs num pileup vertices
    for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
      if (npu == vertcountintervals[f]) {
        allRecoPhoton_npu[counter][f]++;
        if(isSigMatched){
          signalRecoPhoton_npu[counter][f]++;
        }
        break;
      }    
    }// End for (unsigned int f=0; f<vertcountintervals.size()-1; f++){

  }

  for(PFCandidateConstIterator refPho = RefPhotons->begin(); refPho != RefPhotons->end(); refPho++){

    bool isSigMatched = false;
    bool isRemoved = true;

    for(vector<PhotonMCTruth>::const_iterator mcPho=mcPhotons.begin(); mcPho !=mcPhotons.end(); mcPho++){

      float mcPhi= (*mcPho).fourMomentum().phi();
      float mcEta= (*mcPho).fourMomentum().pseudoRapidity();
      mcEta = etaTransformation(mcEta, (*mcPho).primaryVertex().z() );
      float PX = (*mcPho).fourMomentum().px();
      float PY = (*mcPho).fourMomentum().py();
      float PT = sqrt(PX*PX + PY*PY);

      if(photonMatching(mcPhi,mcEta,PT,*refPho)){

        if((mcPho->primaryVertex().x()==MainInt.position().x()) &&
           (mcPho->primaryVertex().y()==MainInt.position().y()) &&
           (mcPho->primaryVertex().z()==MainInt.position().z())) isSigMatched=true;
        break;

      }

    }

    for(PFCandidateConstIterator iPho = photons->begin(); iPho != photons->end(); iPho++){

      if(photonMatching(refPho->eta(),refPho->eta(),refPho->pt(),*iPho)){

        isRemoved=false;
        break;

      }

    }

    // vs eta
    for (unsigned int f=0; f<etaintervals.size()-1; f++){
      if (refPho->eta()>etaintervals[f]&&
          refPho->eta()<etaintervals[f+1]) {
        if(isSigMatched){
          allSignalRecoPhoton_eta[counter][f]++; 
         }else{
          allPURecoPhoton_eta[counter][f]++; 
        }    
        if(isRemoved){
          allRemovedRecoPhoton_eta[counter][f]++;
          if(isSigMatched){
            removedSignalRecoPhoton_eta[counter][f]++; 
          }else{
            removedPURecoPhoton_eta[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(refPho->eta()>etaintervals[f]&&refPho->eta()<etaintervals[f+1]){
    } // END for(unsigned int f=0; f<etaintervals.size()-1; f++){

    // vs pt
    for (unsigned int f=0; f<ptintervals.size()-1; f++){
      if (sqrt(refPho->momentum().perp2())>ptintervals[f]&&
          sqrt(refPho->momentum().perp2())<ptintervals[f+1]){
        if(isSigMatched){
          allSignalRecoPhoton_pt[counter][f]++; 
         }else{
          allPURecoPhoton_pt[counter][f]++; 
        }    
        if(isRemoved){
          allRemovedRecoPhoton_pt[counter][f]++;
          if(isSigMatched){
            removedSignalRecoPhoton_pt[counter][f]++; 
          }else{
            removedPURecoPhoton_pt[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(sqrt(refPho->momentum().perp2())>ptintervals[f]&&sqrt(refPho->momentum().perp2())<ptintervals[f+1]]){
    } // END for(unsigned int f=0; f<ptintervals.size()-1; f++){

    // vs npu
    for (unsigned int f=0; f<vertcountintervals.size()-1; f++){
      if (npu == vertcountintervals[f]) {
        if(isSigMatched){
          allSignalRecoPhoton_npu[counter][f]++; 
         }else{
          allPURecoPhoton_npu[counter][f]++; 
        }    
        if(isRemoved){
          allRemovedRecoPhoton_npu[counter][f]++;
          if(isSigMatched){
            removedSignalRecoPhoton_npu[counter][f]++; 
          }else{
            removedPURecoPhoton_npu[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(npu == vertcountintervals[counter][f]){
    } // END for(unsigned int f=0; f<vertcountintervals.size()-1; f++){

  }

  num_photon_reco.at(counter)->Fill(photons->size());

}

bool 
TrackValidatorAlgos::photonSelector(PhotonMCTruth pho, SimVertex MainInt)
{

  float PX = pho.fourMomentum().px();
  float PY = pho.fourMomentum().py();
  float PT = sqrt(PX*PX + PY*PY);

  return ((PT>=MCphotonPtMin_) &&
         (fabs(pho.vertex().z())<=MCphotonTip_) &&
         (sqrt(pho.vertex().perp2())<=MCphotonLip_) &&
         (pho.primaryVertex().x()==MainInt.position().x()) &&
         (pho.primaryVertex().y()==MainInt.position().y()) &&
         (pho.primaryVertex().z()==MainInt.position().z()));

}

bool 
TrackValidatorAlgos::photonMatching(float mcPhi, float mcEta, float mcPt, PFCandidate reco_photon)
{
 
  float recoPt = reco_photon.pt();

  float phiClu=reco_photon.phi();
  float etaClu=reco_photon.eta();
  float deltaPhi = phiClu-mcPhi;
  float deltaEta = etaClu-mcEta;

  if(deltaPhi > Geom::pi()) deltaPhi -= Geom::twoPi();
  if(deltaPhi < -Geom::pi()) deltaPhi += Geom::twoPi();

  deltaPhi=pow(deltaPhi,2);
  deltaEta=pow(deltaEta,2);
  float delta =  deltaPhi+deltaEta ;

  return ((delta < etaPhiDistance) &&
          (fabs(1.-recoPt/mcPt)<=0.1));

}


float 
TrackValidatorAlgos::etaTransformation(  float EtaParticle , float Zvertex)
{

  //Definitions
  const float PI    = 3.1415927;

  //Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479;

  //ETA correction

  float Theta = 0.0  ;
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+PI ;
  float ETA = - log(tan(0.5*Theta));

  if( fabs(ETA) > etaBarrelEndcap ){
    
    float Zend = Z_Endcap ;
    if(EtaParticle<0.0 )  Zend = -Zend ;

    float Zlen = Zend - Zvertex ;
    float RR = Zlen/sinh(EtaParticle);

    Theta = atan(RR/Zend);
    if(Theta<0.0) Theta = Theta+PI ;
    ETA = - log(tan(0.5*Theta));

  }

  return ETA;

}

void
TrackValidatorAlgos::fillFractionHistosFromVectors(int counter)
{

  fillFractionHisto(effic_eta[counter],assSignalTP_eta[counter],allSignalTP_eta[counter],"effic");
  fillFractionHisto(effic_pt[counter],assSignalTP_pt[counter],allSignalTP_pt[counter],"effic");
  fillFractionHisto(effic_npu[counter],assSignalTP_npu[counter],allSignalTP_npu[counter],"effic");

  fillFractionHisto(fakerate_eta[counter],assSignalRT_eta[counter],allRT_eta[counter],"fakerate");
  fillFractionHisto(fakerate_pt[counter],assSignalRT_pt[counter],allRT_pt[counter],"fakerate");
  fillFractionHisto(fakerate_npu[counter],assSignalRT_npu[counter],allRT_npu[counter],"fakerate");

  fillFractionHisto(PU_effic_eta[counter],removedPURT_eta[counter],allPURT_eta[counter],"effic");
  fillFractionHisto(PU_effic_pt[counter],removedPURT_pt[counter],allPURT_pt[counter],"effic");
  fillFractionHisto(PU_effic_npu[counter],removedPURT_npu[counter],allPURT_npu[counter],"effic");

  fillFractionHisto(PU_fakerate_1_eta[counter],removedSigRT_eta[counter],allRemovedRT_eta[counter],"effic");
  fillFractionHisto(PU_fakerate_1_pt[counter],removedSigRT_pt[counter],allRemovedRT_pt[counter],"effic");
  fillFractionHisto(PU_fakerate_1_npu[counter],removedSigRT_npu[counter],allRemovedRT_npu[counter],"effic");

  fillFractionHisto(PU_fakerate_2_eta[counter],removedSigRT_eta[counter],allSigRT_eta[counter],"effic");
  fillFractionHisto(PU_fakerate_2_pt[counter],removedSigRT_pt[counter],allSigRT_pt[counter],"effic");
  fillFractionHisto(PU_fakerate_2_npu[counter],removedSigRT_npu[counter],allSigRT_npu[counter],"effic");


  effic_npu_Contr[counter]->Divide(num_assoc_npu_Contr[counter],num_simul_tracks_npu_Contr[counter]);

  fakerate_npu_Contr_help[counter]->Divide(num_assoc2_npu_Contr[counter],num_reco_tracks_npu_Contr[counter]);

  for(int x_ite=1; x_ite<=nintVertcount; x_ite++){
    for(int y_ite=1; y_ite<=nintContribution; y_ite++){
      fakerate_npu_Contr[counter]->SetBinContent(x_ite,y_ite, 1. - fakerate_npu_Contr_help[counter]->GetBinContent(x_ite,y_ite));
      fakerate_npu_Contr[counter]->SetBinError(x_ite,y_ite,fakerate_npu_Contr_help[counter]->GetBinError(x_ite,y_ite));
    }
  }

}

void
TrackValidatorAlgos::fillFractionHistosFromVectorsPF(int counter)
{

  fillFractionHisto(photon_effic_eta[counter],assSignalPhoton_eta[counter],allSignalPhoton_eta[counter],"effic");
  fillFractionHisto(photon_effic_pt[counter],assSignalPhoton_pt[counter],allSignalPhoton_pt[counter],"effic");
  fillFractionHisto(photon_effic_npu[counter],assSignalPhoton_npu[counter],allSignalPhoton_npu[counter],"effic");

  fillFractionHisto(photon_fakerate_eta[counter],signalRecoPhoton_eta[counter],allRecoPhoton_eta[counter],"fakerate");
  fillFractionHisto(photon_fakerate_pt[counter],signalRecoPhoton_pt[counter],allRecoPhoton_pt[counter],"fakerate");
  fillFractionHisto(photon_fakerate_npu[counter],signalRecoPhoton_npu[counter],allRecoPhoton_npu[counter],"fakerate");

  fillFractionHisto(photon_PU_effic_eta[counter],removedPURecoPhoton_eta[counter],allPURecoPhoton_eta[counter],"effic");
  fillFractionHisto(photon_PU_effic_pt[counter],removedPURecoPhoton_pt[counter],allPURecoPhoton_pt[counter],"effic");
  fillFractionHisto(photon_PU_effic_npu[counter],removedPURecoPhoton_npu[counter],allPURecoPhoton_npu[counter],"effic");

  fillFractionHisto(photon_PU_fakerate_1_eta[counter],removedSignalRecoPhoton_eta[counter],allRemovedRecoPhoton_eta[counter],"effic");
  fillFractionHisto(photon_PU_fakerate_1_pt[counter],removedSignalRecoPhoton_pt[counter],allRemovedRecoPhoton_pt[counter],"effic");
  fillFractionHisto(photon_PU_fakerate_1_npu[counter],removedSignalRecoPhoton_npu[counter],allRemovedRecoPhoton_npu[counter],"effic");

  fillFractionHisto(photon_PU_fakerate_2_eta[counter],removedSignalRecoPhoton_eta[counter],allSignalRecoPhoton_eta[counter],"effic");
  fillFractionHisto(photon_PU_fakerate_2_pt[counter],removedSignalRecoPhoton_pt[counter],allSignalRecoPhoton_pt[counter],"effic");
  fillFractionHisto(photon_PU_fakerate_2_npu[counter],removedSignalRecoPhoton_npu[counter],allSignalRecoPhoton_npu[counter],"effic");

}

void
TrackValidatorAlgos::fillHistosFromVectors(int counter)
{

  fillPlotFromVector(num_track_simul_eta[counter],allSignalTP_eta[counter]);
  fillPlotFromVector(num_assoc_eta[counter],assSignalTP_eta[counter]);

  fillPlotFromVector(num_track_simul_pt[counter],allSignalTP_pt[counter]);
  fillPlotFromVector(num_assoc_pt[counter],assSignalTP_pt[counter]);

  fillPlotFromVector(num_track_simul_npu[counter],allSignalTP_npu[counter]);
  fillPlotFromVector(num_assoc_npu[counter],assSignalTP_npu[counter]);


  fillPlotFromVector(num_track_reco_eta[counter],allRT_eta[counter]);
  fillPlotFromVector(num_assoc2_eta[counter],assSignalRT_eta[counter]);

  fillPlotFromVector(num_track_reco_pt[counter],allRT_pt[counter]);
  fillPlotFromVector(num_assoc2_pt[counter],assSignalRT_pt[counter]);

  fillPlotFromVector(num_track_reco_npu[counter],allRT_npu[counter]);
  fillPlotFromVector(num_assoc2_npu[counter],assSignalRT_npu[counter]);


  fillPlotFromVector(num_removed_reco_signal_eta[counter],removedSigRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_eta[counter],allRemovedRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_PU_eta[counter],removedPURT_eta[counter]);

  fillPlotFromVector(num_removed_reco_signal_pt[counter],removedSigRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_pt[counter],allRemovedRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_PU_pt[counter],removedPURT_pt[counter]);

  fillPlotFromVector(num_removed_reco_signal_npu[counter],removedSigRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_npu[counter],allRemovedRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_PU_npu[counter],removedPURT_npu[counter]);


  fillPlotFromVector(num_track_reco_signal_eta[counter],allSigRT_eta[counter]);
  fillPlotFromVector(num_track_reco_PU_eta[counter],allPURT_eta[counter]);

  fillPlotFromVector(num_track_reco_signal_pt[counter],allSigRT_pt[counter]);
  fillPlotFromVector(num_track_reco_PU_pt[counter],allPURT_pt[counter]);

  fillPlotFromVector(num_track_reco_signal_npu[counter],allSigRT_npu[counter]);
  fillPlotFromVector(num_track_reco_PU_npu[counter],allPURT_npu[counter]);


  fillPlotFromVector(num_reco_PU_eta[counter],allAssPURT_eta[counter]);

}

void
TrackValidatorAlgos::fillHistosFromVectorsPF(int counter)
{

  fillPlotFromVector(num_photon_simul_eta[counter],allSignalPhoton_eta[counter]);
  fillPlotFromVector(num_assoc_photon_eta[counter],assSignalPhoton_eta[counter]);

  fillPlotFromVector(num_photon_simul_pt[counter],allSignalPhoton_pt[counter]);
  fillPlotFromVector(num_assoc_photon_pt[counter],assSignalPhoton_pt[counter]);

  fillPlotFromVector(num_photon_simul_npu[counter],allSignalPhoton_npu[counter]);
  fillPlotFromVector(num_assoc_photon_npu[counter],assSignalPhoton_npu[counter]);


  fillPlotFromVector(num_photon_reco_eta[counter],allRecoPhoton_eta[counter]);
  fillPlotFromVector(num_assoc2_photon_eta[counter],signalRecoPhoton_eta[counter]);

  fillPlotFromVector(num_photon_reco_pt[counter],allRecoPhoton_pt[counter]);
  fillPlotFromVector(num_assoc2_photon_pt[counter],signalRecoPhoton_pt[counter]);

  fillPlotFromVector(num_photon_reco_npu[counter],allRecoPhoton_npu[counter]);
  fillPlotFromVector(num_assoc2_photon_npu[counter],signalRecoPhoton_npu[counter]);


  fillPlotFromVector(num_removedPURecoPhoton_eta[counter],removedPURecoPhoton_eta[counter]);
  fillPlotFromVector(num_allPURecoPhoton_eta[counter],allPURecoPhoton_eta[counter]);

  fillPlotFromVector(num_removedPURecoPhoton_pt[counter],removedPURecoPhoton_pt[counter]);
  fillPlotFromVector(num_allPURecoPhoton_pt[counter],allPURecoPhoton_pt[counter]);

  fillPlotFromVector(num_removedPURecoPhoton_npu[counter],removedPURecoPhoton_npu[counter]);
  fillPlotFromVector(num_allPURecoPhoton_npu[counter],allPURecoPhoton_npu[counter]);


  fillPlotFromVector(num_removedSignalRecoPhoton_eta[counter],removedSignalRecoPhoton_eta[counter]);
  fillPlotFromVector(num_removedSignalRecoPhoton_pt[counter],removedSignalRecoPhoton_pt[counter]);
  fillPlotFromVector(num_removedSignalRecoPhoton_npu[counter],removedSignalRecoPhoton_npu[counter]);

  fillPlotFromVector(num_allRemovedRecoPhoton_eta[counter],allRemovedRecoPhoton_eta[counter]);
  fillPlotFromVector(num_allRemovedRecoPhoton_pt[counter],allRemovedRecoPhoton_pt[counter]);
  fillPlotFromVector(num_allRemovedRecoPhoton_npu[counter],allRemovedRecoPhoton_npu[counter]);

  fillPlotFromVector(num_allSignalRecoPhoton_eta[counter],allSignalRecoPhoton_eta[counter]);
  fillPlotFromVector(num_allSignalRecoPhoton_pt[counter],allSignalRecoPhoton_pt[counter]);
  fillPlotFromVector(num_allSignalRecoPhoton_npu[counter],allSignalRecoPhoton_npu[counter]);

}

void
TrackValidatorAlgos::fillPlotFromVector(TH1F* histo,vector<double> vIn)
{

  for(unsigned ite=0; ite<vIn.size();ite++){
    histo->SetBinContent(ite+1,vIn.at(ite));
  } 

}

void
TrackValidatorAlgos::fillFractionHisto(TH1F* histo,vector<double> num,vector<double> denum,string type)
{

  for (unsigned int j=0; j<num.size(); j++){

    double val = 0.;
    double err = 0.;

    if (denum[j]!=0){
      if (type=="effic"){
	val = ((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1-val)/(double) denum[j] );
      } else if (type=="fakerate"){
	val = 1-((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1-val)/(double) denum[j] );
      } else if (type=="pileup"){
	val = ((double) num[j])*1./((double) denum[j]);
        err = sqrt( val*(1+val)/(double) denum[j] );
      } else return;
    }

    histo->SetBinContent(j+1, val);
    histo->SetBinError(j+1, err);

  }

}

bool
TrackValidatorAlgos::findRefTrack(const Track& refTrack,const Track& track)
{

	return (
	  refTrack.eta() == track.eta() &&
	  refTrack.phi() == track.phi() &&
	  refTrack.chi2() == track.chi2() &&
	  refTrack.ndof() == track.ndof() &&
	  refTrack.p() == track.p()
	);

}

double
TrackValidatorAlgos::getTrackWeight(const TrackingParticle* tp, const GenJetCollection* genJetColl)
{


  	double trackP = tp->p();

  	for(unsigned jet_ite=0; jet_ite<genJetColl->size(); jet_ite++){

    	  GenJetRef jetRef(genJetColl,jet_ite); 

    	  double jetP = jetRef->p();

 	  vector<const GenParticle*> constituents = jetRef->getGenConstituents();

	  for (unsigned const_ite=0; const_ite<constituents.size(); const_ite++) {

	    if(isGenPart(tp,constituents[const_ite])) return trackP*1./jetP;

  	  }
	
  	}

  	return 0.0001;

}

double
TrackValidatorAlgos::getTrackWeightReco(const Track& track, const PFJetCollection* jetColl)
{

  	double trackP = track.p();

  	for(unsigned jet_ite=0; jet_ite<jetColl->size(); jet_ite++){

    	  PFJetRef jetRef(jetColl,jet_ite); 

    	  double jetP = jetRef->p();

  	  for (unsigned daugth_ite=0; daugth_ite<jetRef->numberOfDaughters(); daugth_ite++){

    	    const PFCandidatePtr pfcand = jetRef->getPFConstituent(daugth_ite);
    	    TrackRef trackref = pfcand->trackRef();
    	    if((trackref.isAvailable()) && (trackref.isNonnull())){
              if(findRefTrack(track,*trackref)) return trackP*1./jetP;
    	    }
  	  }
	
  	}

  	return 0.0001;

}

bool
TrackValidatorAlgos::isGenPart(const TrackingParticle* tp, const GenParticle* gp)
{

	return((fabs(tp->momentum().eta() - gp->eta())<0.001) &&
	       (fabs(tp->momentum().phi() - gp->phi())<0.001) &&
	       (fabs(1.-(tp->pt() / gp->pt()))<0.01));

}