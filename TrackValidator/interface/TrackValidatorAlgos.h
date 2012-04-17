#ifndef TrackValidatorAlgos_h
#define TrackValidatorAlgos_h

// -*- C++ -*-
//
// Package:    TrackValidator
// Class:      TrackValidator
// 
/**\class TrackValidator TrackValidator.cc MGeisler/TrackValidator/src/TrackValidator.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Fri Feb  3 13:57:40 CET 2012
// $Id: TrackValidatorAlgos.h,v 1.6 2012/03/22 15:10:16 mgeisler Exp $
//
//

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

#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"

// ROOT include files
#include <TH1F.h>

using namespace std;
using namespace edm;
using namespace reco;

typedef vector<pair<TrackingParticleRef, double> > TpDoubV;

//
// class declaration
//

class TrackValidatorAlgos{
 public:

    TrackValidatorAlgos(const edm::ParameterSet&);

    void CreateIntervalVectors();

    void GetInputCollections(const Event&);

    void initialize(){setUpVectors();};

    void initializePF(){setUpVectorsPF();};

    void BookHistos(TFileDirectory); 

    void BookHistosPF(TFileDirectory); 

    void fill_independent_histos(int,int,int);

    void fill_generic_simTrack_histos(int,TrackingParticle*);

    void fill_recoAssociated_simTrack_histos(int,TrackingParticle*,const Track*,int);

    void fill_simAssociated_recoTrack_histos(int,const Track&,bool,bool,int,TpDoubV);

    void fill_removedRecoTrack_histos(int,const Track&,bool,bool,int,TpDoubV);

    void fill_photon_related_histos(int,vector<PhotonMCTruth>,auto_ptr<PFCandidateCollection>,
                                    auto_ptr<PFCandidateCollection>,SimVertex,int);

    void fillFractionHistosFromVectors(int);
    void fillFractionHistosFromVectorsPF(int);

    void fillHistosFromVectors(int);
    void fillHistosFromVectorsPF(int);

    bool findRefTrack(const Track&,const Track&);

    double getTrackWeight(const TrackingParticle*, const GenJetCollection*);

    bool isGenPart(const TrackingParticle*, const GenParticle*);

 protected:
  //protected functions 

 private: 

  // private methods for internal usage
  void setUpVectors();
  void setUpVectorsPF();

  void fillPlotFromVector(TH1F*,vector<int>);

  void fillFractionHisto(TH1F*,vector<int>,vector<int>,string);

  bool photonMatching(float, float, float, PFCandidate);

  bool photonSelector(PhotonMCTruth,SimVertex);

  float etaTransformation(float, float);


  //private data members  

  bool useLogpt_;

  bool useJetWeighting_;
  string genJetCollLabel_;

  Handle<GenJetCollection> genJetCollH;

  TrackingParticleSelector* generalTpSignalSelector;
  TrackingParticleSelector* generalTpPUSelector;

  //parameters for the histograms

  double minEta, maxEta;  int nintEta;
  double minpt, maxpt;  int nintpt;
  double minVertcount, maxVertcount;  int nintVertcount;
  double minTrackcount, maxTrackcount;  int nintTrackcount;

  //parameters for the photon selector

  float etaPhiDistance, BARL, END_LO, END_HI;
  float MCphotonPtMin_, MCphotonLip_, MCphotonTip_;


  // ###########
  // histograms 
  // ###########

  // track collection

  vector<TH1F*> PU_effic_eta; vector<TH1F*> PU_effic_pt; vector<TH1F*> PU_effic_npu;

  vector<TH1F*> PU_fakerate_1_eta; vector<TH1F*> PU_fakerate_1_pt; 
  vector<TH1F*> PU_fakerate_1_npu;
  vector<TH1F*> PU_fakerate_2_eta; vector<TH1F*> PU_fakerate_2_pt; 
  vector<TH1F*> PU_fakerate_2_npu;

  vector<TH1F*> effic_eta; vector<TH1F*> effic_pt; 
  vector<TH1F*> effic_npu;
  vector<TH1F*> fakerate_eta; vector<TH1F*> fakerate_pt; 
  vector<TH1F*> fakerate_npu;

  vector<TH1F*> num_simul_tracks; vector<TH1F*> num_track_simul_eta;
  vector<TH1F*> num_track_simul_pt; vector<TH1F*> num_track_simul_npu; 
  vector<TH1F*> num_simul_vertex; 
  vector<TH1F*> num_reco_tracks; vector<TH1F*> num_track_reco_eta;
  vector<TH1F*> num_track_reco_pt; vector<TH1F*> num_track_reco_npu;

  vector<TH1F*> num_removed_reco_signal_eta; vector<TH1F*> num_removed_reco_eta;
  vector<TH1F*> num_removed_reco_PU_eta; vector<TH1F*> num_reco_PU_eta;

  vector<TH1F*> num_removed_reco_signal_pt; vector<TH1F*> num_removed_reco_pt;
  vector<TH1F*> num_removed_reco_PU_pt;

  vector<TH1F*> num_removed_reco_signal_npu; vector<TH1F*> num_removed_reco_npu;
  vector<TH1F*> num_removed_reco_PU_npu;

  vector<TH1F*> num_assoc_eta; vector<TH1F*> num_assoc2_eta;
  vector<TH1F*> num_assoc_pt; vector<TH1F*> num_assoc2_pt;
  vector<TH1F*> num_assoc_npu; vector<TH1F*> num_assoc2_npu;

  vector<TH1F*> num_track_reco_PU_eta; vector<TH1F*> num_track_reco_signal_eta;
  vector<TH1F*> num_track_reco_PU_pt; vector<TH1F*> num_track_reco_signal_pt;
  vector<TH1F*> num_track_reco_PU_npu; vector<TH1F*> num_track_reco_signal_npu;


  // particle flow collection

  vector<TH1F*> num_photon_simul; vector<TH1F*> num_photon_simul_eta;
  vector<TH1F*> num_photon_simul_pt; vector<TH1F*> num_photon_simul_npu;

  vector<TH1F*> num_photon_reco; vector<TH1F*> num_photon_reco_eta;
  vector<TH1F*> num_photon_reco_pt; vector<TH1F*> num_photon_reco_npu;

  vector<TH1F*> num_assoc_photon_eta; vector<TH1F*> num_assoc2_photon_eta;
  vector<TH1F*> num_assoc_photon_pt; vector<TH1F*> num_assoc2_photon_pt;
  vector<TH1F*> num_assoc_photon_npu; vector<TH1F*> num_assoc2_photon_npu;

  vector<TH1F*> photon_effic_eta; vector<TH1F*> photon_effic_pt; 
  vector<TH1F*> photon_effic_npu;
  vector<TH1F*> photon_fakerate_eta; vector<TH1F*> photon_fakerate_pt; 
  vector<TH1F*> photon_fakerate_npu;

  vector<TH1F*> num_removedPURecoPhoton_eta, num_removedPURecoPhoton_npu;
  vector<TH1F*> num_removedPURecoPhoton_pt;
  vector<TH1F*> num_allPURecoPhoton_eta, num_allPURecoPhoton_pt;
  vector<TH1F*> num_allPURecoPhoton_npu;

  vector<TH1F*> photon_PU_effic_eta; vector<TH1F*> photon_PU_effic_pt; 
  vector<TH1F*> photon_PU_effic_npu;

  vector<TH1F*> num_removedSignalRecoPhoton_eta, num_removedSignalRecoPhoton_npu;
  vector<TH1F*> num_removedSignalRecoPhoton_pt;

  vector<TH1F*> num_allRemovedRecoPhoton_eta, num_allRemovedRecoPhoton_npu;
  vector<TH1F*> num_allRemovedRecoPhoton_pt;

  vector<TH1F*> num_allSignalRecoPhoton_eta, num_allSignalRecoPhoton_npu;
  vector<TH1F*> num_allSignalRecoPhoton_pt;

  vector<TH1F*> photon_PU_fakerate_1_eta; vector<TH1F*> photon_PU_fakerate_1_pt; 
  vector<TH1F*> photon_PU_fakerate_1_npu;
  vector<TH1F*> photon_PU_fakerate_2_eta; vector<TH1F*> photon_PU_fakerate_2_pt; 
  vector<TH1F*> photon_PU_fakerate_2_npu;


  // ###########
  // vectors 
  // ###########

  vector<double> etaintervals;
  vector<double> ptintervals;
  vector<double> vertcountintervals;

  // track collection

  vector< vector<int> > allSignalTP_eta, allRT_eta;
  vector< vector<int> > allSignalTP_npu, allRT_npu;
  vector< vector<int> > allSignalTP_pt,  allRT_pt;

  vector< vector<int> > assSignalTP_eta, assSignalTP_npu, assSignalTP_pt;
  vector< vector<int> > assSignalRT_eta, assSignalRT_npu, assSignalRT_pt;

  vector< vector<int> > allSigRT_eta, allPURT_eta, allAssPURT_eta;
  vector< vector<int> > allSigRT_pt,  allPURT_pt;
  vector< vector<int> > allSigRT_npu, allPURT_npu;

  vector< vector<int> > allRemovedRT_eta, allRemovedRT_npu, allRemovedRT_pt;
  vector< vector<int> > removedSigRT_eta, removedSigRT_npu, removedSigRT_pt;
  vector< vector<int> > removedPURT_eta,  removedPURT_npu,  removedPURT_pt;

  vector<int> sim_tracks;


  // particle flow collection

  vector< vector<int> > allSignalPhoton_eta, allRecoPhoton_eta;
  vector< vector<int> > allSignalPhoton_npu, allRecoPhoton_npu;
  vector< vector<int> > allSignalPhoton_pt, allRecoPhoton_pt;

  vector< vector<int> > assSignalPhoton_eta, signalRecoPhoton_eta;
  vector< vector<int> > assSignalPhoton_npu, signalRecoPhoton_npu;
  vector< vector<int> > assSignalPhoton_pt, signalRecoPhoton_pt;

  vector< vector<int> > removedPURecoPhoton_eta, removedPURecoPhoton_npu, removedPURecoPhoton_pt;
  vector< vector<int> > allPURecoPhoton_eta, allPURecoPhoton_pt, allPURecoPhoton_npu;

  vector< vector<int> > removedSignalRecoPhoton_eta, removedSignalRecoPhoton_npu;
  vector< vector<int> > removedSignalRecoPhoton_pt;
  vector< vector<int> > allRemovedRecoPhoton_eta, allRemovedRecoPhoton_npu, allRemovedRecoPhoton_pt;
  vector< vector<int> > allSignalRecoPhoton_eta, allSignalRecoPhoton_npu, allSignalRecoPhoton_pt;

};

#endif
