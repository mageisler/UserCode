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
// $Id$
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

// ROOT include files
#include <TH1F.h>

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class TrackValidatorAlgos{
 public:

    TrackValidatorAlgos(const edm::ParameterSet&);

    void initialize(){setUpVectors();};

    void BookHistos(TFileDirectory); 

    void fill_independent_histos(int,int,int);

    void fill_generic_simTrack_histos(int,TrackingParticle*);

    void fill_recoAssociated_simTrack_histos(int,TrackingParticle*,const Track*,int);

    void fill_simAssociated_recoTrack_histos(int,const Track&,vector<pair<TrackingParticleRef, double> >,int);

    void fill_removedRecoTrack_histos(int,const Track&,vector<pair<TrackingParticleRef, double> >,bool,int);

    void fillFractionHistosFromVectors(int);

    void fillHistosFromVectors(int);

    void fillPlotFromVector(TH1F*,vector<int>);

    void fillFractionHisto(TH1F*,vector<int>,vector<int>,string);

    bool findRefTrack(const Track&,const Track&);

 protected:
  //protected functions 

 private: 

  // private methods for internal usage
  void setUpVectors();


  //private data members  

  bool useLogpt_;

  TrackingParticleSelector* generalTpSignalSelector;
  TrackingParticleSelector* generalTpPUSelector;

  //parameters for the histograms

  double minEta, maxEta;  int nintEta;
  double minpt, maxpt;  int nintpt;
  double minVertcount, maxVertcount;  int nintVertcount;
  double minTrackcount, maxTrackcount;  int nintTrackcount;


  //histograms

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
  vector<TH1F*> num_simul_vertex; 
  vector<TH1F*> num_reco_tracks; vector<TH1F*> num_track_reco_eta;

  vector<TH1F*> num_removed_reco_signal_eta; vector<TH1F*> num_removed_reco_eta;
  vector<TH1F*> num_removed_reco_PU_eta; vector<TH1F*> num_reco_PU_eta;

  vector<TH1F*> num_removed_reco_signal_pt; vector<TH1F*> num_removed_reco_pt;
  vector<TH1F*> num_removed_reco_PU_pt;

  vector<TH1F*> num_removed_reco_signal_npu; vector<TH1F*> num_removed_reco_npu;
  vector<TH1F*> num_removed_reco_PU_npu;

  vector<TH1F*> num_assoc_eta; vector<TH1F*> num_assoc2_eta;


  //vectors

  vector< vector<double> > etaintervals;
  vector< vector<double> > ptintervals;
  vector< vector<double> > vertcountintervals;

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

};

#endif
