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

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

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
  maxpt  = 1000.;
  nintpt = 40;
  
  //parameters for Pileup plots
  minVertcount  = 0.;
  maxVertcount  = 50.;
  nintVertcount = 50;
  
  //parameters for track number plots
  minTrackcount  = -0.5;
  maxTrackcount  = 499.5;
  nintTrackcount = 500;

  //configure TP selectors

  using namespace reco::modules;

  ParameterSet generalTpSignalSelectorPSet = iConfig.getParameter<ParameterSet>("generalTpSelector");

  generalTpSignalSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpSignalSelectorPSet));

  ParameterSet generalTpPUSelectorPSet = generalTpSignalSelectorPSet;
  Entry name("signalOnly",false,true);
  generalTpPUSelectorPSet.insert(true,"signalOnly",name);

  generalTpPUSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(generalTpPUSelectorPSet));

  // fix for the LogScale by Ryan
  useLogpt_ = iConfig.getParameter<bool>("UseLogPt");

  if(useLogpt_){
    maxpt=log10(maxpt);
    if(minpt > 0){
      minpt=log10(minpt);
    }else{
      minpt=log10(0.1);
    }
  }

}


void 
TrackValidatorAlgos::setUpVectors()
{

  // eta vectors
  vector<double> etaintervalsv;
  vector<int> etaintervalsh;

  double eta_step=(maxEta-minEta)/nintEta;
  etaintervalsv.push_back(minEta);

  for (int k=1;k<nintEta+1;k++) {

    double d=minEta+k*eta_step;
    etaintervalsv.push_back(d);
    etaintervalsh.push_back(0);

  }

  etaintervals.push_back(etaintervalsv);
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
  vector<double> ptintervalsv;
  vector<int> ptintervalsh;

  double pt_step=(maxpt-minpt)/nintpt;
  ptintervalsv.push_back(minpt);

  for (int k=1;k<nintpt+1;k++) {

    double d;
    if(useLogpt_){
      d=pow(10,minpt+k*pt_step);
    }else{
      d=minpt+k*pt_step;
    }
    ptintervalsv.push_back(d);
    ptintervalsh.push_back(0);

  }

  ptintervals.push_back(ptintervalsv);
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
  vector<double> vertcountintervalsv;
  vector<int> vertcountintervalsh;

  double stepVertcount=(maxVertcount-minVertcount)/nintVertcount;
  vertcountintervalsv.push_back(minVertcount);

  for (int k=1;k<nintVertcount+1;k++) {

    double d=minVertcount+k*stepVertcount;
    vertcountintervalsv.push_back(d);
    vertcountintervalsh.push_back(0);

  }   

  vertcountintervals.push_back(vertcountintervalsv);
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
TrackValidatorAlgos::BookHistos(TFileDirectory subDir) 
{

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
  num_simul_vertex.push_back(subDir.make<TH1F>("num_simul_vertex", "Number of simulated vertices", nintVertcount, minVertcount, maxVertcount));

  //Book reconstruction related histograms

  num_reco_tracks.push_back(subDir.make<TH1F>("num_reco_tracks", "Number of reconstructed tracks", nintTrackcount, minTrackcount, maxTrackcount));
  num_track_reco_eta.push_back(subDir.make<TH1F>("num_track_reco_eta", "Number of reconstructed tracks vs eta", nintEta, minEta, maxEta));

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

  //Book association related histograms

  num_assoc_eta.push_back(subDir.make<TH1F>("num_assoc(simToReco)_eta", "Number of associated simulated tracks vs eta", nintEta, minEta, maxEta));
  num_assoc2_eta.push_back(subDir.make<TH1F>("num_assoc(recoToSim)_eta", "Number of associated reconstructed tracks vs eta", nintEta, minEta, maxEta));

}

void 
TrackValidatorAlgos::fill_independent_histos(int counter, int npu, int rt)
{

  num_simul_vertex.at(counter)->Fill(npu);
  num_simul_tracks.at(counter)->Fill(sim_tracks[counter]);
  num_reco_tracks.at(counter)->Fill(rt);

}

void 
TrackValidatorAlgos::fill_generic_simTrack_histos(int counter, TrackingParticle* tp)
{

  if((*generalTpSignalSelector)(*tp)){
    num_track_simul_eta.at(counter)->Fill(tp->momentum().eta());
  }

}

void
TrackValidatorAlgos::fill_recoAssociated_simTrack_histos(int counter, TrackingParticle* tp, const Track* track, int npu)
{

  bool isMatched = track;
  double tp_eta = tp->momentum().eta();

  if((*generalTpSignalSelector)(*tp)){

    sim_tracks[counter]++;

    //effic vs eta
    for (unsigned int f=0; f<etaintervals[counter].size()-1; f++){
      if (tp_eta>etaintervals[counter][f]&&
	  tp_eta<etaintervals[counter][f+1]) {
	allSignalTP_eta[counter][f]++;
	if (isMatched) {
	  assSignalTP_eta[counter][f]++;
	}
        break;
      }
    } // END for (unsigned int f=0; f<etaintervals[w].size()-1; f++){

    //effic vs pt
    for (unsigned int f=0; f<ptintervals[counter].size()-1; f++){
      if (sqrt(tp->momentum().perp2())>ptintervals[counter][f]&&
	  sqrt(tp->momentum().perp2())<ptintervals[counter][f+1]) {
        allSignalTP_pt[counter][f]++; 
        if (isMatched) {
	  assSignalTP_pt[counter][f]++;
        }	
        break;      
      }
    } // End for (unsigned int f=0; f<ptintervals[count].size()-1; f++){

    //effic vs num pileup vertices
    for (unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){
      if (npu == vertcountintervals[counter][f]) {
        allSignalTP_npu[counter][f]++;
        if (isMatched) {
          assSignalTP_npu[counter][f]++;
        }
        break;
      }    
    }// End for (unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){

  }


}

void 
TrackValidatorAlgos::fill_simAssociated_recoTrack_histos(int counter,const Track& track,vector<pair<TrackingParticleRef, double> > tp,int npu)
{

  bool isMatched = false;
  bool isSigMatched = false;

  if (tp.size()!=0) {
    
    isMatched = true;
    for (unsigned int tp_ite=0;tp_ite<tp.size();++tp_ite){ 

      TrackingParticle trackpart = *(tp[tp_ite].first);
      if ((*generalTpSignalSelector)(trackpart)){
        isSigMatched = true;
        break;
      }

    }
  }

  //fake rate vs eta
  for (unsigned int f=0; f<etaintervals[counter].size()-1; f++){
    if (track.eta()>etaintervals[counter][f]&&
        track.eta()<etaintervals[counter][f+1]) {
      allRT_eta[counter][f]++;
      if (isSigMatched) {
	assSignalRT_eta[counter][f]++;
      }
      break;
    }
  } // END for (unsigned int f=0; f<etaintervals[w].size()-1; f++){

  //fake rate vs pt
  for (unsigned int f=0; f<ptintervals[counter].size()-1; f++){
    if (sqrt(track.momentum().perp2())>ptintervals[counter][f]&&
        sqrt(track.momentum().perp2())<ptintervals[counter][f+1]) {
      allRT_pt[counter][f]++; 
      if (isSigMatched) {
	assSignalRT_pt[counter][f]++;
      }	  
      break;    
    }
  } // End for (unsigned int f=0; f<ptintervals[count].size()-1; f++){

  //fake rate vs num pileup vertices
  for (unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){
    if (npu == vertcountintervals[counter][f]) {
      allRT_npu[counter][f]++;
      if (isSigMatched) {
        assSignalRT_npu[counter][f]++;
      }
      break;
    }    
  }// End for (unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){

  //pileup related

  if(isMatched){
 
    TrackingParticle trackpart = *(tp[0].first);

    if(((*generalTpPUSelector)(trackpart)) && (!(*generalTpSignalSelector)(trackpart))){

      // vs eta
      for (unsigned int f=0; f<etaintervals[counter].size()-1; f++){
        if (track.eta()>etaintervals[counter][f]&&
            track.eta()<etaintervals[counter][f+1]) {
          allAssPURT_eta[counter][f]++;
        }
      } // END for (unsigned int f=0; f<etaintervals[w].size()-1; f++){

    }// END for (((*generalTpPUSelector)(trackpart)) && (!(*generalTpSignalSelector)(trackpart)))

  } // END for (isMatched)


}

void 
TrackValidatorAlgos::fill_removedRecoTrack_histos(int counter,const Track& refTrack,vector<pair<TrackingParticleRef, double> > refTp,bool isRemoved,int npu)
{

  bool isMatched = false;
  bool isSigMatched = false;

  if (refTp.size()!=0) {
    for (unsigned int tp_ite=0;tp_ite<refTp.size();++tp_ite){ 

      TrackingParticle trackpart = *(refTp[tp_ite].first);
      if ((*generalTpPUSelector)(trackpart)){
    
        isMatched = true;
        if ((*generalTpSignalSelector)(trackpart))isSigMatched = true;
        break;
      }

    }
  }

  if(isMatched){ // means that a TP has been found, no matter of signal or not
 
    TrackingParticle trackpart = *(refTp[0].first);

    // vs eta
    for (unsigned int f=0; f<etaintervals[counter].size()-1; f++){
      if (refTrack.eta()>etaintervals[counter][f]&&
          refTrack.eta()<etaintervals[counter][f+1]) {
        if(isSigMatched){
          allSigRT_eta[counter][f]++; 
        }else{
  	  allPURT_eta[counter][f]++; 
        }   
        if(isRemoved){
          allRemovedRT_eta[counter][f]++;
          if(isSigMatched){
            removedSigRT_eta[counter][f]++; 
          }else{
  	    removedPURT_eta[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(refTrack.eta()>etaintervals[counter][f]&&refTrack.eta()<etaintervals[counter][f+1){
    } // END for(unsigned int f=0; f<etaintervals[w].size()-1; f++){

    // vs pt
    for (unsigned int f=0; f<ptintervals[counter].size()-1; f++){
      if (sqrt(refTrack.momentum().perp2())>ptintervals[counter][f]&&
          sqrt(refTrack.momentum().perp2())<ptintervals[counter][f+1]) {
        if(isSigMatched){
          allSigRT_pt[counter][f]++; 
        }else{
  	  allPURT_pt[counter][f]++; 
        }   
        if(isRemoved){
          allRemovedRT_pt[counter][f]++;
          if(isSigMatched){
            removedSigRT_pt[counter][f]++; 
          }else{
  	    removedPURT_pt[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(sqrt(refTrack.momentum().perp2())>ptintervals[counter][f]&&sqrt(refTrack.momentum().perp2())<ptintervals[counter][f+1]){
    } // END for(unsigned int f=0; f<ptintervals[counter].size()-1; f++){

    // vs npu
    for (unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){
      if (npu == vertcountintervals[counter][f]) {
        if(isSigMatched){
          allSigRT_npu[counter][f]++; 
        }else{
  	  allPURT_npu[counter][f]++; 
        }   
        if(isRemoved){
          allRemovedRT_npu[counter][f]++;
          if(isSigMatched){
            removedSigRT_npu[counter][f]++; 
          }else{
  	    removedPURT_npu[counter][f]++; 
          }
        } // END if(isRemoved){
      } // END if(npu == vertcountintervals[counter][f]){
    } // END for(unsigned int f=0; f<vertcountintervals[counter].size()-1; f++){

  }// END for (isMatched)
 
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

}

void
TrackValidatorAlgos::fillHistosFromVectors(int counter)
{

  fillPlotFromVector(num_track_simul_eta[counter],allSignalTP_eta[counter]);
  fillPlotFromVector(num_assoc_eta[counter],assSignalTP_eta[counter]);

  fillPlotFromVector(num_track_reco_eta[counter],allRT_eta[counter]);
  fillPlotFromVector(num_assoc2_eta[counter],assSignalRT_eta[counter]);

  fillPlotFromVector(num_removed_reco_signal_eta[counter],removedSigRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_eta[counter],allRemovedRT_eta[counter]);
  fillPlotFromVector(num_removed_reco_PU_eta[counter],removedPURT_eta[counter]);

  fillPlotFromVector(num_removed_reco_signal_pt[counter],removedSigRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_pt[counter],allRemovedRT_pt[counter]);
  fillPlotFromVector(num_removed_reco_PU_pt[counter],removedPURT_pt[counter]);

  fillPlotFromVector(num_removed_reco_signal_npu[counter],removedSigRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_npu[counter],allRemovedRT_npu[counter]);
  fillPlotFromVector(num_removed_reco_PU_npu[counter],removedPURT_npu[counter]);

  fillPlotFromVector(num_reco_PU_eta[counter],allAssPURT_eta[counter]);

}

void
TrackValidatorAlgos::fillPlotFromVector(TH1F* histo,vector<int> vInt)
{

  for(unsigned ite=0; ite<vInt.size();ite++){
    histo->SetBinContent(ite+1,vInt.at(ite));
  } 

}

void
TrackValidatorAlgos::fillFractionHisto(TH1F* histo,vector<int> num,vector<int> denum,string type)
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