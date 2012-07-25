// -*- C++ -*-
//
// Package:    JetAnlzr
// Class:      JetAnlzr
// 
/**\class JetAnlzr JetAnlzr.cc MGeisler/JetAnlzr/src/JetAnlzr.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Mon Jan 23 11:53:44 CET 2012
// $Id: JetAnlzr.cc,v 1.3 2012/03/26 12:32:56 mgeisler Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unistd.h>
#include <stdio.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/src/Utilities.cc"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//ROOT include files
#include "TProfile3D.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TFile.h"

using namespace edm;
using namespace std;
using namespace reco;

typedef math::PtEtaPhiMLorentzVectorD TLV; 

//
// class declaration
//

class JetAnlzr : public edm::EDAnalyzer {
   public:
      explicit JetAnlzr(const edm::ParameterSet&);
      ~JetAnlzr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual bool Match(TLV, TLV);

      virtual void endRun(Run const&, EventSetup const&);

      // ----------member data ---------------------------

      string input_GenJets_;
      std::vector<string> input_RecoJets_;
      string label_pileupinfo_;

      string outname;
//       vector<FILE*> datei;

      char dirName_[256];

      vector<TH3F*> ResponseVsPtVsEta_NpuCr, ResponseVsPtVsEta_Npu15;
      vector<TH3F*> ResponseVsPtVsEta_Npu30, ResponseVsPtVsEta_Npu50;

      vector<TH3F*> ResolutionVsPtVsEta_NpuCr, ResolutionVsPtVsEta_Npu15;
      vector<TH3F*> ResolutionVsPtVsEta_Npu30, ResolutionVsPtVsEta_Npu50;

      vector<TH3F*> SimJetsVsPtVsEtaVsNpu, RecoJetsVsPtVsEtaVsNpu;
      vector<TH3F*> AssocJetsVsPtVsEtaVsNpu, Assoc2JetsVsPtVsEtaVsNpu;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetAnlzr::JetAnlzr(const edm::ParameterSet& iConfig)
{
   	//now do what ever initialization is needed

  	outname = iConfig.getParameter<string>("Outname");

       	Service<TFileService> tfserv;
  	vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  	input_GenJets_ = iConfig.getParameter<string>("genJets");
  	input_RecoJets_ = iConfig.getParameter<vector<string> >("recoJets");
  	label_pileupinfo_ = iConfig.getParameter<string>("PileUpInfo");

	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

    	  string dirName = input_RecoJets_[rjc_ite];
	  dirName.erase(dirName.length()-1,1);

          subDir->push_back(tfserv->mkdir(dirName));

	  ResponseVsPtVsEta_NpuCr.push_back(subDir->at(rjc_ite).make<TH3F>("ResponseVsPtVsEta_NpuCr","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 100, 0., 2.0));
	  ResponseVsPtVsEta_Npu15.push_back(subDir->at(rjc_ite).make<TH3F>("ResponseVsPtVsEta_Npu15","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 100, 0., 2.0));

	  ResponseVsPtVsEta_Npu30.push_back(subDir->at(rjc_ite).make<TH3F>("ResponseVsPtVsEta_Npu30","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 100, 0., 2.0));
	  ResponseVsPtVsEta_Npu50.push_back(subDir->at(rjc_ite).make<TH3F>("ResponseVsPtVsEta_Npu50","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 100, 0., 2.0));


	  ResolutionVsPtVsEta_NpuCr.push_back(subDir->at(rjc_ite).make<TH3F>("ResolutionVsPtVsEta_NpuCr","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 50, -1.0, 1.0));
	  ResolutionVsPtVsEta_Npu15.push_back(subDir->at(rjc_ite).make<TH3F>("ResolutionVsPtVsEta_Npu15","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 50, -1.0, 1.0));

	  ResolutionVsPtVsEta_Npu30.push_back(subDir->at(rjc_ite).make<TH3F>("ResolutionVsPtVsEta_Npu30","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 50, -1.0, 1.0));
	  ResolutionVsPtVsEta_Npu50.push_back(subDir->at(rjc_ite).make<TH3F>("ResolutionVsPtVsEta_Npu50","pt response vs pt vs eta", 50, 0., 1000., 50, -5., 5., 50, -1.0, 1.0));

	  SimJetsVsPtVsEtaVsNpu.push_back(subDir->at(rjc_ite).make<TH3F>("SimJetsVsPtVsEtaVsNpu","Simulated jets vs pt vs eta vs npu", 50, 0., 1000., 50, -5., 5., 50, -0.5, 49.5));
	  RecoJetsVsPtVsEtaVsNpu.push_back(subDir->at(rjc_ite).make<TH3F>("RecoJetsVsPtVsEtaVsNpu","Reconstructed jets vs pt vs eta vs npu", 50, 0., 1000., 50, -5., 5., 50, -0.5, 49.5));

	  AssocJetsVsPtVsEtaVsNpu.push_back(subDir->at(rjc_ite).make<TH3F>("AssocJetsVsPtVsEtaVsNpu","Sim To Reco associated jets vs pt vs eta vs npu", 50, 0., 1000., 50, -5., 5., 50, -0.5, 49.5));
	  Assoc2JetsVsPtVsEtaVsNpu.push_back(subDir->at(rjc_ite).make<TH3F>("Assoc2JetsVsPtVsEtaVsNpu","Reco To Sim associated jets vs pt vs eta vs npu", 50, 0., 1000., 50, -5., 5., 50, -0.5, 49.5));


// 	  string help = outname + "_" + dirName + ".txt";

// 	  const char * c = help.c_str();
//    	  datei.push_back(fopen(c , "w"));

	}

}


JetAnlzr::~JetAnlzr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
bool 
JetAnlzr::Match(TLV jet1, TLV jet2)
{

	double minDR = 0.3;

 	double dR = deltaR(jet1,jet2);
	 
	return dR<minDR;

}

// ------------ method called for each event  ------------
void
JetAnlzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    	Service<TFileService> tfserv;

	// get the number of pileup interactions
  
  	Handle< vector<PileupSummaryInfo> > puinfoH;
  	iEvent.getByLabel(label_pileupinfo_,puinfoH);
  	PileupSummaryInfo puinfo;      
  
  	for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    	    if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      	    puinfo=(*puinfoH)[puinfo_ite];
      	    break;
    	  }
  	}

	unsigned npu = puinfo.getPU_NumInteractions();

	// get the information for the gen jets

	Handle<vector<float> > genJetPtH;
    	iEvent.getByLabel( input_GenJets_,"pt",genJetPtH);

	Handle<vector<float> > genJetEtaH;
    	iEvent.getByLabel( input_GenJets_,"eta",genJetEtaH);

	Handle<vector<float> > genJetPhiH;
	iEvent.getByLabel( input_GenJets_,"phi",genJetPhiH);

	string genJetsLabel = input_GenJets_;
        genJetsLabel.erase(genJetsLabel.size()-1,1); 

	Handle<GenJetCollection > gJH;
	iEvent.getByLabel( genJetsLabel,gJH);

    	vector<TLV> genJets;
	vector<GenJet> gJ;

 	for(unsigned gen_ite=0; gen_ite<genJetPtH->size(); gen_ite++){

	  TLV genJet( genJetPtH->at(gen_ite),
		      genJetEtaH->at(gen_ite),
		      genJetPhiH->at(gen_ite),
		      0.0 );
	
	  genJets.push_back(genJet);
	  gJ.push_back(gJH->at(gen_ite));

	}

	// loop over all reco jet collections
	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

	  // get the information for the reco jets

	  vector<double> ptReco1D;

	  Handle<vector<float> > jetPtH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"pt",jetPtH);

	  Handle<vector<float> > jetEtaH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"eta",jetEtaH);

	  Handle<vector<float> > jetPhiH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"phi",jetPhiH);

	  string recoJetsLabel = input_RecoJets_[rjc_ite];
          recoJetsLabel.erase(recoJetsLabel.size()-1,1); 

	  Handle<PFJetCollection> rJH;
	  iEvent.getByLabel( recoJetsLabel,rJH);
 
	  // loop over sim jet collection
	  //  - efficiency calculation
	  //  - response && resolution calculation
    	  for(unsigned gen_ite=0; gen_ite<genJets.size();gen_ite++){

	    TLV genJet = genJets.at(gen_ite);
	    SimJetsVsPtVsEtaVsNpu[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), npu);

	    for(unsigned reco_ite=0; reco_ite<jetPtH->size(); reco_ite++){
	
	      TLV recoJet( jetPtH->at(reco_ite),
		           jetEtaH->at(reco_ite),
		           jetPhiH->at(reco_ite),
		           0.0 );	

	      if(Match(genJet,recoJet)){

	        AssocJetsVsPtVsEtaVsNpu[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), npu);

	        //response && resolution
	        if((genJet.Pt()>30.0) && (gen_ite<2)){

                  ResponseVsPtVsEta_NpuCr[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), recoJet.Pt() *1./genJet.Pt());
                  ResolutionVsPtVsEta_NpuCr[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());

                  if(npu<16){
		    ResponseVsPtVsEta_Npu15[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), recoJet.Pt() *1./genJet.Pt());
		    ResolutionVsPtVsEta_Npu15[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());
	          }

                  if((npu>15)&&(npu<31)){
                    ResponseVsPtVsEta_Npu30[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), recoJet.Pt() *1./genJet.Pt());
                    ResolutionVsPtVsEta_Npu30[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());
	          }

                  if(npu>30){ 
                    ResponseVsPtVsEta_Npu50[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), recoJet.Pt() *1./genJet.Pt());
                    ResolutionVsPtVsEta_Npu50[rjc_ite]->Fill(genJet.Pt(), genJet.Eta(), (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());
		  }

// 		  fprintf(datei[rjc_ite],"%f\t%f\t%f\t%i\n",genJet.Pt(),recoJet.Pt(),genJet.Eta(),npu);

	        }

	        break;

	      }

	    }

	  }

	  // loop over reco jet collection
	  //  - fakerate calculation
 	  for(unsigned reco_ite=0; reco_ite<jetPtH->size(); reco_ite++){

	    TLV recoJet( jetPtH->at(reco_ite),
		         jetEtaH->at(reco_ite),
			 jetPhiH->at(reco_ite),
			 0.0 );	 

	    RecoJetsVsPtVsEtaVsNpu[rjc_ite]->Fill(recoJet.Pt(), recoJet.Eta(), npu);

    	    for(unsigned gen_ite=0; gen_ite<genJets.size();gen_ite++){

	      TLV genJet = genJets.at(gen_ite);

              if(Match(genJet,recoJet)){

	        Assoc2JetsVsPtVsEtaVsNpu[rjc_ite]->Fill(recoJet.Pt(), recoJet.Eta(), npu);

	        break;

	      }

	    }

	  }

	}

}

// ------------ method called when ending the processing of a run  ------------
void 
JetAnlzr::endRun(edm::Run const&, edm::EventSetup const&)
{
// 	for(unsigned i=0; i<datei.size(); i++){
// 	  fclose(datei[i]);
//         }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnlzr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnlzr);
