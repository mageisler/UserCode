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
// $Id: JetAnlzr.cc,v 1.2 2012/01/25 17:00:32 mgeisler Exp $
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

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/src/Utilities.cc"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//ROOT include files
#include "TH3F.h"
#include "TH2F.h"
#include "TMath.h"

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
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual bool findGenJet(TLV, vector<TLV>, TLV*);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      // ----------member data ---------------------------

      string input_GenJets_;
      std::vector<string> input_RecoJets_;
      string label_pileupinfo_;
      string label_JetCorrector_;

      char dirName_[256];

      vector<TH2F*> WithVsPt, WithVsEta, WithVsNpu;
      vector<TH2F*> WithOutVsPt, WithOutVsEta, WithOutVsNpu;

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


       	Service<TFileService> tfserv;
  	vector<TFileDirectory>* subDir(new vector<TFileDirectory>());

  	input_GenJets_ = iConfig.getParameter<string>("genJets");
  	input_RecoJets_ = iConfig.getParameter<vector<string> >("recoJets");
  	label_pileupinfo_ = iConfig.getParameter<string>("PileUpInfo");
  	label_JetCorrector_ = iConfig.getParameter<string>("JetCorrector");

	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

    	  string dirName = input_RecoJets_[rjc_ite];
	  dirName.erase(dirName.length()-1,1);

          subDir->push_back(tfserv->mkdir(dirName));

	  WithVsPt.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithFactorVsPt","pt ratio with factor vs pt", 50, 0., 1000., 200, 0., 2.0));
	  WithVsEta.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithFactorVsEta","pt ratio with factor vs eta", 50, -5., 5., 200, 0., 2.0));
	  WithVsNpu.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithFactorVsNpu","pt ratio with factor vs npu", 50, -0.5, 49.5, 200, 0., 2.0));

	  WithOutVsPt.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithOutFactorVsPt","pt ratio without factor vs pt", 50, 0., 1000., 200, 0., 2.0));
	  WithOutVsEta.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithOutFactorVsEta","pt ratio without factor vs eta", 50, -5., 5., 200, 0., 2.0));
	  WithOutVsNpu.push_back(subDir->at(rjc_ite).make<TH2F>("ptRatioWithOutFactorVsNpu","pt ratio without factor vs npu", 50, -0.5, 49.5, 200, 0., 2.0));

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
JetAnlzr::findGenJet(TLV recoJet, vector<TLV> genJets, TLV* genJet)
{

	double minDR = 0.5;

    	for(unsigned gen_ite=0; gen_ite<genJets.size();gen_ite++){

	  double dR = deltaR(recoJet,genJets.at(gen_ite));
	 
	  if(dR<minDR){	    
	    *genJet = genJets.at(gen_ite);
	    return true;
	  }

	}    

	return false;

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


    	vector<TLV> genJets;

 	for(unsigned gen_ite=0; gen_ite<genJetPtH->size(); gen_ite++){
 	  
	  if(genJetPtH->at(gen_ite)<30.0) continue;

	  TLV genJet( genJetPtH->at(gen_ite),
		      genJetEtaH->at(gen_ite),
		      genJetPhiH->at(gen_ite),
		      0.0 );
	
	  genJets.push_back(genJet);

	}

	// loop over all reco jet collections
	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

	  string jetLabel = input_RecoJets_[rjc_ite];
	  jetLabel.erase(jetLabel.length()-1,1);

	  // get the information for the reco jets
        
	  Handle<PFJetCollection> jets;        			//define input jet collection
	  iEvent.getByLabel( jetLabel, jets);    		//get input jet collection

	  const JetCorrector* corrector = JetCorrector::getJetCorrector(label_JetCorrector_,iSetup);  

	  Handle<vector<float> > jetPtH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"pt",jetPtH);

	  Handle<vector<float> > jetEtaH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"eta",jetEtaH);

	  Handle<vector<float> > jetPhiH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"phi",jetPhiH);

	  unsigned numJets = 2;
	  if(jets->size()<2) numJets = jets->size();

 	  for(unsigned reco_ite=0; reco_ite<numJets; reco_ite++){
 	  
	    if(jetPtH->at(reco_ite)<30.0) continue;  

	    TLV recoJet( jetPtH->at(reco_ite),
		         jetEtaH->at(reco_ite),
			 jetPhiH->at(reco_ite),
			 0.0 );	 

            double factor = corrector->correction(jets->at(reco_ite),iEvent,iSetup);

	    TLV genJet;
            if (findGenJet(recoJet,genJets,&genJet)){

              WithVsPt[rjc_ite]->Fill(genJet.Pt(), factor * recoJet.Pt() /genJet.Pt());
              WithVsEta[rjc_ite]->Fill(genJet.Eta(), factor * recoJet.Pt() /genJet.Pt());
              WithVsNpu[rjc_ite]->Fill(npu, factor * recoJet.Pt() /genJet.Pt());

              WithOutVsPt[rjc_ite]->Fill(genJet.Pt(), recoJet.Pt() /genJet.Pt());
              WithOutVsEta[rjc_ite]->Fill(genJet.Eta(), recoJet.Pt() /genJet.Pt());
              WithOutVsNpu[rjc_ite]->Fill(npu, recoJet.Pt() /genJet.Pt());

	    }

	  }

	}

}


// ------------ method called once each job just before starting event loop  ------------
void 
JetAnlzr::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetAnlzr::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
JetAnlzr::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}

// ------------ method called when ending the processing of a run  ------------
void 
JetAnlzr::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
JetAnlzr::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
JetAnlzr::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
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
