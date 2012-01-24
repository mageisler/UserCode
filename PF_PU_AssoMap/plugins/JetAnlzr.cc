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
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
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

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//ROOT include files
#include "TH3F.h"
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

      char dirName_[256];
      FactorizedJetCorrector* jec;

//       vector<TH3F*> ptRatioHistos;
//       vector<TH3F*> etaRatioHistos;
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

  	input_GenJets_ = iConfig.getParameter<string>("genJets");
  	input_RecoJets_ = iConfig.getParameter<vector<string> >("recoJets");
  	label_pileupinfo_ = iConfig.getParameter<string>("PileUpInfo");

    	Service<TFileService> tfserv;


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

	  char etaName[256];
	  sprintf(etaName,"etaRecoRatio_%s",input_RecoJets_[rjc_ite].c_str());
  	  TH3F* etaRatioHisto = tfserv->make<TH3F>(etaName,etaName, 50, -5.0, 5.0, 200, 0., 2.0, 5, 0, 25);

	  char ptName[256];
	  sprintf(ptName,"ptRecoRatio_%s",input_RecoJets_[rjc_ite].c_str());
  	  TH3F* ptRatioHisto = tfserv->make<TH3F>(ptName,ptName, 50, 0., 500., 200, 0., 2.0, 5, 0, 25);

	  // get the information for the reco jets

	  Handle<vector<float> > jetPtH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"pt",jetPtH);

	  Handle<vector<float> > jetEtaH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"eta",jetEtaH);

	  Handle<vector<float> > jetPhiH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"phi",jetPhiH);

	  Handle<vector<float> > jetAreaH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"jetArea",jetAreaH);

	  Handle<vector<float> > jetRhoH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"rho",jetRhoH);
     	  float rho = jetRhoH->at(0);

 	  for(unsigned reco_ite=0; reco_ite<2; reco_ite++){
 	  
	    if(jetPtH->at(reco_ite)<30.0) continue;

	    TLV recoJet( jetPtH->at(reco_ite),
		         jetEtaH->at(reco_ite),
			 jetPhiH->at(reco_ite),
			 0.0 );	    

            jec->setJetEta(jetEtaH->at(reco_ite));
            jec->setJetPt(jetPtH->at(reco_ite));
            jec->setJetA(jetAreaH->at(reco_ite));
            jec->setRho(rho);
            double factor = jec->getCorrection();

	    TLV genJet;
            if (findGenJet(recoJet,genJets,&genJet)){
              ptRatioHisto->Fill(genJet.Pt(), factor * recoJet.Pt() /genJet.Pt(),  npu);
              etaRatioHisto->Fill(genJet.Eta(), factor * recoJet.Pt() /genJet.Pt(),  npu);
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

	ifstream jecStr("Jec11_V3_All.txt");

	char line[256];

	vector<JetCorrectorParameters> jecPars;

 	if (jecStr.is_open()){
    	  while (jecStr.good()){
            if (jecStr.eof()) break;
      	    jecStr.getline(line,256);
    	    JetCorrectorParameters ijec(line);
    	    jecPars.push_back(ijec);
    	  }
	}
    
	jec = new FactorizedJetCorrector(jecPars); 
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
