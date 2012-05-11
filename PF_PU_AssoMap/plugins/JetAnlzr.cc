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
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual bool findGenJet(TLV, vector<TLV>, TLV*);
      virtual vector<int> GetSortedList(vector<float>);

      // ----------member data ---------------------------

      string input_GenJets_;
      std::vector<string> input_RecoJets_;
      string label_pileupinfo_;

      char dirName_[256];

      vector<TH3F*> ResponseVsPt, ResponseVsEta;
      vector<TH3F*> ResolutionVsPt, ResolutionVsEta;

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

	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

    	  string dirName = input_RecoJets_[rjc_ite];
	  dirName.erase(dirName.length()-1,1);

          subDir->push_back(tfserv->mkdir(dirName));

	  ResponseVsPt.push_back(subDir->at(rjc_ite).make<TH3F>("ptResponseVsPtVsNpu","pt response vs pt vs npu", 50, 0., 1000., 50, -0.5, 49.5, 100, 0., 2.0));
	  ResponseVsEta.push_back(subDir->at(rjc_ite).make<TH3F>("ptResponseVsEtaVsNpu","pt response vs eta vs npu", 50, -5., 5., 50, -0.5, 49.5, 100, 0., 2.0));

	  ResolutionVsPt.push_back(subDir->at(rjc_ite).make<TH3F>("ptResolutionVsPtVsNpu","pt resolution vs pt vs npu", 50, 0., 1000., 50, -0.5, 49.5, 50, -1.0, 1.0));
	  ResolutionVsEta.push_back(subDir->at(rjc_ite).make<TH3F>("ptResolutionVsEtaVsNpu","pt resolution vs eta vs npu", 50, -5., 5., 50, -0.5, 49.5, 50, -1.0, 1.0));

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
 
vector<int> 
JetAnlzr::GetSortedList(vector<float> Pt)
{

  vector<int> output;

  float max1 = 1000000.;

  for(unsigned i=0;i<Pt.size();i++){  

    float max2 = 0.;
    int max_pos = 0;

    for(unsigned ite=0; ite<Pt.size();ite++){

      if((Pt.at(ite)<max1) && (Pt.at(ite)>max2)){
        max2=Pt.at(ite);
        max_pos=ite;
      }

    }

    max1=max2;
    output.push_back(max_pos);

  }

  return output;

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
 	  
	  if(genJetPtH->at(gen_ite)<0.0) continue;

	  TLV genJet( genJetPtH->at(gen_ite),
		      genJetEtaH->at(gen_ite),
		      genJetPhiH->at(gen_ite),
		      0.0 );
	
	  genJets.push_back(genJet);

	}

	// loop over all reco jet collections
	for(unsigned rjc_ite=0; rjc_ite<input_RecoJets_.size(); rjc_ite++){

	  // get the information for the reco jets

	  Handle<vector<float> > jetPtH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"pt",jetPtH);

          vector<int> SortedPt = GetSortedList(*jetPtH);

	  Handle<vector<float> > jetEtaH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"eta",jetEtaH);

	  Handle<vector<float> > jetPhiH;
    	  iEvent.getByLabel( input_RecoJets_[rjc_ite],"phi",jetPhiH);

 	  for(unsigned reco_ite=0; reco_ite<SortedPt.size(); reco_ite++){

            int pos = SortedPt.at(reco_ite);
 	  
// 	    if((jetPtH->at(pos)<30.0) || (reco_ite==2)) break;
	    if(jetPtH->at(pos)<0.0) break; 

	    TLV recoJet( jetPtH->at(pos),
		         jetEtaH->at(pos),
			 jetPhiH->at(pos),
			 0.0 );	 

	    TLV genJet;
            if (findGenJet(recoJet,genJets,&genJet)){

              ResponseVsPt[rjc_ite]->Fill(genJet.Pt(), npu, recoJet.Pt() *1./genJet.Pt());
              ResponseVsEta[rjc_ite]->Fill(genJet.Eta(), npu, recoJet.Pt() *1./genJet.Pt());

              ResolutionVsPt[rjc_ite]->Fill(genJet.Pt(), npu, (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());
              ResolutionVsEta[rjc_ite]->Fill(genJet.Eta(), npu, (recoJet.Pt() - genJet.Pt())*1./genJet.Pt());

	    }

	  }

	}

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
