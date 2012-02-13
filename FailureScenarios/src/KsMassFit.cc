// -*- C++ -*-
//
// Package:    KsMassFit
// Class:      KsMassFit
// 
/**\class KsMassFit KsMassFit.cc KsAnlz/KsMassFit/src/KsMassFit.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler
//         Created:  Tue Oct 18 16:08:19 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

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

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

//C++ include files
#include <vector>
#include <string>
#include <numeric>
#include <map>
#include <list>
#include <cctype>

//ROOT include files
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TMath.h"
#include "TStyle.h"

//
// class declaration
//

using namespace edm;
using namespace std;
using namespace reco;

class KsMassFit : public edm::EDAnalyzer {
   public:
      explicit KsMassFit(const edm::ParameterSet&);
      ~KsMassFit();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void CreateCanvas();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      string fileName_;
      TH1F* ksmasshisto;
      TH1F* lambdamasshisto;
      char thpname[256];
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
KsMassFit::KsMassFit(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed

	fileName_ = iConfig.getUntrackedParameter<string>("FileName");

    	Service<TFileService> tfserv;

	sprintf(thpname,"%s_KsMassDistribution",fileName_.c_str());
  	ksmasshisto = tfserv->make<TH1F>(thpname,"K_{s} mass distribution; K_{s} mass / MeV; # entries",100,400.,600.);

	sprintf(thpname,"%s_LambdaMassDistribution",fileName_.c_str());
  	lambdamasshisto = tfserv->make<TH1F>(thpname,"#Lambda mass distribution; #Lambda mass / GeV; # entries",100,1.05,1.2);

}


KsMassFit::~KsMassFit()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
KsMassFit::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    	Service<TFileService> tfserv;

	//get the vertex composite candidate collection for the Kshort's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
	iEvent.getByLabel("generalV0Candidates","Kshort", vertCompCandCollKshortH);

	for(VertexCompositeCandidateCollection::const_iterator iKS=vertCompCandCollKshortH->begin(); iKS!=vertCompCandCollKshortH->end(); iKS++){

	  ksmasshisto->Fill(1000.*iKS->mass());

	}

	//get the vertex composite candidate collection for the Lambda's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
	iEvent.getByLabel("generalV0Candidates","Lambda", vertCompCandCollLambdaH);

	for(VertexCompositeCandidateCollection::const_iterator iLa=vertCompCandCollLambdaH->begin(); iLa!=vertCompCandCollLambdaH->end(); iLa++){

	  lambdamasshisto->Fill(iLa->mass());

	}
}


// ------------ method called once each job just before starting event loop  ------------
void 
KsMassFit::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
KsMassFit::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
KsMassFit::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
KsMassFit::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
KsMassFit::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
KsMassFit::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
KsMassFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KsMassFit);
