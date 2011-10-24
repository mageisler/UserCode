// -*- C++ -*-
//
// Package:    PU_counter
// Class:      PU_counter
// 
/**\class PU_counter PU_counter.cc MGeisler/PU_counter/src/PU_counter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Fri Oct 21 16:12:04 CEST 2011
// $Id: PU_counter.cc,v 1.1 2011/10/24 14:43:09 mgeisler Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


//
// class declaration
//
   
using namespace edm;
using namespace std;

class PU_counter : public edm::EDProducer {
   public:
      explicit PU_counter(const edm::ParameterSet&);
      ~PU_counter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

      string input_OutName_;

      TFile* outfile_MC;
      TH1F* pileup_MC;

      TFile* outfile_Data;

      // True 2011 luminosity distribution for the 6 May 2011 GoodMuon JSON file.
      TH1F* pileup_Data_May2011;
  
      // Flat10+Tail distribution taken directly from MixingModule input:  
      // (Can be used for Spring11 and Summer11 if you don't worry about small shifts in the mean) 
      // SHOULD be used for 3-D Reweighting, as this is the "true" input for all Summer11 samples.
      TH1F* pileup_Data_Flat10Tail;

      // Summer11 PU_S4, distribution obtained by only looking at the in-time crossing.  
      // This is the "spike+smear" distribution, RECOMMENDED FOR REWEIGHTING.
      TH1F* pileup_Data_SpikeSmear;

      // Summer11 PU_S4, distribution obtained by averaging the number of interactions
      // in each beam crossing to estimate the true mean.  
      TH1F* pileup_Data_PoissonInt;
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
PU_counter::PU_counter(const edm::ParameterSet& iConfig)
{
  const uint NumberOfBins=iConfig.getUntrackedParameter<int> ("NumberOfBins");
  input_OutName_ = iConfig.getUntrackedParameter<string>("OutName");

  outfile_MC = new TFile(input_OutName_.c_str(), "RECREATE");  
  outfile_Data = new TFile("MCReweightingInput.root", "RECREATE");  

  pileup_MC  = new TH1F("pileup","pileup;pileup;events",NumberOfBins,-0.5,NumberOfBins-0.5);
  pileup_Data_May2011  = new TH1F("pileup_Data_May2011","pileup;pileup;events",NumberOfBins,-0.5,NumberOfBins-0.5);
  pileup_Data_Flat10Tail  = new TH1F("pileup_Data_Flat10Tail","pileup;pileup;events",NumberOfBins,-0.5,NumberOfBins-0.5);
  pileup_Data_SpikeSmear  = new TH1F("pileup_Data_SpikeSmear","pileup;pileup;events",NumberOfBins,-0.5,NumberOfBins-0.5);
  pileup_Data_PoissonInt  = new TH1F("pileup_Data_PoissonInt","pileup;pileup;events",NumberOfBins,-0.5,NumberOfBins-0.5);

  float TrueDist2011_f[35] = {
    0.019091,
    0.0293974,
    0.0667931,
    0.108859,
    0.139533,
    0.149342,
    0.138629,
    0.114582,
    0.0859364,
    0.059324,
    0.0381123,
    0.0229881,
    0.0131129,
    0.00711764,
    0.00369635,
    0.00184543,
    0.000889604,
    0.000415683,
    0.000188921,
    0.000146288,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };

  float probdistFlat10[35] = {
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0698146584,
    0.0630151648,
    0.0526654164,
    0.0402754482,
    0.0292988928,
    0.0194384503,
    0.0122016783,
    0.007207042,
    0.004003637,
    0.0020278322,
    0.0010739954,
    0.0004595759,
    0.0002229748,
    0.0001028162,
    4.58337152809607E-05,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };

  float PoissonOneXDist[35] = {
    1.45346E-01,
    6.42802E-02,
    6.95255E-02,
    6.96747E-02,
    6.92955E-02,
    6.84997E-02,
    6.69528E-02,
    6.45515E-02,
    6.09865E-02,
    5.63323E-02,
    5.07322E-02,
    4.44681E-02,
    3.79205E-02,
    3.15131E-02,
    2.54220E-02,
    2.00184E-02,
    1.53776E-02,
    1.15387E-02,
    8.47608E-03,
    6.08715E-03,
    4.28255E-03,
    2.97185E-03,
    2.01918E-03,
    1.34490E-03,
    8.81587E-04,
    5.69954E-04,
    3.61493E-04,
    2.28692E-04,
    1.40791E-04,
    8.44606E-05,
    5.10204E-05,
    3.07802E-05,
    1.81401E-05,
    1.00201E-05,
    5.80004E-06
  };

  float PoissonIntDist[35] = {
    0.104109,
    0.0703573,
    0.0698445,
    0.0698254,
    0.0697054,
    0.0697907,
    0.0696751,
    0.0694486,
    0.0680332,
    0.0651044,
    0.0598036,
    0.0527395,
    0.0439513,
    0.0352202,
    0.0266714,
    0.019411,
    0.0133974,
    0.00898536,
    0.0057516,
    0.00351493,
    0.00212087,
    0.00122891,
    0.00070592,
    0.000384744,
    0.000219377,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
  };

  for( uint i=0; i<NumberOfBins; ++i) {
    if(i<sizeof(TrueDist2011_f))    pileup_Data_May2011->SetBinContent(i+1,TrueDist2011_f[i]);
    else  pileup_Data_May2011->SetBinContent(i+1,0.0);
    if(i<sizeof(probdistFlat10))    pileup_Data_Flat10Tail->SetBinContent(i+1,probdistFlat10[i]);
    else pileup_Data_Flat10Tail->SetBinContent(i+1,0.0);
    if(i<sizeof(PoissonOneXDist))    pileup_Data_SpikeSmear->SetBinContent(i+1,PoissonOneXDist[i]);
    else pileup_Data_SpikeSmear->SetBinContent(i+1,0.0);
    if(i<sizeof(PoissonIntDist))    pileup_Data_PoissonInt->SetBinContent(i+1,PoissonIntDist[i]);
    else  pileup_Data_PoissonInt->SetBinContent(i+1,0.0);
  }
  
}


PU_counter::~PU_counter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PU_counter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   	using namespace edm;
   	using namespace std;
  	
	//get the number of pileup vertices from the PUInfoSummary
  	edm::Handle< vector<PileupSummaryInfo> > puinfoH;
  	iEvent.getByLabel("addPileupInfo",puinfoH);
  	PileupSummaryInfo puinfo;      
	  
  	for (unsigned int puinfo_ite=0;puinfo_ite<(*puinfoH).size();++puinfo_ite){ 
    	  if ((*puinfoH)[puinfo_ite].getBunchCrossing()==0){
      	    puinfo=(*puinfoH)[puinfo_ite];
      	    break;
    	  }
  	}

	int punumber = puinfo.getPU_NumInteractions();

   	pileup_MC->Fill(punumber);
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
PU_counter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PU_counter::endJob() {

  outfile_MC->cd();
  pileup_MC->Write();
  outfile_MC->Write();
  outfile_MC->Close();

  outfile_Data->cd();
  pileup_Data_May2011->Write();
  pileup_Data_Flat10Tail->Write();
  pileup_Data_SpikeSmear->Write();
  pileup_Data_PoissonInt->Write();
  outfile_Data->Write();
  outfile_Data->Close();

}

// ------------ method called when starting to processes a run  ------------
void 
PU_counter::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
PU_counter::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
PU_counter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
PU_counter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PU_counter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PU_counter);
