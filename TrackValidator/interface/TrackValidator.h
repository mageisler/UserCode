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
#include "MGeisler/TrackValidator/interface/TrackValidatorAlgos.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/InputTag.h"

using namespace std;
using namespace edm;
// using namespace reco;

//
// class declaration
//

class TrackValidator : public EDAnalyzer, private TrackValidatorAlgos {

  public:
    explicit TrackValidator(const ParameterSet&);
    ~TrackValidator();

    static void fillDescriptions(ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const Event&, const EventSetup&);

    virtual void beginRun(Run const&, EventSetup const&);
    virtual void endRun(Run const&, EventSetup const&);

  // ----------member data ---------------------------

  InputTag tcRefLabel_;
  vector<InputTag> tcLabels_;

  InputTag puLabel_;
  InputTag tpLabel_;

  bool ignoremissingtkcollection_;


//   string dirName_;
//   InputTag associatormap;
//   bool UseAssociators;
// 
//   bool useGsf;
//   bool runStandalone;
//   TrackingParticleSelector tpSelector;				      
//   CosmicTrackingParticleSelector cosmictpSelector;
//   MTVHistoProducerAlgo* histoProducerAlgo_;
// 
//   vector<string> associators;
//   string sim;
//   string parametersDefiner;
// 
// 
//   InputTag bsSrc;
// 
//   std::string out;
// 
//   InputTag m_dEdx1Tag;
//   InputTag m_dEdx2Tag;
// 
//   ESHandle<MagneticField> theMF;
//   vector<const TrackAssociatorBase*> associator;

};
