// -*- C++ -*-
//
// Package:    FirstVertexAnlzr
// Class:      FirstVertexAnlzr
// 
/**\class FirstVertexAnlzr FirstVertexAnlzr.cc MGeisler/FirstVertexAnlzr/src/FirstVertexAnlzr.cc

 Description: Analyze how often the simulated first vertex is the first one in the VertexCollection and AssociationMap

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Geisler,32 4-B20,+41227676487,
//         Created:  Wed Aug 24 11:08:27 CEST 2011
// $Id$
//
//
#include "MGeisler/PF_PU_AssoMap/interface/PF_PU_AssoMap.h"

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "DataFormats/FWLite/interface/EventBase.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "CommonTools/RecoAlgos/interface/TrackingParticleSelector.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

//
// class declaration
//
   
using namespace edm;
using namespace std;
using namespace reco;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

  typedef pair<TrackRef, float> TrackQualityPair;
  typedef vector<TrackQualityPair > TrackQualityPairVector;

  typedef vector<VertexRef > VertexRefV;

class FirstVertexAnlzr : public edm::EDAnalyzer {
   public:
      explicit FirstVertexAnlzr(const edm::ParameterSet&);
      ~FirstVertexAnlzr();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

     typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

     typedef vector<pair<TrackRef, float> > TrackQualityPairVector;
     typedef pair<TrackRef, float> TrackQualityPair;


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

      TrackingParticleSelector* tpSelector;
      DQMStore* dbe_;

      LumiReWeighting LumiWeights_;

      MonitorElement *h_simulated, *h_found_T_VC_50, *h_found_T_AM_50, *h_found_T_VC_25, *h_found_T_AM_25;
      MonitorElement *h_found_V_VC, *h_found_V_AM, *h_found_V_QV;

      MonitorElement *h_num_T_TP, *h_num_T_VC, *h_num_T_AM;
      MonitorElement *h_num_V_VC, *h_num_V_AM, *h_num_V_QV;

      MonitorElement *h_vertexdistanceVC, *h_vertexdistanceAM;
      MonitorElement *h_vertexdistanceMinVC, *h_vertexdistanceMinAM;

      MonitorElement *h_simHiggs, *h_foundHiggs_VC, *h_foundHiggs_AM, *h_foundHiggs_QV;
      MonitorElement *h_numHiggs_VC, *h_numHiggs_AM, *h_numHiggs_QV;

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
FirstVertexAnlzr::FirstVertexAnlzr(const edm::ParameterSet& iConfig):LumiWeights_("Lumi_GluGluToHiggsToGG.root","Lumi_Data.root","pileup","pileup_Data_Flat10Tail")

{
   //now do what ever initialization is needed
	using namespace reco::modules;

      	dbe_ = Service<DQMStore>().operator->();

	LumiWeights_.weight3D_init();

  	ParameterSet TpSelectorPSet = iConfig.getParameter<ParameterSet>("TrackingParticleSelection");
	tpSelector = new TrackingParticleSelector(ParameterAdapter<TrackingParticleSelector>::make(TpSelectorPSet));

}


FirstVertexAnlzr::~FirstVertexAnlzr()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FirstVertexAnlzr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ 

	const edm::EventBase* iEventB = &iEvent;
  	double MyWeight3D = LumiWeights_.weight3D( (*iEventB) );
  
	//get the input vertex<->general track association map
  	Handle<TrackVertexAssMap> GTassomap;
  	iEvent.getByLabel("Tracks2Vertex",GTassomap);
	
	int am_size = GTassomap->size();

 	const VertexRef assomap_vertexref = GTassomap->begin()->key;
	double vertAM_x = assomap_vertexref->x();
	double vertAM_y = assomap_vertexref->y();
	double vertAM_z = assomap_vertexref->z();
 
	//get the input vertex collection
  	Handle<VertexCollection> vtxcoll;
  	iEvent.getByLabel("offlinePrimaryVertices",vtxcoll);
	
	int vc_size = vtxcoll->size();

	//get the reconstructed  first vertex
        const VertexRef firstvertexref(vtxcoll,0);
	double vertVC_x = firstvertexref->x();
	double vertVC_y = firstvertexref->y();
	double vertVC_z = firstvertexref->z();

	//get the vertex composite candidate collection for the Kshort's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH;
	iEvent.getByLabel("generalV0Candidates","Kshort", vertCompCandCollKshortH);

	//get the vertex composite candidate collection for the Lambda's
	Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH;
	iEvent.getByLabel("generalV0Candidates","Lambda", vertCompCandCollLambdaH);

	//get the displaced vertex collection for nuclear interactions
	Handle<PFDisplacedVertexCollection> displVertexCollH;
	iEvent.getByLabel("particleFlowDisplacedVertex", displVertexCollH);
 
	//get the qualified vertex collection
 	VertexRefV* qualvtxcoll = PF_PU_AssoMapAlgos::QualifiedVertices(vtxcoll,true,4.,displVertexCollH,vertCompCandCollKshortH,vertCompCandCollLambdaH);
	
	int qv_size = qualvtxcoll->size();
  	
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

	//get the Tracking Particles
        Handle<TrackingParticleCollection>  TPCollectionH;
        iEvent.getByLabel("mergedtruth","MergedTrackTruth",TPCollectionH);

	//get the Tracking Vertices
        Handle<TrackingVertexCollection>  TVCollectionH;
        iEvent.getByLabel("mergedtruth","MergedTrackTruth",TVCollectionH);

	int anzahl = 0;
	vector <TrackingVertex>* TV_sim(new vector <TrackingVertex>() );

	for (TrackingVertexCollection::const_iterator v = TVCollectionH -> begin(); v != TVCollectionH -> end(); ++v) {
	  if ((v->eventId().bunchCrossing()==0) && (v->eventId().event()>anzahl)){	
// 	    cout << v->eventId().event() << endl;
	    anzahl=v->eventId().event();
	    TV_sim->push_back(*(v));
	  }
	}

// 	cout << punumber << " = " << TV_sim->size() << endl;

	//get the main interaction vertex  
	Handle<SimVertexContainer> simVtxs;
  	iEvent.getByLabel("g4SimHits", simVtxs);

	//get simulated first vertex
	SimVertexContainer::const_iterator firstsimvtx=simVtxs->begin();
	SimVertex maininteraction = simVtxs->at(0);
	double firstvertx_x = firstsimvtx->position().x();
   	double firstvertx_y = firstsimvtx->position().y();
	double firstvertx_z = firstsimvtx->position().z(); 

	double MinDistanceAM = 100000.;

	for(TrackVertexAssMap::const_iterator am_ite=GTassomap->begin(); am_ite!=GTassomap->end();am_ite++){

	  VertexRef assomap_vertexref_tmp = am_ite->key;
	  double vertAM_x_tmp = assomap_vertexref_tmp->x();
	  double vertAM_y_tmp = assomap_vertexref_tmp->y();
	  double vertAM_z_tmp = assomap_vertexref_tmp->z();

	  double distance = (firstvertx_x-vertAM_x_tmp)*(firstvertx_x-vertAM_x_tmp) + (firstvertx_y-vertAM_y_tmp)*(firstvertx_y-vertAM_y_tmp) + (firstvertx_z-vertAM_z_tmp)*(firstvertx_z-vertAM_z_tmp);

	  if(distance<MinDistanceAM) MinDistanceAM=distance;	

	}

	double MinDistanceVC = 100000.;

	for(VertexCollection::const_iterator vc_ite=vtxcoll->begin(); vc_ite!=vtxcoll->end();vc_ite++){

	  double vertVC_x_tmp = vc_ite->x();
	  double vertVC_y_tmp = vc_ite->y();
	  double vertVC_z_tmp = vc_ite->z();

	  double distance = (firstvertx_x-vertVC_x_tmp)*(firstvertx_x-vertVC_x_tmp) + (firstvertx_y-vertVC_y_tmp)*(firstvertx_y-vertVC_y_tmp) + 	(firstvertx_z-vertVC_z_tmp)*(firstvertx_z-vertVC_z_tmp);

	  if(distance<MinDistanceVC) MinDistanceVC=distance;	

	}

	//get the input VC track collection
  	Handle<View<Track> >  VCtrks;
  	iEvent.getByLabel("FirstVertexTrackCollection","VCTracks",VCtrks);
  	Handle<TrackCollection >  VCtrksH;
  	iEvent.getByLabel("FirstVertexTrackCollection","VCTracks",VCtrksH);
 
	//get the input AM track collection
  	Handle<View<Track> >  AMtrks;
  	iEvent.getByLabel("FirstVertexTrackCollection",AMtrks);
  	Handle<TrackCollection >  AMtrksH;
  	iEvent.getByLabel("FirstVertexTrackCollection",AMtrksH);


	//associate reco tracks to tracking particles
	ESHandle<TrackAssociatorBase> theAssociator;
	iSetup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theAssociator);
	TrackAssociatorBase* theTrackAssociator_ = (TrackAssociatorBase *) theAssociator.product(); 

        RecoToSimCollection recSimCollVC;
	recSimCollVC=theTrackAssociator_->associateRecoToSim(VCtrks,TPCollectionH,&iEvent,&iSetup);

        RecoToSimCollection recSimCollAM;
	recSimCollAM=theTrackAssociator_->associateRecoToSim(AMtrks,TPCollectionH,&iEvent,&iSetup);

	int simTracks = 0;
	double higgs_x = firstvertx_x;
	double higgs_y = firstvertx_y;
	double higgs_z = firstvertx_z;
	bool foundHiggs = false;

        for (TrackingParticleCollection::size_type i=0; i<TPCollectionH->size(); i++){ 
	  
 	  TrackingParticleRef tpr(TPCollectionH, i);
	  TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());

	  if((*tpSelector)(*tp)) simTracks++;	  

	  //find the Higgs vertex
	  if((tp->eventId().event() == 0) && (tp->eventId().bunchCrossing() == 0) && ((tp->pdgId()==25) || (tp->pdgId()==37))){

	    higgs_x = tp->vertex().x();
	    higgs_y = tp->vertex().y();
	    higgs_z = tp->vertex().z();
	    foundHiggs=true;
	    break;    

	  }

	}

	int foundAM = 0;
	int foundVC = 0;

	int recoTracksVC = VCtrksH->size();
	int recoTracksAM = AMtrksH->size();
	
	for(unsigned int am_ite=0; am_ite<AMtrksH->size();am_ite++){

	  const TrackBaseRef am_trackbaseref(AMtrks,am_ite);
	  
	  if(recSimCollAM.find(am_trackbaseref) == recSimCollAM.end()) continue;
// 	  cout << "YEAH 1a" << endl;

	  vector<pair<TrackingParticleRef, double> > tp;
	  tp = recSimCollAM[am_trackbaseref];

	  if(tp.size()==0) continue;
	  
	  for(unsigned int tp_ite=0; tp_ite<tp.size();tp_ite++){

	    if ((*tpSelector)(*(tp.at(tp_ite).first))){
	      foundAM++;
// 	      cout << "YEAH 2a" << endl;
	      break;
	    }

	  }

	}
	
	for(unsigned int vc_ite=0; vc_ite<VCtrksH->size();vc_ite++){

	  const TrackBaseRef vc_trackbaseref(VCtrks,vc_ite);
	  
	  if(recSimCollVC.find(vc_trackbaseref) == recSimCollVC.end()) continue;
// 	  cout << "YEAH 1b" << endl;

	  vector<pair<TrackingParticleRef, double> > tp;
	  tp = recSimCollVC[vc_trackbaseref];

	  if(tp.size()==0) continue;
	  
	  for(unsigned int tp_ite=0; tp_ite<tp.size();tp_ite++){

	    if ((*tpSelector)(*(tp.at(tp_ite).first))){
	      foundVC++;
// 	      cout << "YEAH 2b" << endl;
	      break;
	    }

	  }

	}

	int foundVC_V = 0;
	int foundAM_V = 0;
	int foundQV_V = 0;
 
	for(unsigned int tv_ite=0; tv_ite<TV_sim->size();tv_ite++){

	  double tv_x = TV_sim->at(tv_ite).position().x();
	  double tv_y = TV_sim->at(tv_ite).position().y();
	  double tv_z = TV_sim->at(tv_ite).position().z();

	  for(VertexCollection::const_iterator vtx_ite=vtxcoll->begin(); vtx_ite!=vtxcoll->end(); vtx_ite++){
	    if((fabs(vtx_ite->position().x()-tv_x)<=0.5) && (fabs(vtx_ite->position().y()-tv_y)<=0.5) && (fabs(vtx_ite->position().z()-tv_z)<=0.5)){
 	      foundVC_V++;
	      break;
	    }
	  }

	  for(TrackVertexAssMap::const_iterator assomap_ite=GTassomap->begin(); assomap_ite!=GTassomap->end(); assomap_ite++){
	    VertexRef vtx_ite = assomap_ite->key;
	    if((fabs(vtx_ite->position().x()-tv_x)<=0.5) && (fabs(vtx_ite->position().y()-tv_y)<=0.5) && (fabs(vtx_ite->position().z()-tv_z)<=0.5)){
              foundAM_V++; 
	      break;
	    }
	  }

	  for(unsigned ite=0; ite<qualvtxcoll->size(); ite++){
	    VertexRef vtx_ite = qualvtxcoll->at(ite);
	    if((fabs(vtx_ite->position().x()-tv_x)<=0.5) && (fabs(vtx_ite->position().y()-tv_y)<=0.5) && (fabs(vtx_ite->position().z()-tv_z)<=0.5)){
	      foundQV_V++;
	      break;
	    }
	  }

	}

	if(simTracks>0.) h_simulated->Fill(punumber,MyWeight3D);
	if(foundVC >= (recoTracksVC*1./2.)) h_found_T_VC_50->Fill(punumber,MyWeight3D);
	if(foundAM >= (recoTracksAM*1./2.)) h_found_T_AM_50->Fill(punumber,MyWeight3D);

	if(foundVC >= (recoTracksVC*1./4.)) h_found_T_VC_25->Fill(punumber,MyWeight3D);
	if(foundAM >= (recoTracksAM*1./4.)) h_found_T_AM_25->Fill(punumber,MyWeight3D);

	h_num_T_TP->Fill(punumber,simTracks,MyWeight3D);
	h_num_T_VC->Fill(punumber,recoTracksVC,MyWeight3D);
	h_num_T_AM->Fill(punumber,recoTracksAM,MyWeight3D);
	
	double distanceVC = (firstvertx_x-vertVC_x)*(firstvertx_x-vertVC_x) + (firstvertx_y-vertVC_y)*(firstvertx_y-vertVC_y) + (firstvertx_z-vertVC_z)*(firstvertx_z-vertVC_z);
	double distanceAM = (firstvertx_x-vertAM_x)*(firstvertx_x-vertAM_x) + (firstvertx_y-vertAM_y)*(firstvertx_y-vertAM_y) + (firstvertx_z-vertAM_z)*(firstvertx_z-vertAM_z);

	h_vertexdistanceVC->Fill(punumber,sqrt(distanceVC),MyWeight3D);
	h_vertexdistanceAM->Fill(punumber,sqrt(distanceAM),MyWeight3D);

	if(distanceVC==MinDistanceVC) h_vertexdistanceMinVC->Fill(punumber,MyWeight3D);
	if(distanceAM==MinDistanceAM) h_vertexdistanceMinAM->Fill(punumber,MyWeight3D);

	h_num_V_VC->Fill(punumber,vc_size,MyWeight3D);
 	h_num_V_AM->Fill(punumber,am_size,MyWeight3D);
	h_num_V_QV->Fill(punumber,qv_size,MyWeight3D);

	h_found_V_VC->Fill(punumber,foundVC_V,MyWeight3D);
	h_found_V_AM->Fill(punumber,foundAM_V,MyWeight3D);
	h_found_V_QV->Fill(punumber,foundQV_V,MyWeight3D);

	int Hpos_VC = 0;
	for(VertexCollection::const_iterator vtx_ite=vtxcoll->begin(); vtx_ite!=vtxcoll->end(); vtx_ite++,Hpos_VC++){
  	  double x_dist = vtx_ite->position().x() - higgs_x;
	  double y_dist = vtx_ite->position().y() - higgs_y;
	  double z_dist = vtx_ite->position().z() - higgs_z;
	  double Higgs_distance = x_dist*x_dist + y_dist*y_dist + z_dist*z_dist; 
	  if(sqrt(Higgs_distance)<=0.05){
   	    h_foundHiggs_VC->Fill(punumber,MyWeight3D);
	    h_numHiggs_VC->Fill(punumber,Hpos_VC,MyWeight3D);
	    break;
	  }
	}

	int Hpos_AM = 0;
	for(TrackVertexAssMap::const_iterator assomap_ite=GTassomap->begin(); assomap_ite!=GTassomap->end(); assomap_ite++,Hpos_AM++){
 	  VertexRef vtx_ite = assomap_ite->key;
  	  double x_dist = vtx_ite->position().x() - higgs_x;
	  double y_dist = vtx_ite->position().y() - higgs_y;
	  double z_dist = vtx_ite->position().z() - higgs_z;
	  double Higgs_distance = x_dist*x_dist + y_dist*y_dist + z_dist*z_dist; 
	  if(sqrt(Higgs_distance)<=0.05){
   	    h_foundHiggs_AM->Fill(punumber,MyWeight3D);
	    h_numHiggs_AM->Fill(punumber,Hpos_AM,MyWeight3D);
	    break;
	  }
	}

	int Hpos_QV = 0;
	for(unsigned ite=0; ite<qualvtxcoll->size(); ite++,Hpos_QV++){
	  VertexRef vtx_ite = qualvtxcoll->at(ite);
  	  double x_dist = vtx_ite->position().x() - higgs_x;
	  double y_dist = vtx_ite->position().y() - higgs_y;
	  double z_dist = vtx_ite->position().z() - higgs_z;
	  double Higgs_distance = x_dist*x_dist + y_dist*y_dist + z_dist*z_dist; 
	  if(sqrt(Higgs_distance)<=0.05){
   	    h_foundHiggs_QV->Fill(punumber,MyWeight3D);
	    h_numHiggs_QV->Fill(punumber,Hpos_QV,MyWeight3D);
	    break;
	  }
	}

// 	cout << "Simulated First Vertex Tracks: " << simTracks << endl;
// 	cout << "Found First Vertex Tracks in VC: " << foundVC << " out of " << recoTracksVC << endl;
// 	cout << "Found First Vertex Tracks in AM: " << foundAM << " out of " << recoTracksAM << endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
FirstVertexAnlzr::beginJob()
{

	h_simulated = dbe_->book1D("h_simulated","# events with a simulated first vertex; number of pileup vertices",51,-0.5,50.5);

	h_found_T_VC_50 = dbe_->book1D("h_found_T_VC_50","# events with found first vertex in VertexCollection (50%); number of pileup vertices",51,-0.5,50.5);
	h_found_T_AM_50 = dbe_->book1D("h_found_T_AM_50","# events with found first vertex in AssociationMap (50%); number of pileup vertices",51,-0.5,50.5);

	h_found_T_VC_25 = dbe_->book1D("h_found_T_VC_25","# events with found first vertex in VertexCollection (25%); number of pileup vertices",51,-0.5,50.5);
	h_found_T_AM_25 = dbe_->book1D("h_found_T_AM_25","# events with found first vertex in AssociationMap (25%); number of pileup vertices",51,-0.5,50.5);

	h_found_V_VC = dbe_->book2D("h_found_V_VC","number of found Tracking Vertices vs number of pileup vertices; number of pileup vertices; number of found Tracking Vertices in the VertexCollection",51,-0.5,50.5,51,-0.5,50.5);
	h_found_V_AM = dbe_->book2D("h_found_V_AM","number of found Tracking Vertices vs number of pileup vertices; number of pileup vertices; number of found Tracking Vertices in the AssociationMap",51,-0.5,50.5,51,-0.5,50.5);
	h_found_V_QV = dbe_->book2D("h_found_V_QV","number of found Tracking Vertices vs number of pileup vertices; number of pileup vertices; number of found Tracking Vertices in the qualified VertexCollection",51,-0.5,50.5,51,-0.5,50.5);

	h_num_T_TP = dbe_->book2D("h_num_T_TP","number of Tracking Particles vs number of pileup vertices; number of pileup vertices; number of Tracking Particles",51,-0.5,50.5,501,-0.5,500.5);
	h_num_T_VC = dbe_->book2D("h_num_T_VC","number of VC Tracks vs number of pileup vertices; number of pileup vertices; number of tracks from the Vertex Collection",51,-0.5,50.5,251,-0.5,250.5);
	h_num_T_AM = dbe_->book2D("h_num_T_AM","number of AM Tracks vs number of pileup vertices; number of pileup vertices; number of tracks from the Association Map",51,-0.5,50.5,251,-0.5,250.5);

	h_num_V_VC = dbe_->book2D("h_num_V_VC","number of VC Vertices vs number of pileup vertices; number of pileup vertices; number of vertices from the Vertex Collection",51,-0.5,50.5,51,-0.5,50.5);
	h_num_V_AM = dbe_->book2D("h_num_V_AM","number of AM Vertices vs number of pileup vertices; number of pileup vertices; number of vertices from the Association Map",51,-0.5,50.5,51,-0.5,50.5);
	h_num_V_QV = dbe_->book2D("h_num_V_QV","number of qualified AM Vertices vs number of pileup vertices; number of pileup vertices; number of qualified VertexCollection",51,-0.5,50.5,51,-0.5,50.5);

	h_vertexdistanceVC = dbe_->book2D("h_vertexdistanceVC","3D distance between main interaction and reco first vertex; number of pileup vertices; distance",51,-0.5,50.5,101,-0.05,10.05);
	h_vertexdistanceAM = dbe_->book2D("h_vertexdistanceAM","3D distance between main interaction and reco first vertex; number of pileup vertices; distance",51,-0.5,50.5,101,-0.05,10.05);

	h_vertexdistanceMinVC = dbe_->book1D("h_vertexdistanceMinVC","# events with found first vertex in VertexCollection (MinDist); number of pileup vertices",51,-0.5,50.5);
	h_vertexdistanceMinAM = dbe_->book1D("h_vertexdistanceMinAM","# events with found first vertex in AssociationMap (MinDist); number of pileup vertices",51,-0.5,50.5);

	h_simHiggs = dbe_->book1D("h_simHiggs","# events with a simulated Higgs; number of pileup vertices",51,-0.5,50.5);

	h_foundHiggs_VC = dbe_->book1D("h_foundHiggs_VC","# events with found Higgs vertex in VertexCollection; number of pileup vertices",51,-0.5,50.5);
	h_foundHiggs_AM = dbe_->book1D("h_foundHiggs_AM","# events with found Higgs vertex in AssociationMap; number of pileup vertices",51,-0.5,50.5);
	h_foundHiggs_QV = dbe_->book1D("h_foundHiggs_QV","# events with found Higgs vertex in qualified VertexCollection; number of pileup vertices",51,-0.5,50.5);

	h_numHiggs_VC = dbe_->book2D("h_numHiggs_VC","position of the main interaction vs number of pileup vertices; number of pileup vertices; position of the main interaction in the Vertex Collection",51,-0.5,50.5,51,-0.5,50.5);
	h_numHiggs_AM = dbe_->book2D("h_numHiggs_AM","position of the main interaction vs number of pileup vertices; number of pileup vertices; position of the main interaction in the Association Map",51,-0.5,50.5,51,-0.5,50.5);
	h_numHiggs_QV = dbe_->book2D("h_numHiggs_QV","position of the main interaction vs number of pileup vertices; number of pileup vertices; position of the main interaction in the qualified VertexCollection",51,-0.5,50.5,51,-0.5,50.5);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
FirstVertexAnlzr::endJob() 
{

	dbe_->save("FirstVertexAnalysis.root");

}

// ------------ method called when starting to processes a run  ------------
void 
FirstVertexAnlzr::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
FirstVertexAnlzr::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
FirstVertexAnlzr::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
FirstVertexAnlzr::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FirstVertexAnlzr::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FirstVertexAnlzr);
