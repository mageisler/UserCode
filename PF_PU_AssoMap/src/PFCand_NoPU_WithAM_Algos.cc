#include "MGeisler/PF_PU_AssoMap/interface/PFCand_NoPU_WithAM_Algos.h"

#include <vector>
#include <string>

#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OneToManyWithQuality.h"
#include "DataFormats/Common/interface/OneToManyWithQualityGeneric.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"

#include "TMath.h"
   
using namespace edm;
using namespace std;
using namespace reco;

const double eMass = 0.000511;
const double kMass = 0.497;
const double piMass = 0.1396;

  typedef AssociationMap<OneToManyWithQuality< VertexCollection, TrackCollection, float> > TrackVertexAssMap;

  typedef vector<VertexRef > VertexRefV;

  typedef math::XYZTLorentzVector LorentzVector;


/*******************************************************************************************/
/* function to create a vertex collection containing only the vertices that pass the cuts  */ 
/*******************************************************************************************/

VertexRefV* 
PFCand_NoPU_WithAM_Algos::QualifiedVertices(Handle<VertexCollection> input_vtxcoll, bool input_VertexQuality_, double input_VertexMinNdof_) 
{

	VertexRefV* qualvtxcoll(new VertexRefV() );

  	for(unsigned int index_vtx=0;  index_vtx<input_vtxcoll->size(); ++index_vtx){

          const VertexRef vertexref(input_vtxcoll,index_vtx);

	  if(!(input_VertexQuality_) ||  
	     ((vertexref->ndof()>=input_VertexMinNdof_) 
	      && !(vertexref->isFake())
	      && (vertexref->position().z()<24.))) qualvtxcoll->push_back(vertexref);

	}

	return qualvtxcoll;

}


/*************************************************************************************/
/* function to find the closest vertex in z for a certain point                      */ 
/*************************************************************************************/

VertexRef 
PFCand_NoPU_WithAM_Algos::FindClosestInZ(double ztrack, VertexRefV* qualvtxcoll)
{

	VertexRef bestvertexref;

	double dzmin = 5.;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<qualvtxcoll->size(); ++index_vtx){

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
 
	  //find and store the closest vertex in z
          double dz = fabs(ztrack - vertexref->z());
          if(dz<dzmin) {
            dzmin = dz; 
            bestvertexref = vertexref;
          }
	
	}

	if(dzmin<5.) return bestvertexref;
	  else return qualvtxcoll->at(0);
}


/*************************************************************************************/
/* function to compare two pfcandidates                                              */ 
/*************************************************************************************/

bool
PFCand_NoPU_WithAM_Algos::Match(const PFCandidatePtr pfc, const RecoCandidate* rc)
{

	return (
	  (fabs(pfc->eta()-rc->eta())<0.1) &&
	  (fabs(pfc->phi()-rc->phi())<0.1) &&
	  (fabs(pfc->vertexChi2()-rc->vertexChi2())<0.1) &&
	  (fabs(pfc->vertexNdof()-rc->vertexNdof())<0.1) &&
	  (fabs(pfc->p()-rc->p())<0.1) &&
	  (pfc->charge() == rc->charge())
	);
}


/*************************************************************************************/
/* function to find out if the track comes from a gamma conversion                   */ 
/*************************************************************************************/

bool
PFCand_NoPU_WithAM_Algos::ComesFromConversion(const PFCandidatePtr candptr, Handle<ConversionCollection> convCollH, 
					      VertexRefV* qualvtxcoll, VertexRef* primVtxRef)
{

	Conversion* gamma = new Conversion();

	if(candptr->gsfElectronRef().isNull()) return false;

	for(unsigned int convcoll_ite=0; convcoll_ite<convCollH->size(); convcoll_ite++){
	
	  if(ConversionTools::matchesConversion(*(candptr->gsfElectronRef()),convCollH->at(convcoll_ite))){
	
	    *gamma = convCollH->at(convcoll_ite);

	    double ztrackfirst = gamma->conversionVertex().z();
	    double radius = gamma->conversionVertex().position().rho();
	    double tracktheta = candptr->theta();
	    if(gamma->nTracks()==2) tracktheta = gamma->pairMomentum().theta();

	    double ztrack = ztrackfirst - (radius/tan(tracktheta));

	    *primVtxRef = FindClosestInZ(ztrack,qualvtxcoll);

	  }

	}

	return false;
}


/*************************************************************************************/
/* function to find the best vertex for a pf candidate                               */ 
/*************************************************************************************/

VertexRef
PFCand_NoPU_WithAM_Algos::FindPFCandVertex(const PFCandidatePtr candptr, VertexRefV* qualvtxcoll, const edm::EventSetup& iSetup)
{

	VertexRef bestvertexref = qualvtxcoll->at(qualvtxcoll->size()-1);

	ESHandle<TransientTrackBuilder> theTTB;
  	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

  	ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	GlobalVector trkMomentum = GlobalVector( candptr->momentum().x(),
		          		         candptr->momentum().y(),
				     	   	 candptr->momentum().z() );

	const math::XYZPoint interactionPoint = candptr->vertex();
	const math::XYZVector primaryVector( trkMomentum.x(),trkMomentum.y(),trkMomentum.z() );
	const TrackBase::CovarianceMatrix covMat;

    	float chi = 0.;				      
    	float ndf = 0.;

	Track primary = Track(chi,ndf,interactionPoint,primaryVector,candptr->charge(),covMat);
	TransientTrack genTrkTT(primary, &(*bFieldHandle) );

  	double IpMin = 3.;
	bool VertexAssUseAbsDistance = true;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<qualvtxcoll->size(); ++index_vtx){

          VertexRef vertexref = qualvtxcoll->at(index_vtx);
	        
	  double genTrk3DIp = 10001.;
	  double genTrk3DIpSig = 10001.;
	  pair<bool,Measurement1D> genTrk3DIpPair = IPTools::absoluteImpactParameter3D(genTrkTT, *vertexref);

	  if(genTrk3DIpPair.first){
    	    genTrk3DIp = genTrk3DIpPair.second.value();
	    genTrk3DIpSig = genTrk3DIpPair.second.significance();
	  }
 
	  //find and store the closest vertex
	  if(VertexAssUseAbsDistance){
            if(genTrk3DIp<IpMin){
              IpMin = genTrk3DIp; 
              bestvertexref = vertexref;
            }
	  }else{
            if(genTrk3DIpSig<IpMin){
              IpMin = genTrk3DIpSig; 
              bestvertexref = vertexref;
            }
	  }

        }

	return bestvertexref;
}


/*************************************************************************************/
/* function to find out if the track comes from a V0 decay                           */ 
/*************************************************************************************/

bool
PFCand_NoPU_WithAM_Algos::ComesFromV0Decay(const PFCandidatePtr candptr, Handle<VertexCompositeCandidateCollection> vertCompCandCollKshortH, 
	 	 	  	           Handle<VertexCompositeCandidateCollection> vertCompCandCollLambdaH, VertexRefV* qualvtxcoll, VertexRef* primVtxRef)
{

	//the part for the reassociation of particles from Kshort decays
	for(VertexCompositeCandidateCollection::const_iterator iKS=vertCompCandCollKshortH->begin(); iKS!=vertCompCandCollKshortH->end(); iKS++){

	  const RecoCandidate *dauCand1 = dynamic_cast<const RecoCandidate*>(iKS->daughter(0));
	  const RecoCandidate *dauCand2 = dynamic_cast<const RecoCandidate*>(iKS->daughter(1));

	  if(PFCand_NoPU_WithAM_Algos::Match(candptr,dauCand1) || PFCand_NoPU_WithAM_Algos::Match(candptr,dauCand2)){

            double ztrackfirst = iKS->vertex().z();
	    double radius = iKS->vertex().rho();
	    double tracktheta = iKS->p4().theta();

	    double ztrack = ztrackfirst - (radius/tan(tracktheta));

     	    *primVtxRef = FindClosestInZ(ztrack,qualvtxcoll);

	    return true;

	  }

	}

	//the part for the reassociation of particles from Lambda decays
	for(VertexCompositeCandidateCollection::const_iterator iLambda=vertCompCandCollLambdaH->begin(); iLambda!=vertCompCandCollLambdaH->end(); iLambda++){

	  const RecoCandidate *dauCand1 = dynamic_cast<const RecoCandidate*>(iLambda->daughter(0));
	  const RecoCandidate *dauCand2 = dynamic_cast<const RecoCandidate*>(iLambda->daughter(1));

	  if(PFCand_NoPU_WithAM_Algos::Match(candptr,dauCand1) || PFCand_NoPU_WithAM_Algos::Match(candptr,dauCand2)){

            double ztrackfirst = iLambda->vertex().z();
	    double radius = iLambda->vertex().rho();
	    double tracktheta = iLambda->p4().theta();

	    double ztrack = ztrackfirst - (radius/tan(tracktheta));

     	    *primVtxRef = FindClosestInZ(ztrack,qualvtxcoll);

	    return true;

	  }

	}

	return false;
}


/*************************************************************************************/
/* function to find the closest vertex in z for a track from a nuclear interaction   */ 
/*************************************************************************************/

VertexRef
PFCand_NoPU_WithAM_Algos::FindNIVertex(const PFCandidatePtr candptr, PFDisplacedVertex displVtx, VertexRefV* qualvtxcoll, 
				       bool oneDim, const edm::EventSetup& iSetup)
{

	TrackCollection refittedTracks = displVtx.refittedTracks();

	if((displVtx.isTherePrimaryTracks()) || (displVtx.isThereMergedTracks())){

	  for(TrackCollection::const_iterator trkcoll_ite=refittedTracks.begin(); trkcoll_ite!=refittedTracks.end(); trkcoll_ite++){
	
	    const TrackBaseRef retrackbaseref = displVtx.originalTrack(*trkcoll_ite); 

	    if(displVtx.isIncomingTrack(retrackbaseref)){

	      double ztrackfirst = candptr->vertex().z();
	      double radius = candptr->vertex().rho();
	      double tracktheta = candptr->theta();

              double ztrack = ztrackfirst - (radius/tan(tracktheta));

	      return FindClosestInZ(ztrack,qualvtxcoll);	      

	    }

	  }

	}
	
	math::XYZTLorentzVector mom_sec = displVtx.secondaryMomentum((string) "PI", true);

        double ztrackfirst = displVtx.position().z();
	double radius = displVtx.position().rho();     
	double tracktheta = mom_sec.theta();	

	double ztrack = ztrackfirst - (radius/tan(tracktheta));

	if(oneDim) return FindClosestInZ(ztrack,qualvtxcoll);
	else{

	  VertexRef bestvertexref = qualvtxcoll->at(0);

	  ESHandle<TransientTrackBuilder> theTTB;
  	  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

	  ESHandle<MagneticField> bFieldHandle;
	  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	  GlobalVector trkMomentum = GlobalVector( mom_sec.x(),
					           mom_sec.y(),
				     	   	   mom_sec.z() );

	  const math::XYZPoint interactionPoint = displVtx.position();
	  const math::XYZVector primaryVector( trkMomentum.x(),trkMomentum.y(),trkMomentum.z() );
	  const TrackBase::CovarianceMatrix covMat;

	  Track primary = Track(displVtx.chi2(),displVtx.ndof(),interactionPoint,primaryVector,0,covMat);
	  TransientTrack genPhoTT(primary, &(*bFieldHandle) );
    
          KinematicParticleFactoryFromTransientTrack pFactory;
	  vector<RefCountedKinematicParticle> PhoParticles;

    	  float chi = 0.;				      
    	  float ndf = 0.;
	  float piMassSigma = piMass*1.e-6;

	  //loop over all tracks from the nuclear interaction	
	  for(unsigned convtrk_ite=0; convtrk_ite<refittedTracks.size(); convtrk_ite++){

    	    // get kinematic particles          		   
    	    TransientTrack ConvDau = (*theTTB).build(refittedTracks.at(convtrk_ite)); 
            PhoParticles.push_back(pFactory.particle(ConvDau, piMass, chi, ndf, piMassSigma));

	  }

	  KinematicParticleVertexFitter fitter;    

	  RefCountedKinematicTree PhoVertexFitTree = fitter.fit(PhoParticles);

	  if(PhoVertexFitTree->isValid()){

      	    PhoVertexFitTree->movePointerToTheTop();						       
            RefCountedKinematicParticle PhoFitKinematicParticle = PhoVertexFitTree->currentParticle();

	    KinematicState theCurrentKinematicState = PhoFitKinematicParticle->currentState();
            FreeTrajectoryState thePhoFTS = theCurrentKinematicState.freeTrajectoryState();
            genPhoTT = (*theTTB).build(thePhoFTS);

	  }

  	  double IpMin = 10.;
          
	  //loop over all vertices with a good quality in the vertex collection
  	  for(unsigned int index_vtx=0;  index_vtx<qualvtxcoll->size(); ++index_vtx){

            VertexRef vertexref = qualvtxcoll->at(index_vtx);
	        
	    double genPho3DIpSig = 10001.;
	    pair<bool,Measurement1D> genPho3DIpPair = IPTools::signedImpactParameter3D(genPhoTT, trkMomentum, *vertexref);

	    if(genPho3DIpPair.first){
	      genPho3DIpSig = fabs(genPho3DIpPair.second.significance());
	    }
 
	    //find and store the closest vertex
            if(genPho3DIpSig<IpMin){
              IpMin = genPho3DIpSig; 
              bestvertexref = vertexref;
            }
          }

	return bestvertexref;

        } 

}


/*************************************************************************************/
/* function to find out if the track comes from a nuclear interaction                */ 
/*************************************************************************************/

bool
PFCand_NoPU_WithAM_Algos::ComesFromNI(const PFCandidatePtr candptr, Handle<PFDisplacedVertexCollection> displVertexCollH, 
				      PFDisplacedVertex* displVtx, const edm::EventSetup& iSetup)
{

	ESHandle<TransientTrackBuilder> theTTB;
  	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

	ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	GlobalVector trkMomentum = GlobalVector( candptr->momentum().x(),
		          		         candptr->momentum().y(),
				     	   	 candptr->momentum().z() );

	const math::XYZPoint interactionPoint = candptr->vertex();
	const math::XYZVector primaryVector( trkMomentum.x(),trkMomentum.y(),trkMomentum.z() );
	const TrackBase::CovarianceMatrix covMat;

    	float chi = 0.;				      
    	float ndf = 0.;

	Track primary = Track(chi,ndf,interactionPoint,primaryVector,candptr->charge(),covMat);
	TransientTrack genTrkTT(primary, &(*bFieldHandle) );

  	double IpMin = 100.;
	bool VertexAssUseAbsDistance = true;
          
	//loop over all vertices with a good quality in the vertex collection
  	for(unsigned int index_vtx=0;  index_vtx<displVertexCollH->size(); ++index_vtx){

          PFDisplacedVertexRef pfvertexref(displVertexCollH,index_vtx);
	        
	  double genTrk3DIp = 10001.;
	  double genTrk3DIpSig = 10001.;
	  pair<bool,Measurement1D> genTrk3DIpPair = IPTools::absoluteImpactParameter3D(genTrkTT, *pfvertexref);

	  if(genTrk3DIpPair.first){
    	    genTrk3DIp = genTrk3DIpPair.second.value();
	    genTrk3DIpSig = genTrk3DIpPair.second.significance();
	  }
 
	  //find and store the closest vertex
	  if(VertexAssUseAbsDistance){
            if(genTrk3DIp<IpMin){
              IpMin = genTrk3DIp; 
              *displVtx = *pfvertexref;
            }
	  }else{
            if(genTrk3DIpSig<IpMin){
              IpMin = genTrk3DIpSig; 
              *displVtx = *pfvertexref;
            }
	  }

        }

	if(IpMin<0.5) return true;

	return false;
}