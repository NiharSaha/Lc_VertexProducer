// -*- C++ -*-

// Package:    VertexCompositeProducer
// Class:      LamC3PFitter

/**\class LamC3PFitter LamC3PFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/LamC3PFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/


#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//Nihar
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
//#include "RecoVertex/KinematicFit/interface/RefCountedKinematicTree.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//#include "DataFormats/VertexReco/interface/VertexExtra.h"
//#include "DataFormats/VertexReco/interface/Error4D.h"



#include <vector>
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

const float piMassLamC3P = 0.13957018;
const float piMassLamC3PSquared = piMassLamC3P*piMassLamC3P;
const float kaonMassLamC3P = 0.493677;
const float kaonMassLamC3PSquared = kaonMassLamC3P*kaonMassLamC3P;
const float protonMassLamC3P = 0.938272013; 
const float protonMassLamC3PSquared = protonMassLamC3P*protonMassLamC3P;
const float lamCMassLamC3P = 2.28646;
float piMassLamC3P_sigma = 3.5E-7f;
float kaonMassLamC3P_sigma = 1.6E-5f;
float protonMassLamC3P_sigma = 1.6E-5f;
float lamCMassLamC3P_sigma = lamCMassLamC3P*1.e-6;

float cand1Mass[2] = {piMassLamC3P, protonMassLamC3P};
float cand2Mass[2] = {protonMassLamC3P, piMassLamC3P};
float cand1Mass_sigma[2] = {piMassLamC3P_sigma, protonMassLamC3P_sigma};
float cand2Mass_sigma[2] = {protonMassLamC3P_sigma, piMassLamC3P_sigma};



#define PROTON_MASS 0.9383
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677

using std::vector;
//using std::cout;
//using std::endl;
using std::string;
using namespace reco;
using namespace edm;
using namespace std;







// Constructor and (empty) destructor
//LamC3PFitter::LamC3PFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC) :
//    bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>())
LamC3PFitter::LamC3PFitter(const edm::ParameterSet& theParameters, edm::ConsumesCollector && iC, std::vector<std::vector<int>>& selectedTkhidxSetIn): selectedTkhidxSet(selectedTkhidxSetIn) , bField_esToken_(iC.esConsumes<MagneticField, IdealMagneticFieldRecord>()) 
{

  // Get the track reco algorithm from the ParameterSet
  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));

  token_packedCandidates = iC.consumes<std::vector<pat::PackedCandidate>>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));

  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));



  // Second, initialize post-fit cuts
  mKPCutMin = theParameters.getParameter<double>(string("mKPCutMin"));
  mKPCutMax = theParameters.getParameter<double>(string("mKPCutMax"));
  mPiKPCutMin = theParameters.getParameter<double>(string("mPiKPCutMin"));
  mPiKPCutMax = theParameters.getParameter<double>(string("mPiKPCutMax"));
  tkDCACut = theParameters.getParameter<double>(string("tkDCACut"));
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  tkPtCut = theParameters.getParameter<double>(string("tkPtCut"));
  tkPtErrCut = theParameters.getParameter<double>(string("tkPtErrCut"));
  tkEtaCut = theParameters.getParameter<double>(string("tkEtaCut"));
//  tkPtSumCut = theParameters.getParameter<double>(string("tkPtSumCut"));
//  tkEtaDiffCut = theParameters.getParameter<double>(string("tkEtaDiffCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut2D = theParameters.getParameter<double>(string("collinearityCut2D"));
  collinCut3D = theParameters.getParameter<double>(string("collinearityCut3D"));
  lamCMassCut = theParameters.getParameter<double>(string("lamCMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  VtxChiProbCut = theParameters.getParameter<double>(string("VtxChiProbCut"));
  dPt3Cut = theParameters.getParameter<double>(string("dPt3Cut"));
  alphaCut = theParameters.getParameter<double>(string("alphaCut"));
  alpha2DCut = theParameters.getParameter<double>(string("alpha2DCut"));
  isWrongSign = theParameters.getParameter<bool>(string("isWrongSign"));

  //Nihar
  dPtCut_= theParameters.getParameter<double>(string("dPtCut"));
  dRapidityCut_= theParameters.getParameter<double>(string("dRapidityCut"));
  
  useAnyMVA_ = false;
  forestLabel_ = "LamC3PInpPb";
  std::string type = "BDT";
  useForestFromDB_ = true;
  dbFileName_ = "";

  forest_ = nullptr;

  if(theParameters.exists("useAnyMVA")) useAnyMVA_ = theParameters.getParameter<bool>("useAnyMVA");

  if(useAnyMVA_){
    if(theParameters.exists("mvaType"))type = theParameters.getParameter<std::string>("mvaType");
    if(theParameters.exists("GBRForestLabel"))forestLabel_ = theParameters.getParameter<std::string>("GBRForestLabel");
    if(theParameters.exists("GBRForestFileName")){
      dbFileName_ = theParameters.getParameter<std::string>("GBRForestFileName");
      useForestFromDB_ = false;
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeProducer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;
    mvaToken_ = iC.esConsumes<GBRForest, GBRWrapperRcd>(edm::ESInputTag("", forestLabel_));
  }

  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }
}

LamC3PFitter::~LamC3PFitter() {
  delete forest_;
}

// Method containing the algorithm for vertex reconstruction



void LamC3PFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  using namespace edm;
  using namespace reco;

  //std::cout << "[fitAll] Starting fitAll()" << std::endl;
  // Prepare custom track vectors
  std::vector<Track> lst;
  std::vector<TrackXYZP2> lstXYZP2;

  // Handles
  Handle<pat::PackedCandidateCollection> packedHandle;
  Handle<VertexCollection> vertexHandle;
  Handle<BeamSpot> beamSpotHandle;

  // Get data
  iEvent.getByToken(token_packedCandidates, packedHandle);
  iEvent.getByToken(token_vertices, vertexHandle);
  iEvent.getByToken(token_beamSpot, beamSpotHandle);




  if (!packedHandle.isValid() || packedHandle->empty()) return;

  // Magnetic field (needed for other things in your module)
  auto bFieldHandle = iSetup.getHandle(bField_esToken_);
  magField = bFieldHandle.product();

  // Get primary vertex or fallback to beamspot
  reco::Vertex thePrimaryV;
  //bool isVtxPV = false;
  
  if (!vertexHandle->empty() && !vertexHandle->begin()->isFake() && vertexHandle->begin()->tracksSize() >= 2) {
    thePrimaryV = vertexHandle->front();
    //isVtxPV = true;
  } else {
    thePrimaryV = reco::Vertex(beamSpotHandle->position(), reco::Vertex::Error(), 1, 1, 1);
    //isVtxPV = false;
  }



  reco::Vertex theBeamSpotV = reco::Vertex(beamSpotHandle->position(), reco::Vertex::Error(), 1, 1, 1);

  /*reco::Vertex theBeamSpotV(
    reco::Vertex::Point(beamSpotHandle->x(), beamSpotHandle->y(), beamSpotHandle->z()),
    reco::Vertex::Error(), 0.0, 0, 0
    );*/

  // Build input track list
  std::vector<pat::PackedCandidate> input_tracks;
  for (size_t i = 0; i < packedHandle->size(); ++i) {


    const auto& cand = (*packedHandle)[i];
    //std::cout << ">>> trackHandle size: " << packedHandle->size() << std::endl;


    
    if (!cand.hasTrackDetails()) continue;  // skip candidates without track info
    if (cand.pt() < tkPtCut || fabs(cand.eta()) > tkEtaCut) continue;

    // Optional quality filters (you can expand this)
    if (cand.pseudoTrack().normalizedChi2() >= tkChi2Cut ||
        cand.pseudoTrack().numberOfValidHits() < tkNhitsCut ||
	(cand.bestTrack() && (cand.bestTrack()->ptError() / cand.pt()) >= tkPtErrCut)) continue;


    //(cand.ptError() / cand.pt()) >= tkPtErrCut) continue;

    double dzvtx = cand.dz(thePrimaryV.position());
    double dxyvtx = cand.dxy(thePrimaryV.position());
    double dzerror = sqrt(cand.dzError() * cand.dzError() + thePrimaryV.zError() * thePrimaryV.zError());
    double dxyerror = sqrt(cand.dxyError() * cand.dxyError() +
                           thePrimaryV.xError() * thePrimaryV.yError());

    if (fabs(dzvtx / dzerror) <= dauLongImpactSigCut) continue;
    if (fabs(dxyvtx / dxyerror) <= dauTransImpactSigCut) continue;

    // Fill track structures
    input_tracks.push_back(cand);
    //input_tracks.push_back(cand);
    lst.emplace_back(cand.pt(), cand.eta(), cand.phi(), cand.charge(), i);
    lstXYZP2.emplace_back(lst.back());  // from Track to TrackXYZP2

    //std::cout << "Number of input_tracks: " << input_tracks.size() << std::endl;


  }

  

  // Mass window
  float mass_window[2] = {2.1, 2.4};  // Example for Lambda_c — adjust as needed


  // Call the updated fit function
  fitLamCCandidates(input_tracks, thePrimaryV, theBeamSpotV, lst, lstXYZP2, mass_window, iSetup);
  //fitLamCCandidates(packedHandle, thePrimaryV, theBeamSpotV, lst, lstXYZP2, mass_window, iSetup);
}




//long long int cnt0=0;
/*
  void LamC3PFitter::TkCombinationPermutation_Lc_v3(
						  const std::vector<pat::PackedCandidate> input_tracks, 
						  std::vector<Track> lst,
						  std::vector<TrackXYZP2> lstXYZP2,
						  float *mass_window,
						  //double tktkRes_mass,
						  //double tktkRes_mass_window,
						  std::vector<TrackSelected>& selectedTkhidxSet
						  //int Dchannel_number
						  ){
  
	int tk1_hindex = -1;
	int tk2_hindex = -1;
	int tk3_hindex = -1;

	mass_window[0] = 2.1; //Check!!!
	mass_window[1] = 2.5; //Check!!!

	
	//auto t0 = std::chrono::high_resolution_clock::now();	
	int number_NeededTrack = (int) lst.size();

	std::cout<<"number_NeededTrack="<<number_NeededTrack<<std::endl;
	
	for(int tk1idx = 0; tk1idx < number_NeededTrack; tk1idx++){
	  const TrackXYZP2& tr1 = lstXYZP2[tk1idx];
	  tk1_hindex = tr1.index;
	  int perm1 = tr1.q;  // q1+1
	  double p1sq = tr1.p2;
	  for(int tk2idx = tk1idx+1; tk2idx < number_NeededTrack; tk2idx++){
	    const TrackXYZP2& tr2 = lstXYZP2[tk2idx];
	    tk2_hindex = tr2.index;
	    int perm2 = (perm1 << 1) + tr2.q; // 2(q1+1) + (q2+1)
	    double p2sq = tr2.p2;
	    P3 p12(tr2);
	    p12 += tr1;
	    for(int tk3idx = tk2idx+1; tk3idx < number_NeededTrack; tk3idx++){
	      const TrackXYZP2& tr3 = lstXYZP2[tk3idx];
	      tk3_hindex = tr3.index;
	      int perm3 = (perm2 << 1) + tr3.q; // 4(q1+1) + 2(q2+1) + (q3+1)
	      if (perm3 == 0 || perm3 == 14) continue; //this remove the useless permutation
	      double p3sq = tr3.p2;
	      P3 pD(tr3);
	      pD += p12;
	      
	      // pT cut (independent of permutations)
	      double ptD2 = pD.px * pD.px + pD.py * pD.py;
	      //if(ptD2 < (dPtCut_[Dchannel_number-1])*(dPtCut_[Dchannel_number-1]))continue;
	      if(ptD2 < dPtCut_*dPtCut_)continue;//Nihar
	      double pzD2 = pD.pz * pD.pz;
	      for (int p = 0; p < 2; p ++) {
		double p0 = Functs.totE(perm3 + p, p1sq, p2sq, p3sq);
		double mD2 = p0 * p0 - ptD2 - pzD2;
		if(mD2 <(mass_window[0])*(mass_window[0]) || mD2 >(mass_window[1])*(mass_window[1])) continue;
		double mtD2 = mD2 + ptD2;
		//if (pzD2 > shymax2 * mtD2) continue; //this rapdity needs to be changed.
		//if (pzD2 > sinh(dRapidityCut_[Dchannel_number-1])*sinh(dRapidityCut_[Dchannel_number-1]) * mtD2) continue;

		if (pzD2 > sinh(dRapidityCut_)*sinh(dRapidityCut_) * mtD2) continue;//Nihar

		

		//cnt0++;
		
 
		selectedTkhidxSet.push_back(TrackSelected(tk1_hindex,tk2_hindex,tk3_hindex,perm3+p));//here also need to store the permutation number//here already have all the permutation		
		
		continue;
	      }//p
	    }//tk3id
	  }//tk2id
	}//tk1id
	return;
}
*/

void LamC3PFitter::TkCombinationPermutation_Lc_v3(
						  const std::vector<pat::PackedCandidate> input_tracks, 
						  std::vector<Track> lst,
						  std::vector<TrackXYZP2> lstXYZP2,
						  float *mass_window,
						  std::vector<TrackSelected>& selectedTkhidxSet
						  
						  ){
  
  int tk1_hindex = -1;
  int tk2_hindex = -1;
  int tk3_hindex = -1;

	mass_window[0] = 2.1; //Check!!!
	mass_window[1] = 2.5; //Check!!!

	
	//auto t0 = std::chrono::high_resolution_clock::now();	
	int number_NeededTrack = (int) lst.size();

	//std::cout<<"number_NeededTrack="<<number_NeededTrack<<std::endl;
	
	for (int tk1idx = 0; tk1idx < number_NeededTrack; tk1idx++) {
	  const TrackXYZP2& tr1 = lstXYZP2[tk1idx];
	  tk1_hindex = tr1.index;
	  int perm1 = tr1.q;
	  
	  for (int tk2idx = tk1idx + 1; tk2idx < number_NeededTrack; tk2idx++) {
	    const TrackXYZP2& tr2 = lstXYZP2[tk2idx];
	    tk2_hindex = tr2.index;
		int perm2 = (perm1 << 1) + tr2.q;
	    P3 p12(tr2);
	    p12 += tr1;
	    
	    for (int tk3idx = tk2idx + 1; tk3idx < number_NeededTrack; tk3idx++) {
	      const TrackXYZP2& tr3 = lstXYZP2[tk3idx];
	      tk3_hindex = tr3.index;
		  int perm3 = (perm2 << 1) + tr3.q;
	      if (perm3 == 0 || perm3 == 14) continue;
	      
	      P3 pD(tr3);
	      pD += p12;
	      
	      double ptD2 = pD.px * pD.px + pD.py * pD.py;
	      if (ptD2 < dPtCut_ * dPtCut_) continue;
	      
	      double pzD2 = pD.pz * pD.pz;
	      for (int p = 0; p < 2; p++) {
                double p0 = Functs.totE(perm3 + p, tr1.p2, tr2.p2, tr3.p2);
                double mD2 = p0 * p0 - ptD2 - pzD2;
                if (mD2 < mass_window[0] * mass_window[0] ||
                    mD2 > mass_window[1] * mass_window[1]) continue;
		
                double mtD2 = mD2 + ptD2;
                if (pzD2 > sinh(dRapidityCut_) * sinh(dRapidityCut_) * mtD2) continue;
		
		        selectedTkhidxSet.push_back(TrackSelected(tk1_hindex,tk2_hindex,tk3_hindex, perm3 + p));
                //selectedTkhidxSet.emplace_back(tk1idx, tk2idx, tk3idx, perm3 + p);
			  continue;
	      }
	    }
	  }
	}
	std::cout<<"TkCombinationPermutation, selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
	return;
}


//Define all 

double Mass_in_permutation[16][3] = { 
			   {  0.,              0.,                   0.   }, // 0 (- - -)
			   {  0.,              0.,                   0.   },
			   {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 2 (- - +)  perm2+p=2
			   {  PION_MASS,       PROTON_MASS,          KAON_MASS  },//perm2+p=3
			   {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 4 (- + -)
			   {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
			   {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 6 (- + +)
			   {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
			   {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 8 (+ - -)
			   {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
			   {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 10 (+ - +)
			   {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
			   {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 12 (+ + -)
			   {  PION_MASS,       PROTON_MASS,          KAON_MASS  },
			   {  0.,              0.,                   0.   }, // 14 (+ + +)
			   {  0.,              0.,                   0.   }
		};


  void LamC3PFitter::fitLamCCandidates(
				       //DInfoBranches &DInfo, 
				       //std::vector<pat::PackedCandidate> input_tracks,
				       const std::vector<pat::PackedCandidate> input_tracks, 
				       reco::Vertex thePrimaryV,
				       reco::Vertex theBeamSpotV,
				       vector<Track> lst,
				       vector<TrackXYZP2> lstXYZP2,
				       //std::vector<int> &D_counter,
				       float * mass_window,
				       //double tktkRes_mass,
				       //double tktkRes_mass_window,
				       //bool doConstrainFit,
				       //bool SequentialFit,
				       //int Dchannel_number,
				       //int TkCombinationMethod
				       const edm::EventSetup &iSetup
				       ){

    vector<TrackSelected> selectedTkhidxSet;

    TkCombinationPermutation_Lc_v3( input_tracks,lst,lstXYZP2, mass_window, selectedTkhidxSet);

    float chi = 0.;
    float ndf = 0.;
    
    //particle factory: produce transient tracks
    KinematicParticleFactoryFromTransientTrack pFactory;
    VirtualKinematicParticleFactory vFactory;
    //fitter for D
    KinematicParticleVertexFitter   tktk_fitter;
    RefCountedKinematicTree         tktk_VFT;
    RefCountedKinematicParticle     tktk_VFP;
    RefCountedKinematicVertex       tktk_VFPvtx;
    //constrain fit fitter
    KinematicConstrainedVertexFitter kcv_tktk_fitter;
    //fitter for Res
    KinematicParticleVertexFitter   tktkRes_fitter;
    RefCountedKinematicTree         tktkRes_VFT;
    RefCountedKinematicParticle     tktkRes_VFP;
    RefCountedKinematicVertex       tktkRes_VFPvtx;

    //const edm::EventSetup& iSetup
    //edm::ESHandle<MagneticField> bFieldHandle;
    //iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
    //const MagneticField *field = bFieldHandle.product();

    const MagneticField& bField = iSetup.getData(bField_esToken_);
    const MagneticField* field = &bField;
  
    //const MagneticField *field = bField.product();
    //AnalyticalImpactPointExtrapolator extrapolator(field);
    TransverseImpactPointExtrapolator extrapolator(field);
    TrajectoryStateOnSurface tsos;
    
    TLorentzVector v4_tk;
    std::vector<TLorentzVector> tktk_4vecs;//fitted tks
    TLorentzVector tktk_4vec;//fitted D
    TLorentzVector unfitted_tktk_4vec;//unfitted D
    std::vector<TLorentzVector> tktkRes_4vecs;//fitted Res tks
    TLorentzVector tktkRes_4vec;//fitted Res
    TLorentzVector unfitted_tktkRes_4vec;//unfitted Res
    std::vector<RefCountedKinematicParticle> tktk_candidate;//input tracks to D fitter
    std::vector<RefCountedKinematicParticle> tktkRes_candidate;//input tracks to Res fitter
    std::vector<RefCountedKinematicParticle> tktkCands;//output tracks from D fitter
    std::vector<RefCountedKinematicParticle> tktkResCands;//output tracks from Res fitter
    TLorentzVector temp_vec;//for temporary usage

    //std::cout<<"UPTO THIS IS OKAY----1"<<std::endl;    
    
    for(int i = 0; i < int(selectedTkhidxSet.size()); i++){

      //std::cout<< "selected candidates="<<int(selectedTkhidxSet.size())<<std::endl;


      
      //clear before using
      v4_tk.Clear();
      tktk_4vecs.clear();
      tktk_4vec.Clear();
      unfitted_tktk_4vec.Clear();
      tktkRes_4vecs.clear();
      tktkRes_4vec.Clear();
      unfitted_tktkRes_4vec.Clear();
      tktk_candidate.clear();
      tktkRes_candidate.clear();
      tktkCands.clear();
      tktkResCands.clear();
      
      unfitted_tktk_4vec.SetPxPyPzE(0., 0., 0., 0.);
      unfitted_tktkRes_4vec.SetPxPyPzE(0., 0., 0., 0.);
      
      //push back the Res tracks as first tracks
      ParticleMass tk_mass;
      //keep track of the push_back track index
      std::vector<int> pushbackTrkIdx;
      std::vector<int> pushbackResTrkIdx;
      std::vector<float> pushbackTrkMassHypo;
      std::vector<float> pushbackResTrkMassHypo;
      float tk_sigma;
      
      pushbackTrkMassHypo.push_back(0);
      pushbackTrkMassHypo.push_back(0);
      pushbackTrkMassHypo.push_back(0);
      



      int idx1 = selectedTkhidxSet[i].index_tk1;
      int idx2 = selectedTkhidxSet[i].index_tk2;
      int idx3 = selectedTkhidxSet[i].index_tk3;
      int permIdx = selectedTkhidxSet[i].permutation_number;

      //std::cout<<"idx1="<<idx1<<"  "<<"idx2="<<idx2<<"  "<<"idx3="<<idx3<<std::endl;
      //std::cout<<"input_tracks size="<<(int)input_tracks.size()<<std::endl;
      //std::cout<<"permIdx="<<permIdx<<std::endl;

      /*std::cout << "Checking triplet candidate #" << i << std::endl;
      std::cout << "permIdx = " << permIdx 
		<< ", idx1 = " << idx1 
		<< ", idx2 = " << idx2 
		<< ", idx3 = " << idx3 
		<< ", input_tracks.size() = " << input_tracks.size() 
		<< std::endl;
      */
      
      //if (permIdx < 0 || permIdx >= 16) continue;  // Skip this candidate
      //if (idx1 < 0 || idx1 >= (int)input_tracks.size()) continue;
      //if (idx2 < 0 || idx2 >= (int)input_tracks.size()) continue;
      //if (idx3 < 0 || idx3 >= (int)input_tracks.size()) continue;
      
      
      
/*temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk1].pt(), input_tracks[selectedTkhidxSet[i].index_tk1].eta(), input_tracks[selectedTkhidxSet[i].index_tk1].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0]);
      unfitted_tktk_4vec += temp_vec;
      
      temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk2].pt(), input_tracks[selectedTkhidxSet[i].index_tk2].eta(), input_tracks[selectedTkhidxSet[i].index_tk2].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1]);
      unfitted_tktk_4vec += temp_vec;
      
      temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk3].pt(), input_tracks[selectedTkhidxSet[i].index_tk3].eta(), input_tracks[selectedTkhidxSet[i].index_tk3].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2]);
      unfitted_tktk_4vec += temp_vec;
  
*/    

      

      //std::cout<<"UPTO THIS IS OKAY----2"<<std::endl;          
      reco::TransientTrack tkTT1(input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack(), &(*field) );
      if (tkTT1.isValid())
	{
	  //tk_mass = fabs(TkMassCharge[0].first);
	  tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT1,tk_mass,chi,ndf,tk_sigma));
	  pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk1);
	}
      reco::TransientTrack tkTT2(input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack(), &(*field) );
      if (tkTT2.isValid())
	{
	  //tk_mass = fabs(TkMassCharge[1].first);
	  tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT2,tk_mass,chi,ndf,tk_sigma));
	  pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk2);
	}
      reco::TransientTrack tkTT3(input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack(), &(*field) );
      if (tkTT3.isValid())
	{
	  //tk_mass = fabs(TkMassCharge[2].first);
	  tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2];
	  tk_sigma = Functs.getParticleSigma(tk_mass);
	  tktk_candidate.push_back(pFactory.particle(tkTT3,tk_mass,chi,ndf,tk_sigma));
	  pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk3);
	}

      //std::cout<<"UPTO THIS IS OKAY----3"<<std::endl;          
      double px1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().px();
      double py1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().py();
      double pz1 = input_tracks[selectedTkhidxSet[i].index_tk1].pseudoTrack().pz();
      double m1 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0];
      double E1 = sqrt(m1*m1 + px1*px1 + py1*py1 + pz1*pz1);
      
      double px2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().px();
      double py2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().py();
      double pz2 = input_tracks[selectedTkhidxSet[i].index_tk2].pseudoTrack().pz();
      double m2 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1];
      //double E2 = sqrt(m1*m1 + px1*px1 + py1*py1 + pz1*pz1);
      double E2 = sqrt(m2*m2 + px2*px2 + py2*py2 + pz2*pz2);
      
      double px3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().px();
      double py3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().py();
      double pz3 = input_tracks[selectedTkhidxSet[i].index_tk3].pseudoTrack().pz();
      double m3 = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2];
      //double E3 = sqrt(m1*m1 + px1*px1 + py1*py1 + pz1*pz1);
      double E3 = sqrt(m3*m3 + px3*px3 + py3*py3 + pz3*pz3);
      
      double invmasssquare = (E1+E2+E3)*(E1+E2+E3) - (px1+px2+px3)*(px1+px2+px3) - (py1+py2+py3)*(py1+py2+py3) - (pz1+pz2+pz3)*(pz1+pz2+pz3);
      double invmass = sqrt(invmasssquare);
      //std::cout<<"The inv mass : "<<invmass<<std::endl;
      
      //double MaximumDoca = Functs.getMaxDoca(tktk_candidate);
      tktk_VFT = tktk_fitter.fit(tktk_candidate);

      //if (!tktk_VFT || !tktk_VFT->isValid()) continue;
      if (!tktk_VFT || !tktk_VFT->isValid()) {
	//std::cout << "[ERROR] Kinematic fit failed!" << std::endl;
	continue;
      }

      
      //if(!tktk_VFT->isValid()) continue;
      
      //tktk_VFT->movePointerToTheTop(); // KinematicTree.cc , make the Tree accessible, pointer to particle, and daughters
      //tktk_VFP   = tktk_VFT->currentParticle();
      //tktk_VFPvtx = tktk_VFT->currentDecayVertex();
      //if (!tktk_VFPvtx->vertexIsValid()) continue;

      tktk_VFT->movePointerToTheTop();
      RefCountedKinematicParticle lamcCand = tktk_VFT->currentParticle();
      RefCountedKinematicVertex vtx = tktk_VFT->currentDecayVertex();
      //if (!vtx->vertexIsValid()) continue;

      if (!vtx || !vtx->vertexIsValid()) {
	std::cout << "[DEBUG] vtx is null or invalid!" << std::endl;
	continue;
      }

      
      const pat::PackedCandidate& cand1 = input_tracks[idx1];
      const pat::PackedCandidate& cand2 = input_tracks[idx2];
      const pat::PackedCandidate& cand3 = input_tracks[idx3];

      //Need to comment later
      //double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
      //if(chi2_prob_tktk < VtxChiProbCut) continue; 
      


      /*const GlobalError & cov3d = vtx->error();
      //Error4D cov4d(cov3d);  // or construct manually if needed
      
      double time = 0.;
      double chi2 = vtx->chiSquared();
      double ndof = vtx->degreesOfFreedom();
      size_t size = 1;

      //reco::Vertex lamCVtx(vtx->position(), cov4d, time, chi2, ndof, size);
      //const reco::Vertex lamCVtx(vtx->vertexState().position(), vtx->vertexState().error(), vtx->chiSquared(), vtx->degreesOfFreedom());
      
      
      const reco::Vertex lamCVtx(vtx->vertexState().position(),
                          vtx->vertexState().error(),
                          0.0,
                          vtx->chiSquared(),
                          vtx->degreesOfFreedom(),
                          1);

      //const reco::Vertex lamCVtx(vtx->position(), vtx->error(), vtx->chiSquared(), vtx->degreesOfFreedom());
      const reco::Particle::Point lamCVtxPos(vtx->position().x(), vtx->position().y(), vtx->position().z());
      */

      // Convert GlobalPoint (float-based) to reco::Vertex::Point (double-based)
      reco::Vertex::Point position(vtx->vertexState().position().x(),
				   vtx->vertexState().position().y(),
				   vtx->vertexState().position().z());
      
      // Convert GlobalError to reco::Vertex::Error (double-based 3x3 matrix)
      reco::Vertex::Error error;
      AlgebraicSymMatrix33 covMatrix = vtx->vertexState().error().matrix();
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  error(i, j) = covMatrix(i, j);
      
      // Now construct the reco::Vertex using the 5-argument constructor
      const reco::Vertex lamCVtx(position, error, vtx->chiSquared(), vtx->degreesOfFreedom(), 1);  // size
      
      
      const reco::Particle::LorentzVector lamCP4(lamcCand->currentState().globalMomentum().x(),
                                           lamcCand->currentState().globalMomentum().y(),
                                           lamcCand->currentState().globalMomentum().z(),
                                           lamcCand->currentState().mass());

      
      // Extract charges and build daughter candidates using original PackedCandidates
      int lamCCharge = cand1.charge() + cand2.charge() + cand3.charge();


      if (abs(lamCCharge) != 1) continue;  

      
      reco::Particle::LorentzVector p4_cand1(cand1.px(), cand1.py(), cand1.pz(), cand1.energy());
      RecoChargedCandidate theCand1(cand1.charge(), p4_cand1, lamCVtx.position());
      //theCand1.setTrack(cand1.bestTrack());
      
      reco::Particle::LorentzVector p4_cand2(cand2.px(), cand2.py(), cand2.pz(), cand2.energy());
      RecoChargedCandidate theCand2(cand2.charge(), p4_cand2, lamCVtx.position());
      //theCand2.setTrack(cand2.bestTrack());
      
      reco::Particle::LorentzVector p4_cand3(cand3.px(), cand3.py(), cand3.pz(), cand3.energy());
      RecoChargedCandidate theCand3(cand3.charge(), p4_cand3, lamCVtx.position());
      //theCand3.setTrack(cand3.bestTrack());

      
      // Build the Lambda_c candidate
      
      math::XYZPoint decayVertex(vtx->position().x(), vtx->position().y(), vtx->position().z());
      reco::VertexCompositeCandidate* theLamC3P = new reco::VertexCompositeCandidate(lamCCharge, lamCP4, decayVertex, static_cast<reco::VertexCompositeCandidate::CovarianceMatrix>(vtx->error().matrix()), vtx->chiSquared(), vtx->degreesOfFreedom());


      /*theLamC3P->addDaughter(theCand1);
      theLamC3P->addDaughter(theCand2);
      theLamC3P->addDaughter(theCand3);*/

      theLamC3P->addDaughter(input_tracks[idx1]);  // same as cand1
      theLamC3P->addDaughter(input_tracks[idx2]);
      theLamC3P->addDaughter(input_tracks[idx3]);


      int pdg_id = 4122;
      theLamC3P->setPdgId(pdg_id);
      
      // Update four-momentum using daughters
      AddFourMomenta addp4;
      addp4.set(*theLamC3P);


      //std::cout << "[DEBUG] LamC3P candidate mass: " << theLamC3P->mass() << std::endl;

      
      // Apply mass cut and store
      if (theLamC3P->mass() < lamCMassLamC3P + lamCMassCut &&
	  theLamC3P->mass() > lamCMassLamC3P - lamCMassCut) {
	theLamC3Ps.push_back(*theLamC3P);
	}
      
      if (theLamC3P) delete theLamC3P;

      //std::cout << "[DEBUG] Number of LamC3P candidates stored: " << theLamC3Ps.size() << std::endl;

      //std::cout<<"UPTO THIS IS OKAY----last"<<std::endl;          
      
    }
    
  }// End function
  









/*void LamC3PFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;
  using namespace std; 

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs_pos;
  std::vector<TrackRef> theTrackRefs_neg;
  std::vector<TransientTrack> theTransTracks_pos;
  std::vector<TransientTrack> theTransTracks_neg;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;


  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
  iEvent.getByToken(token_dedx, dEdxHandle);

  if( !theTrackHandle->size() ) return;
  bFieldHandle = iSetup.getHandle(bField_esToken_);

  magField = bFieldHandle.product();



  
  // Setup TMVA
  //  mvaValValueMap = auto_ptr<edm::ValueMap<float> >(new edm::ValueMap<float>);
  //  edm::ValueMap<float>::Filler mvaFiller(*mvaValValueMap);

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2) 
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
  math::XYZPoint bestvtxError(xVtxError,yVtxError,zVtxError);

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {// Track loop begin
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	if (tmpRef->quality(qualities[ndx_])){
	  quality_ok = true;
	  break;          
	}
      }
    }
    if( !quality_ok ) continue;

    //Nihar
    lst.push_back(
		  Track(
			tmpRef->pt(),
			tmpRef->eta(),
			tmpRef->phi(),
			tmpRef->charge() + 1,  
			indx                  
			)
		  );

    //int n = (int) lst.size();
    for (int i = 0; i < (int) lst.size(); i ++) {
      TrackXYZP2 tr = lst[i];
      lstXYZP2.push_back(tr);
    }

    
    if( tmpRef->normalizedChi2() < tkChi2Cut &&
        tmpRef->numberOfValidHits() >= tkNhitsCut &&
        tmpRef->ptError() / tmpRef->pt() < tkPtErrCut &&
        tmpRef->pt() > tkPtCut && fabs(tmpRef->eta()) < tkEtaCut ) {
      //TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
      TransientTrack tmpTk( *tmpRef, magField );

      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        if(tmpRef->charge()>0.0)
        {
          theTrackRefs_pos.push_back( tmpRef );
          theTransTracks_pos.push_back( tmpTk );
        }
        if(tmpRef->charge()<0.0)
        {
          theTrackRefs_neg.push_back( tmpRef );
          theTransTracks_neg.push_back( tmpTk );
        }
      }
    }
  }// Track loop end

  if(!isWrongSign)
  {
    fitLamCCandidates(theTrackRefs_pos,theTrackRefs_neg,theTransTracks_pos,theTransTracks_neg,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,4122);
    fitLamCCandidates(theTrackRefs_neg,theTrackRefs_pos,theTransTracks_neg,theTransTracks_pos,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,-4122);
  }
  else 
  {
    fitLamCCandidates(theTrackRefs_pos,theTrackRefs_pos,theTransTracks_pos,theTransTracks_pos,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,4122);
    fitLamCCandidates(theTrackRefs_neg,theTrackRefs_neg,theTransTracks_neg,theTransTracks_neg,isVtxPV,vtxPrimary,theBeamSpotHandle,bestvtx,bestvtxError,-4122);    
  }
}
*/

  
  /*void LamC3PFitter::fitLamCCandidates(
                                  std::vector<reco::TrackRef> theTrackRefs_sgn1,
                                  std::vector<reco::TrackRef> theTrackRefs_sgn2,
                                  std::vector<reco::TransientTrack> theTransTracks_sgn1,
                                  std::vector<reco::TransientTrack> theTransTracks_sgn2,
                                  bool isVtxPV,
                                  reco::VertexCollection::const_iterator vtxPrimary,
				  edm::Handle<reco::BeamSpot> theBeamSpotHandle,
                                  math::XYZPoint bestvtx,
				  math::XYZPoint bestvtxError,
                                  int pdg_id
                                )
{


  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  int lamCCharge = pdg_id/abs(pdg_id);

  int tk1_hindex = -1;
  int tk2_hindex = -1;
  int tk3_hindex = -1;
  
  std::vector<int> isNeededTrackIdx;
  
  // Loop over tracks and vertex good charged track pairs
  for(unsigned int trdx1 = 0; trdx1 < theTrackRefs_sgn1.size(); trdx1++) {


    tk1_hindex = isNeededTrackIdx[trdx1];
    
    for(unsigned int trdx2 = trdx1 + 1; trdx2 < theTrackRefs_sgn1.size(); trdx2++) {
      //if( (theTrackRefs[trdx1]->pt() + theTrackRefs[trdx2]->pt()) < tkPtSumCut) continue;
      //if( abs(theTrackRefs[trdx1]->eta() - theTrackRefs[trdx2]->eta()) > tkEtaDiffCut) continue;


    tk2_hindex = isNeededTrackIdx[trdx2];
      
      //This vector holds the pair of oppositely-charged tracks to be vertexed
      std::vector<TransientTrack> transTracks;

      TrackRef trackRef1 = theTrackRefs_sgn1[trdx1];
      TrackRef trackRef2 = theTrackRefs_sgn1[trdx2];
      TransientTrack* transTkPtr1 = &theTransTracks_sgn1[trdx1];
      TransientTrack* transTkPtr2 = &theTransTracks_sgn1[trdx2];


      // Fill the vector of TransientTracks to send to KVF

      //Nihar commented!
      //transTracks.push_back(*transTkPtr1);
      //transTracks.push_back(*transTkPtr2);

      // Trajectory states to calculate DCA for the 2 tracks
      FreeTrajectoryState trkState1 = transTkPtr1->impactPointTSCP().theState();
      FreeTrajectoryState trkState2 = transTkPtr2->impactPointTSCP().theState();

      if( !transTkPtr1->impactPointTSCP().isValid() || !transTkPtr2->impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(trkState1, trkState2);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );
      GlobalPoint cxPt = cApp.crossingPoint();

      if (dca < 0. || dca > tkDCACut) continue;

      // Get trajectory states for the tracks at POCA for later cuts
      TrajectoryStateClosestToPoint trkTSCP1 =
        transTkPtr1->trajectoryStateClosestToPoint( cxPt );
      TrajectoryStateClosestToPoint trkTSCP2 =
        transTkPtr2->trajectoryStateClosestToPoint( cxPt );

      if( !trkTSCP1.isValid() || !trkTSCP2.isValid() ) continue;

      double totalE1 = sqrt( trkTSCP1.momentum().mag2() + protonMassLamC3PSquared ) +
                      sqrt( trkTSCP2.momentum().mag2() + piMassLamC3PSquared );
      double totalE1Sq = totalE1*totalE1;

      double totalE2 = sqrt( trkTSCP1.momentum().mag2() + piMassLamC3PSquared ) +
                      sqrt( trkTSCP2.momentum().mag2() + protonMassLamC3PSquared );
      double totalE2Sq = totalE2*totalE2;

      double totalPSq =
        ( trkTSCP1.momentum() + trkTSCP2.momentum() ).mag2();

//      double totalPt =
//        ( trkTSCP1.momentum() + trkTSCP2.momentum() ).perp();

      double mass1 = sqrt( totalE1Sq - totalPSq);
      double mass2 = sqrt( totalE2Sq - totalPSq);

      if( (mass1 > mKPCutMax || mass1 < mKPCutMin) && (mass2 > mKPCutMax || mass2 < mKPCutMin)) continue;
//      if( totalPt < dPtCut ) continue;

      for(unsigned int trdx3 = 0; trdx3 < theTrackRefs_sgn2.size(); trdx3++) {


	tk3_hindex = isNeededTrackIdx[trdx3];
	
        TrackRef trackRef3 = theTrackRefs_sgn2[trdx3];
        TransientTrack* transTkPtr3 = &theTransTracks_sgn2[trdx3];
  
        transTracks.push_back(*transTkPtr3);
        FreeTrajectoryState trkState3 = transTkPtr3->impactPointTSCP().theState();
        if( !transTkPtr3->impactPointTSCP().isValid() ) continue;

        // Measure distance between tracks at their closest approach
        ClosestApproachInRPhi cApp13;
        cApp13.calculate(trkState1, trkState3);
        if( !cApp13.status() ) continue;
        float dca13 = fabs( cApp13.distance() );
        GlobalPoint cxPt13 = cApp13.crossingPoint();
        if (dca13 < 0. || dca13 > tkDCACut) continue;

        // Get trajectory states for the tracks at POCA for later cuts
        TrajectoryStateClosestToPoint trkTSCP31 =
          transTkPtr3->trajectoryStateClosestToPoint( cxPt13 );

        if( !trkTSCP31.isValid() ) continue;

        double totalE31 = sqrt( trkTSCP1.momentum().mag2() + protonMassLamC3PSquared ) +
                          sqrt( trkTSCP2.momentum().mag2() + piMassLamC3PSquared ) + 
                          sqrt( trkTSCP31.momentum().mag2() + kaonMassLamC3PSquared );
        double totalE31Sq = totalE31*totalE31;

        double totalE32 = sqrt( trkTSCP1.momentum().mag2() + piMassLamC3PSquared ) +
                          sqrt( trkTSCP2.momentum().mag2() + protonMassLamC3PSquared ) + 
                          sqrt( trkTSCP31.momentum().mag2() + kaonMassLamC3PSquared );
        double totalE32Sq = totalE32*totalE32;

        double totalP3Sq =
          ( trkTSCP1.momentum() + trkTSCP2.momentum() + trkTSCP31.momentum()).mag2();

        double totalPt3 =
          ( trkTSCP1.momentum() + trkTSCP2.momentum() + trkTSCP31.momentum()).perp();

        double mass31 = sqrt( totalE31Sq - totalP3Sq);
        double mass32 = sqrt( totalE32Sq - totalP3Sq);

        if( (mass31 > mPiKPCutMax || mass31 < mPiKPCutMin) && (mass32 > mPiKPCutMax || mass32 < mPiKPCutMin)) continue;
        if( totalPt3 < dPt3Cut ) continue;

        // Create the vertex fitter object and vertex the tracks
        float cand1TotalE[2]={0.0};
        float cand2TotalE[2]={0.0};
        float lamCTotalE[2]={0.0};



	for(int i=0;i<2;i++)
        {
          //Creating a KinematicParticleFactory
          KinematicParticleFactoryFromTransientTrack pFactory;
        
          float chi = 0.0;
          float ndf = 0.0;

          vector<RefCountedKinematicParticle> lamCParticles;
          lamCParticles.push_back(pFactory.particle(*transTkPtr1,cand1Mass[i],chi,ndf,cand1Mass_sigma[i]));
          lamCParticles.push_back(pFactory.particle(*transTkPtr2,cand2Mass[i],chi,ndf,cand2Mass_sigma[i]));
          lamCParticles.push_back(pFactory.particle(*transTkPtr3,kaonMassLamC3P,chi,ndf,kaonMassLamC3P_sigma));

          KinematicParticleVertexFitter lamCFitter;
          RefCountedKinematicTree lamCVertex;
          lamCVertex = lamCFitter.fit(lamCParticles);

          if( !lamCVertex->isValid() ) continue;

          lamCVertex->movePointerToTheTop();
          RefCountedKinematicParticle lamCCand = lamCVertex->currentParticle();
          if(!lamCCand->currentState().isValid()) continue;

          RefCountedKinematicVertex lamCDecayVertex = lamCVertex->currentDecayVertex();
          if(!lamCDecayVertex->vertexIsValid()) continue;

  	  float lamCC2Prob = TMath::Prob(lamCDecayVertex->chiSquared(),lamCDecayVertex->degreesOfFreedom());
  	  if (lamCC2Prob < VtxChiProbCut) continue;

          lamCVertex->movePointerToTheFirstChild();
          RefCountedKinematicParticle cand1 = lamCVertex->currentParticle();
          lamCVertex->movePointerToTheNextChild();
          RefCountedKinematicParticle cand2 = lamCVertex->currentParticle();
          lamCVertex->movePointerToTheNextChild();
          RefCountedKinematicParticle cand3 = lamCVertex->currentParticle();

          if(!cand1->currentState().isValid() || !cand2->currentState().isValid() || !cand3->currentState().isValid()) continue;

          KinematicParameters cand1KP = cand1->currentState().kinematicParameters();
          KinematicParameters cand2KP = cand2->currentState().kinematicParameters();
          KinematicParameters cand3KP = cand3->currentState().kinematicParameters();

          GlobalVector lamCTotalP = GlobalVector (lamCCand->currentState().globalMomentum().x(),
                                                  lamCCand->currentState().globalMomentum().y(),
                                                  lamCCand->currentState().globalMomentum().z());

          GlobalVector cand1TotalP = GlobalVector(cand1KP.momentum().x(),cand1KP.momentum().y(),cand1KP.momentum().z());
          GlobalVector cand2TotalP = GlobalVector(cand2KP.momentum().x(),cand2KP.momentum().y(),cand2KP.momentum().z());
          GlobalVector cand3TotalP = GlobalVector(cand3KP.momentum().x(),cand3KP.momentum().y(),cand3KP.momentum().z());

          cand1TotalE[i] = sqrt( cand1TotalP.mag2() + cand1Mass[i]*cand1Mass[i] );
          cand2TotalE[i] = sqrt( cand2TotalP.mag2() + cand2Mass[i]*cand2Mass[i] );
          float cand3TotalE = sqrt( cand3TotalP.mag2() + kaonMassLamC3PSquared );

          lamCTotalE[i] = cand1TotalE[i] + cand2TotalE[i] + cand3TotalE;

          const Particle::LorentzVector lamCP4(lamCTotalP.x(), lamCTotalP.y(), lamCTotalP.z(), lamCTotalE[i]);

          Particle::Point lamCVtx((*lamCDecayVertex).position().x(), (*lamCDecayVertex).position().y(), (*lamCDecayVertex).position().z());
          std::vector<double> lamCVtxEVec;
          lamCVtxEVec.push_back( lamCDecayVertex->error().cxx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().cyx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().cyy() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czx() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czy() );
          lamCVtxEVec.push_back( lamCDecayVertex->error().czz() );
          SMatrixSym3D lamCVtxCovMatrix(lamCVtxEVec.begin(), lamCVtxEVec.end());
          const Vertex::CovarianceMatrix lamCVtxCov(lamCVtxCovMatrix);
          double lamCVtxChi2(lamCDecayVertex->chiSquared());
          double lamCVtxNdof(lamCDecayVertex->degreesOfFreedom());
          double lamCNormalizedChi2 = lamCVtxChi2/lamCVtxNdof;

          double rVtxMag = 99999.0; 
          double lVtxMag = 99999.0;
          double sigmaRvtxMag = 999.0;
          double sigmaLvtxMag = 999.0;
          double lamCAngle3D = -100.0;
          double lamCAngle2D = -100.0;

          GlobalVector lamCLineOfFlight = GlobalVector (lamCVtx.x() - bestvtx.x(),
                                                        lamCVtx.y() - bestvtx.y(),
                                                        lamCVtx.z() - bestvtx.z());

          SMatrixSym3D lamCTotalCov;
          if(isVtxPV) lamCTotalCov = lamCVtxCovMatrix + vtxPrimary->covariance();
          else lamCTotalCov = lamCVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

          SVector3 distanceVector3D(lamCLineOfFlight.x(), lamCLineOfFlight.y(), lamCLineOfFlight.z());
          SVector3 distanceVector2D(lamCLineOfFlight.x(), lamCLineOfFlight.y(), 0.0);

          lamCAngle3D = angle(lamCLineOfFlight.x(), lamCLineOfFlight.y(), lamCLineOfFlight.z(),
                          lamCTotalP.x(), lamCTotalP.y(), lamCTotalP.z());
          lamCAngle2D = angle(lamCLineOfFlight.x(), lamCLineOfFlight.y(), (float)0.0,
                          lamCTotalP.x(), lamCTotalP.y(), (float)0.0);

          lVtxMag = lamCLineOfFlight.mag();
          rVtxMag = lamCLineOfFlight.perp();
          sigmaLvtxMag = sqrt(ROOT::Math::Similarity(lamCTotalCov, distanceVector3D)) / lVtxMag;
          sigmaRvtxMag = sqrt(ROOT::Math::Similarity(lamCTotalCov, distanceVector2D)) / rVtxMag;

          if( lamCNormalizedChi2 > chi2Cut ||
              rVtxMag < rVtxCut ||
              rVtxMag / sigmaRvtxMag < rVtxSigCut ||
              lVtxMag < lVtxCut ||
              lVtxMag / sigmaLvtxMag < lVtxSigCut ||
              cos(lamCAngle3D) < collinCut3D || cos(lamCAngle2D) < collinCut2D || lamCAngle3D > alphaCut || lamCAngle2D > alpha2DCut
          ) continue;

          VertexCompositeCandidate* theLamC3P = 0;
          theLamC3P = new VertexCompositeCandidate(lamCCharge, lamCP4, lamCVtx, lamCVtxCov, lamCVtxChi2, lamCVtxNdof);

          RecoChargedCandidate
            theCand1(lamCCharge, Particle::LorentzVector(cand1TotalP.x(),
                                                  cand1TotalP.y(), cand1TotalP.z(),
                                                  cand1TotalE[i]), lamCVtx);
            theCand1.setTrack(trackRef1);

          RecoChargedCandidate
            theCand2(lamCCharge, Particle::LorentzVector(cand2TotalP.x(),
                                                   cand2TotalP.y(), cand2TotalP.z(),
                                                   cand2TotalE[i]), lamCVtx);
            theCand2.setTrack(trackRef2);

          RecoChargedCandidate
            theCand3(-lamCCharge, Particle::LorentzVector(cand3TotalP.x(),
                                                   cand3TotalP.y(), cand3TotalP.z(),
                                                   cand3TotalE), lamCVtx);
            theCand3.setTrack(trackRef3);

          AddFourMomenta addp4;
          theLamC3P->addDaughter(theCand1);
          theLamC3P->addDaughter(theCand2);
          theLamC3P->addDaughter(theCand3);
          theLamC3P->setPdgId(pdg_id);
          addp4.set( *theLamC3P );
          if( theLamC3P->mass() < lamCMassLamC3P + lamCMassCut &&
              theLamC3P->mass() > lamCMassLamC3P - lamCMassCut ) //&&
	     // theLamC3P->pt() > dPtCut ) {
          {
            theLamC3Ps.push_back( *theLamC3P );
          }


          if(theLamC3P) delete theLamC3P;
        } // swap mass
	*/


// Get methods

const reco::VertexCompositeCandidateCollection& LamC3PFitter::getLamC3P() const {
  return theLamC3Ps;
}

const std::vector<float>& LamC3PFitter::getMVAVals() const {
  return mvaVals_;
}


void LamC3PFitter::resetAll() {
    theLamC3Ps.clear();
    mvaVals_.clear();
}
