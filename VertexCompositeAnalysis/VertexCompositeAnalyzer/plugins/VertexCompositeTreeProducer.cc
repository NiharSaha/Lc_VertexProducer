// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "Utils.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


#define PI 3.1416
#define MAXCAN 20000

using namespace std;




class VertexCompositeTreeProducer : public edm::one::EDAnalyzer<> {
	public:
		explicit VertexCompositeTreeProducer(const edm::ParameterSet&);
		~VertexCompositeTreeProducer();

		using MVACollection = std::vector<float>;

	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void fillRECO(const edm::Event&, const edm::EventSetup&) ;
		virtual void endJob() ;
		virtual void initTree();
		void genDecayLength(const uint&, const reco::GenParticle&);


		// ----------member data ---------------------------

		edm::Service<TFileService> fs;

		//CommonFuncts        Functs;


		TTree* VertexCompositeNtuple;
		bool   saveTree_;

		//options
		bool doRecoNtuple_;
		bool dogenntuple_;   
  //bool dogenmatching_;
		bool dogenmatchingtof_;
		bool hasswap_;
		bool decayingen_;
  bool threeProngDecay_;
  int PID_;
		int PID_dau1_;
		int PID_dau2_;
		int PID_dau3_;

		//cut variables
		double multMax_;
		double multMin_;
		double deltaR_; //deltaR for Gen matching


		//tree branches
		//event info
		int centrality;
		int Ntrkoffline;
		int Npixel;
                float HFsumET;//Nihar
		float HFsumETPlus;
		float HFsumETMinus;
		float ZDCPlus;
		float ZDCMinus;
		float bestvx;
		float bestvy;
		float bestvz;
		float bestvxError;
		float bestvyError;
		float bestvzError;
		float BSx;
		float BSy;
		float BSz;
		float BSxerror;
		float BSyerror;
		float BSzerror;
		int candSize;
		float ephfpAngle[3];
		float ephfmAngle[3];
		float ephfpQ[3];
		float ephfmQ[3];
		float ephfpSumW;
		float ephfmSumW;

		//Composite candidate info
		float mva[MAXCAN];
		float pt[MAXCAN];
		float eta[MAXCAN];
		float phi[MAXCAN];
		float flavor[MAXCAN];
		float y[MAXCAN];
		float mass[MAXCAN];
		float VtxProb[MAXCAN];
		float dlos[MAXCAN];
		float dl[MAXCAN];
		float dlerror[MAXCAN];
		float DlxyBS[MAXCAN];
		float DlxyBSErr[MAXCAN];
		float vtxChi2[MAXCAN];
		float ndf[MAXCAN];
		float agl_abs[MAXCAN];
		float Ddca[MAXCAN];
		float agl2D_abs[MAXCAN];
		float dlos2D[MAXCAN];
		float dl2D[MAXCAN];
		float dl2Derror[MAXCAN];
		bool isSwap[MAXCAN];
		bool matchGEN[MAXCAN];
		int idmom_reco[MAXCAN];
		int idd1_reco[MAXCAN];
		int idd2_reco[MAXCAN];
  		int idd3_reco[MAXCAN];

  float gen_agl_abs[MAXCAN];
		float gen_agl2D_abs[MAXCAN];
		float gen_dl[MAXCAN];
		float gen_dl2D[MAXCAN];


		//dau info
		float dzos1[MAXCAN];
		float dzos2[MAXCAN];
  		float dzos3[MAXCAN];
		float dxyos1[MAXCAN];
		float dxyos2[MAXCAN];
  		float dxyos3[MAXCAN];
		float pt1[MAXCAN];
		float pt2[MAXCAN];
  		float pt3[MAXCAN];
		float ptErr1[MAXCAN];
		float ptErr2[MAXCAN];
  		float ptErr3[MAXCAN];
		float p1[MAXCAN];
		float p2[MAXCAN];
  		float p3[MAXCAN];
		float Dtrk1Dz1[MAXCAN];
		float Dtrk2Dz1[MAXCAN];
  		float Dtrk3Dz1[MAXCAN];
		float Dtrk1Dxy1[MAXCAN];
		float Dtrk2Dxy1[MAXCAN];
  		float Dtrk3Dxy1[MAXCAN];
		float Dtrk1DzError1[MAXCAN];
		float Dtrk2DzError1[MAXCAN];
  		float Dtrk3DzError1[MAXCAN];
		float Dtrk1DxyError1[MAXCAN];
		float Dtrk2DxyError1[MAXCAN];
  		float Dtrk3DxyError1[MAXCAN];
		float eta1[MAXCAN];
		float eta2[MAXCAN];
  		float eta3[MAXCAN];

                float phi1[MAXCAN];
		float phi2[MAXCAN];
  		float phi3[MAXCAN];

                int charge1[MAXCAN];
		int charge2[MAXCAN];
  		int charge3[MAXCAN];
		int pid1[MAXCAN];
		int pid2[MAXCAN];
  		int pid3[MAXCAN];
		float tof1[MAXCAN];
		float tof2[MAXCAN];
  		float tof3[MAXCAN];
		float H2dedx1[MAXCAN];
		float H2dedx2[MAXCAN];
  		float H2dedx3[MAXCAN];
		float T4dedx1[MAXCAN];
		float T4dedx2[MAXCAN];
  		float T4dedx3[MAXCAN];
		float trkChi1[MAXCAN];
		float trkChi2[MAXCAN];
  		float trkChi3[MAXCAN];
  
                //vector for gen match                                                                                              
                vector< vector<double> > *pVect;
                vector<double> *Dvector1;
                vector<double> *Dvector2;
                vector<double> *Dvector3;
                vector<int> *pVectIDmom;

                // gen info    
		float pt_gen[MAXCAN];
		float eta_gen[MAXCAN];
		int idmom[MAXCAN];
		float y_gen[MAXCAN];
		float phi_gen[MAXCAN];
		int iddau1[MAXCAN];
		int iddau2[MAXCAN];
		int iddau3[MAXCAN];

		
		bool isSkimMVA_;
		bool isCentrality_;
  //bool isData_cent_;
  //bool isMC_cent_;
		bool doGenNtuple_;
  bool useAnyMVA_;
  bool doGenMatching_;
		bool decayInGen_;

		edm::Handle<int> cbin_;

		//tokens
		edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tok_generalTrk_;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> patCompositeCandidateCollection_Token_;
  //edm::EDGetTokenT<MVACollection> MVAValues_Token_;

  //edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
  //edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
		edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;

		edm::EDGetTokenT<int> tok_centBinLabel_;
		edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

		edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventplaneSrc_;

		//abby
		edm::EDGetTokenT< reco::BeamSpot > bsLabel_;

};


VertexCompositeTreeProducer::VertexCompositeTreeProducer(const edm::ParameterSet& iConfig)

{
	//options
	doRecoNtuple_ = iConfig.getUntrackedParameter<bool>("doRecoNtuple");
	doGenNtuple_ = iConfig.getUntrackedParameter<bool>("doGenNtuple");
	doGenMatching_ = iConfig.getUntrackedParameter<bool>("doGenMatching");
	decayInGen_ = iConfig.getUntrackedParameter<bool>("decayInGen");
	PID_ = iConfig.getUntrackedParameter<int>("PID");
	PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
	PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
	PID_dau3_ = iConfig.getUntrackedParameter<int>("PID_dau3");

	saveTree_ = iConfig.getUntrackedParameter<bool>("saveTree");
	threeProngDecay_ = iConfig.getUntrackedParameter<bool>("threeProngDecay");
	//useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
	isSkimMVA_ = iConfig.getUntrackedParameter<bool>("isSkimMVA"); 
	isCentrality_ = iConfig.getParameter<bool>("isCentrality");
	//isData_cent_=iConfig.getUntrackedParameter<bool>("isData_cent");
	//isMC_cent_=iConfig.getUntrackedParameter<bool>("isMC_cent");

	//useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
	
	  //cut variables
	multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
	multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
	deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.02);


	//input tokens
	tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"));
	tok_generalTrk_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("TrackCollection"));
	patCompositeCandidateCollection_Token_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LamC3P"));
	//MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
	//Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
	//Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
	tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getParameter<edm::InputTag>("GenParticleCollection")));
	//abby
	bsLabel_        = consumes< reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("BSLabel"));

	
	
	if(isCentrality_)
	  {
	    tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
	    tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
	  }
	
	
	//if(useAnyMVA_ && iConfig.exists("MVACollection"))
	//	MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));

	
	

}


VertexCompositeTreeProducer::~VertexCompositeTreeProducer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup&
		iSetup)
{

  std::cout << "[DEBUG] doGenMatching_ = " << doGenMatching_<< std::endl;
  
  
  using std::vector;
  using namespace edm;
  using namespace reco;
  
  if(doRecoNtuple_) fillRECO(iEvent,iSetup);
  
  if(saveTree_) VertexCompositeNtuple->Fill();
}

void VertexCompositeTreeProducer::fillRECO(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //std::cout << "[DEBUG] useAnyMVA_ = " << useAnyMVA_ << "  "<<"doGenMatching_ = "<< doGenMatching_<<std::endl;


#ifdef DEBUG
	using std::cout;
	using std::endl;
#endif
	//get collections
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(tok_offlinePV_,vertices);

	/*if (!vertices.isValid() || vertices->empty()) {
	  edm::LogWarning("LamC3PAna") << "No primary vertices in this event"<< iEvent.id();
	  return; // safely skip this event
	  }*/

	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(bsLabel_, beamSpotHandle);

	/*if (!beamSpotHandle.isValid()) {
	  edm::LogWarning("LamC3PAna") << "No beam spot in this event"<< iEvent.id();
	  return;
	  }*/

	/*edm::Handle<pat::PackedCandidateCollection> tracks;
	if (!tracks.isValid() || tracks->empty()) {
	  edm::LogWarning("LamC3PAna") << "No general tracks in this event!";
	  return;
	}
	iEvent.getByToken(tok_generalTrk_, tracks);
	*/
	
	edm::Handle<pat::CompositeCandidateCollection> lamC3Pcandidates;
	iEvent.getByToken(patCompositeCandidateCollection_Token_,lamC3Pcandidates);

	/*if (!lamC3Pcandidates.isValid() || lamC3Pcandidates->empty()) {
	  edm::LogWarning("LamC3PAna") << "No LamC3P candidates in this event"<<iEvent.id();
	  return;
	  }*/
	
	const pat::CompositeCandidateCollection * lamC3Pcandidates_ = lamC3Pcandidates.product();



	//auto mvavalues = iEvent.getHandle(MVAValues_Token_);
	//edm::Handle<MVACollection> mvavalues;

	/*if (useAnyMVA_) {
	  std::cout<<"TO CHECK!!!!!!!"<<std::endl;
	  
	  iEvent.getByToken(MVAValues_Token_,mvavalues);
	  assert((*mvavalues).size() == lamC3Pcandidates_->size());
	  }*/
	

	edm::Handle<reco::GenParticleCollection> genpars;	
	if (doGenMatching_) {
	  iEvent.getByToken(tok_genParticle_, genpars);
	}

	/*edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
	iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
	if (!dEdxHandle1.isValid()) {
	  edm::LogWarning("LamC3PAna") << "dEdxHandle1 not valid!";
	  // skip or continue safely
	}
 
	edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
	iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
	if (!dEdxHandle2.isValid()) {
	  edm::LogWarning("LamC3PAna") << "dEdxHandle2 not valid!";
	  // skip or continue safely
	  }*/


#ifdef DEBUG
	cout << "Loaded tokens" << endl;
#endif

	centrality=-1;
	if(isCentrality_)
	{
		edm::Handle<reco::Centrality> cent;
		iEvent.getByToken(tok_centSrc_, cent);
		HFsumET = (cent.isValid() ? cent->EtHFtowerSum() : -1.); //Nihar
		HFsumETPlus = (cent.isValid() ? cent->EtHFtowerSumPlus() : -1.);
		HFsumETMinus = (cent.isValid() ? cent->EtHFtowerSumMinus() : -1.);
		Npixel = (cent.isValid() ? cent->multiplicityPixel() : -1);
		ZDCPlus = (cent.isValid() ? cent->zdcSumPlus() : -1.);
		ZDCMinus = (cent.isValid() ? cent->zdcSumMinus() : -1.);
		Ntrkoffline = (cent.isValid() ? cent->Ntracks() : -1);

		edm::Handle<int> cbin;
	  
		iEvent.getByToken(tok_centBinLabel_, cbin);
		centrality = (cbin.isValid() ? *cbin : -1);
		//if(isData_cent_){centrality = (cbin.isValid() ? *cbin : -1);}
		//if(isMC_cent_){centrality = getHiBinFromhiHF(HFsumET);}

		//cout<<"isMC_cent_="<<isMC_cent_<<"isData_cent_="<<isData_cent_<<endl;
	}

	//best vertex
	BSx=-999.9; BSy=-999.9; BSz=-999.9;
	BSxerror=-999.9; BSyerror=-999.9; BSzerror=-999.9;
	//reco::BeamSpot beamSpot;

	float BSdxdz = -999.9;
	float BSdydz = -999.9;
	bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
	bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;

	//const reco::Vertex & vtx;
	//const reco::BeamSpot & beamSpot;
	
	//if (!vertices.isValid() || vertices->empty()) continue; 
	//if (!beamSpotHandle.isValid()) continue;
	
	  const reco::Vertex & vtx = (*vertices)[0];
	//reco::Vertex & vtx = vertices->front();
	  bestvx = vtx.x();
	  bestvy = vtx.y();
	  bestvz = vtx.z();

	  bestvxError = vtx.xError();
	  bestvyError = vtx.yError();
	  bestvzError = vtx.zError();
	
	  float xlxyBS = -999.9;
	  float ylxyBS = -999.9;
	  float vtxYXErr = -999.9;
	  float vtxXErr = -999.9;
	  float vtxYErr = -999.9;


	  reco::BeamSpot beamSpot = *beamSpotHandle;
	  reco::Vertex theBeamSpotV(beamSpot.position(), beamSpot.covariance3D());
	  
	    BSx      = beamSpot.x0();
	    BSy      = beamSpot.y0();
	    BSz      = beamSpot.z0();
	    BSxerror = beamSpot.x0Error();
	    BSyerror = beamSpot.y0Error();
	    BSzerror = beamSpot.z0Error();
	    BSdxdz   = beamSpot.dxdz();
	    BSdydz   = beamSpot.dydz();


	//RECO Candidate info
	candSize = lamC3Pcandidates_->size();
	//std::cout << " Lc candidates in tree producer" <<candSize<< std::endl;

	if (candSize >= MAXCAN){
	  std::cout << "Warning: candSize (" << candSize << ") exceeds MAXCAND (" << MAXCAN << "). Data will be truncated." << std::endl;
	  
	  candSize = MAXCAN;
	}

	
	for(int it=0; it<candSize; ++it){

	  //if (candSize >= MAXCAN) break;
	  
	  
	  //std::cout<<"Getting candidate"<<std::endl;
	  const pat::CompositeCandidate & trk = (*lamC3Pcandidates_)[it];
	  //std::cout<<"Got candidate"<<std::endl;

	  double secvz=-999.9, secvx=-999.9, secvy=-999.9;
	  secvz = trk.userFloat("vtxZ"); secvx = trk.userFloat("vtxX"); secvy = trk.userFloat("vtxY");
	  bestvz = trk.userFloat("bestvtxZ");
	  bestvx = trk.userFloat("bestvtxX");
	  bestvy = trk.userFloat("bestvtxY");
	  bestvzError = trk.userFloat("zVtxError");
	  bestvxError = trk.userFloat("xVtxError");
	  bestvyError = trk.userFloat("yVtxError");
	  
	  //std::cout<<"CHECK-1"<<" for event "<<iEvent.id().event()<<std::endl;
	  reco::Vertex::CovarianceMatrix sec_covariance;
	  for (int i = 0; i < 3; i++)
	    {
	      for (int j = 0; j < 3; j++)
		{
		  sec_covariance(i, j) = trk.userFloat("vertexCovariance_" + std::to_string(i) + "_" + std::to_string(j));
		}
	    }
	  
	  vtxYXErr = sec_covariance(1, 0);
	  vtxXErr = sec_covariance(0, 0);
	  vtxYErr = sec_covariance(1, 1);
	  //std::cout<<"CHECK-2"<<" for event "<<iEvent.id().event()<<std::endl;
	  
	  
	  eta[it] = trk.eta();
	  y[it] = trk.rapidity();
	  pt[it] = trk.pt();
	  
	  phi[it] = trk.phi();
	  flavor[it] = trk.pdgId()/abs(trk.pdgId());
	  
	  //mva[it] = 0.0;
	  //if(useAnyMVA_) mva[it] = (*mvavalues)[it];

	  //std::cout<<"CHECK-3"<<" for event "<<iEvent.id().event()<<std::endl;
	  double px = trk.px();
	  double py = trk.py();
	  double pz = trk.pz();
	  mass[it] = trk.mass();
	  
	  if (mass[it] == 0) {
	    cout << "Error break" << endl;
	  }
	  
	  
	  //std::cout<<"CHECK-4"<<" for event "<<iEvent.id().event()<<std::endl;
	  //const reco::Candidate* cand1 = trk.daughter(0);
	  const pat::PackedCandidate* reco_d1 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(0));
	  
	  //const reco::Candidate* cand2 = trk.daughter(1);
	  const pat::PackedCandidate* reco_d2 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(1));


	const pat::PackedCandidate* reco_d3 = nullptr;
	if (threeProngDecay_) { 
	    reco_d3 = dynamic_cast<const pat::PackedCandidate*>(trk.daughter(2));
	  }

	//std::cout<<"CHECK-5"<<" for event "<<iEvent.id().event()<<std::endl;
	  
	  //Gen match
	  matchGEN[it] = false;
	  isSwap[it] = false;
	  idmom_reco[it] = -77;
	  idd1_reco[it] = -77;
	  idd2_reco[it] = -77;
	  idd3_reco[it] = -77;
	  
	  pt_gen[it] = -999.9;
	  eta_gen[it] = -999.9;
	  idmom[it] = -999;
	  y_gen[it] = -999.9;
	  phi_gen[it] = -999.9;
	  iddau1[it] = -999;
	  iddau2[it] = -999;
	  iddau3[it] = -999;

	  //std::cout<<"CHECK-1"<<std::endl;
	  
	  if(doGenMatching_){
	    
            if(!genpars.isValid()){
	      cout<<"Gen matching cannot be done without Gen collection!!"<<endl;
	      return;
	    }
	    
	    for(unsigned genPair=0; genPair<genpars->size(); ++genPair){
	      
	      const reco::GenParticle & genLamC = (*genpars)[genPair];
	      
	      //std::cout<<"CHECK-2"<<std::endl;
	      //std::cout<<"Lc pdgId="<<genLamC.pdgId()<<std::endl;
	      
	      /*if (fabs(genLamC.pdgId()) == 4122) {
		std::cout << "Found Λc candidate! PDGID = " << genLamC.pdgId()
			  << "   pt=" << genLamC.pt()
			  << "   eta=" << genLamC.eta()
			  << std::endl;
			  }*/

	      if(fabs(genLamC.pdgId()) != PID_) continue;
		
		const reco::Candidate * gen_d1 = genLamC.daughter(0);
		const reco::Candidate * gen_d2 = genLamC.daughter(1);
		const reco::Candidate * gen_d3 = genLamC.daughter(2);
		
		// check PDG IDs (3! = 6 possible permutations)


		//std::cout<<"PID_dau1="<<gen_d1->pdgId()<<std::endl;
		
		bool passPID = (
				(fabs(gen_d1->pdgId())==PID_dau1_ && fabs(gen_d2->pdgId())==PID_dau2_ && fabs(gen_d3->pdgId())==PID_dau3_) ||
				(fabs(gen_d1->pdgId())==PID_dau1_ && fabs(gen_d2->pdgId())==PID_dau3_ && fabs(gen_d3->pdgId())==PID_dau2_) ||
				(fabs(gen_d1->pdgId())==PID_dau2_ && fabs(gen_d2->pdgId())==PID_dau1_ && fabs(gen_d3->pdgId())==PID_dau3_) ||
				(fabs(gen_d1->pdgId())==PID_dau2_ && fabs(gen_d2->pdgId())==PID_dau3_ && fabs(gen_d3->pdgId())==PID_dau1_) ||
				(fabs(gen_d1->pdgId())==PID_dau3_ && fabs(gen_d2->pdgId())==PID_dau1_ && fabs(gen_d3->pdgId())==PID_dau2_) ||
				(fabs(gen_d1->pdgId())==PID_dau3_ && fabs(gen_d2->pdgId())==PID_dau2_ && fabs(gen_d3->pdgId())==PID_dau1_)
				);
		if(!passPID) continue;
		//std::cout<<"CHECK-3"<<std::endl;
		// To check all 6 charge permutations
		struct CandidateCombo { 
		  const reco::Candidate* g1; 
		  const reco::Candidate* g2; 
		  const reco::Candidate* g3; 
		};


		std::vector<CandidateCombo> gen_combos = {
		  {gen_d1, gen_d2, gen_d3},
		  {gen_d1, gen_d3, gen_d2},
		  {gen_d2, gen_d1, gen_d3},
		  {gen_d2, gen_d3, gen_d1},
		  {gen_d3, gen_d1, gen_d2},
		  {gen_d3, gen_d2, gen_d1}
		};
		
		//std::cout<<"CHECK-4"<<std::endl;
		//bool matched = false;
		for(auto &c : gen_combos) {
		  
		  if(reco_d1->charge()==c.g1->charge() &&
		     reco_d2->charge()==c.g2->charge() &&
		     reco_d3->charge()==c.g3->charge()) {
		    
		    // deltaR and deltaPt checks
		    double deltaR = sqrt(pow(reco_d1->eta()-c.g1->eta(),2)+pow(reco_d1->phi()-c.g1->phi(),2));
		    if(deltaR > deltaR_) continue;
		    if(fabs((reco_d1->pt()-c.g1->pt())/reco_d1->pt()) > 0.2) continue;

		    deltaR = sqrt(pow(reco_d2->eta()-c.g2->eta(),2)+pow(reco_d2->phi()-c.g2->phi(),2));
		    if(deltaR > deltaR_) continue;
		    if(fabs((reco_d2->pt()-c.g2->pt())/reco_d2->pt()) > 0.2) continue;
		    
		    deltaR = sqrt(pow(reco_d3->eta()-c.g3->eta(),2)+pow(reco_d3->phi()-c.g3->phi(),2));
		    if(deltaR > deltaR_) continue;
		    if(fabs((reco_d3->pt()-c.g3->pt())/reco_d3->pt()) > 0.2) continue;


		    //std::cout<<"CHECK-5"<<std::endl;
                // matched!
                matchGEN[it] = true;
                if(reco_d1->pdgId() != c.g1->pdgId()) isSwap[it] = true;
                genDecayLength(it, genLamC);

                pt_gen[it] = genLamC.pt();
                eta_gen[it] = genLamC.eta();
                y_gen[it] = genLamC.rapidity();
                phi_gen[it] = genLamC.phi();
                idmom[it] = genLamC.pdgId();

                if(!decayInGen_) continue;

		//std::cout<<"CHECK-6"<<std::endl;
                iddau1[it] = c.g1->pdgId();
                iddau2[it] = c.g2->pdgId();
                iddau3[it] = c.g3->pdgId();
		
                //matched = true;
                break;
		//std::cout<<"CHECK-7"<<std::endl;
		  }
		}
		//if(matched) break; // move to next reco candidate
		idmom_reco[it] = trk.pdgId();                                                                
		idd1_reco[it] = reco_d1->pdgId();
		idd2_reco[it] = reco_d2->pdgId();
		idd3_reco[it] = reco_d3->pdgId();


	    }//genPair loop

	    std::cout<<"End of doGenMatching!!!"<<std::endl;
	  }// End doGenMatching!!
	    
	   //std::cout<<"CHECK-9"<<std::endl;
	    





	  //std::cout<<"CHECK-6"<<" for event "<<iEvent.id().event()<<std::endl;
	  double pxd1 = reco_d1->px();
	  double pyd1 = reco_d1->py();
	  double pzd1 = reco_d1->pz();
	  double pxd2 = reco_d2->px();
	  double pyd2 = reco_d2->py();
	  double pzd2 = reco_d2->pz();
	  
	  
	  TVector3 dauvec1(pxd1,pyd1,pzd1);
	  TVector3 dauvec2(pxd2,pyd2,pzd2);

	  
	  //pt
	  pt1[it] = reco_d1->pt();
	  pt2[it] = reco_d2->pt();
	  
	  //momentum
	  p1[it] = reco_d1->p();
	  p2[it] = reco_d2->p();

	  //eta
	  eta1[it] = reco_d1->eta();
	  eta2[it] = reco_d2->eta();
	  
	  //phi
	  phi1[it] = reco_d1->phi();
	  phi2[it] = reco_d2->phi();
	  
	  //charge
	  charge1[it] = reco_d1->charge();
	  charge2[it] = reco_d2->charge();
	  
	  double pxd3 = -999.9;
	  double pyd3 = -999.9;
	  double pzd3 = -999.9;

	  
	  if(threeProngDecay_ )
	    {
	      //std::cout<<"CHECK-7"<<" for event "<<iEvent.id().event()<<std::endl;
	      pxd3 = reco_d3->px();
	      pyd3 = reco_d3->py();
	      pzd3 = reco_d3->pz();
	      pt3[it] = reco_d3->pt();
	      p3[it] = reco_d3->p();
	      eta3[it] = reco_d3->eta();
	      phi3[it] = reco_d3->phi();
	      charge3[it] = reco_d3->charge();
	      TVector3 dauvec3(pxd3,pyd3,pzd3);
	    }
	  
	  
	  
	  pid1[it] = -99999;
	  pid2[it] = -99999;
	  pid3[it] = -99999;

	  	  
	  //vtxChi2
	  vtxChi2[it] = trk.userFloat("vtxChi2");
	  ndf[it] = trk.userFloat("vtxNdof");
	  VtxProb[it] = TMath::Prob(vtxChi2[it],ndf[it]);

	  //std::cout<<"CHECK-8"<<" for event "<<iEvent.id().event()<<std::endl;
	  //PAngle
	  TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
	  TVector3 secvec(px,py,pz);
	  
	  TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
	  TVector3 secvec2D(px,py,0);
	  
	  
	  agl_abs[it] = secvec.Angle(ptosvec);
	  //Ddca[it] = ptosvec.Mag() * TMath::Sin(agl_abs[it]);
	  agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
	  	  
	  float r2lxyBS = (secvx-BSx-(secvz-BSz)*BSdxdz) * (secvx-BSx-(secvz-BSz)*BSdxdz) + (secvy-BSy-(secvz-BSz)*BSdydz) * (secvy-BSy-(secvz-BSz)*BSdydz);
	  xlxyBS = secvx-BSx - (secvz-BSz)*BSdxdz;
	  ylxyBS = secvy-BSy - (secvz-BSz)*BSdydz;
	  //abby std::cout << "r2lxyBS = " << r2lxyBS << std::endl;
	  DlxyBS[it] = static_cast<float>(TMath::Sqrt(r2lxyBS));
	  DlxyBSErr[it] = static_cast<float>(TMath::Sqrt ((1./r2lxyBS) * ((xlxyBS*xlxyBS)*vtxXErr + (2*xlxyBS*ylxyBS)*vtxYXErr + (ylxyBS*ylxyBS)*vtxYErr)));

	  //std::cout<<"CHECK-9"<<" for event "<<iEvent.id().event()<<std::endl;
	  //Decay length 3D
	  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
	  typedef ROOT::Math::SVector<double, 3> SVector3;
	  typedef ROOT::Math::SVector<double, 6> SVector6;
	  

	  SMatrixSym3D totalCov = vtx.covariance() + sec_covariance;
	  SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
	  
	  dl[it] = ROOT::Math::Mag(distanceVector);
	  dlerror[it] = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl[it];
	  
	  dlos[it] = dl[it]/dlerror[it];
	  

	  //Decay length 2D
	  
	  SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
	  SVector6 v2(sec_covariance(0,0), sec_covariance(0,1),sec_covariance(1,1),0,0,0);

	  SMatrixSym3D sv1(v1);
	  SMatrixSym3D sv2(v2);

	  SMatrixSym3D totalCov2D = sv1 + sv2;
	  SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
	  
	  dl2D[it] = ROOT::Math::Mag(distanceVector2D);
	  double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D[it];
	  
	  dlos2D[it] = dl2D[it]/dl2Derror;
	  
	  //std::cout<<"CHECK-10"<<" for event "<<iEvent.id().event()<<std::endl;
	  //trk info
	  //const pat::PackedCandidate* reco_d1 = dynamic_cast<const pat::PackedCandidate*>(reco_d1);
	  const reco::Track& pseudoTrk1 = reco_d1->pseudoTrack();


	  //trk quality
	  trkChi1[it] = pseudoTrk1.normalizedChi2();
	  ptErr1[it] = pseudoTrk1.ptError();
	  //nhit1[it] = pseudoTrk1.hitPattern().numberOfValidHits();
	  //trkquality1[it] = pseudoTrk1.quality(reco::TrackBase::highPurity);
	  //trkChi1[it] = reco_d1->normalizedChi2();
	  //ptErr1[it] = reco_d1->ptError();
	  
	  
	  //secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
	  	  
	  //DCA
	  math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
	  math::XYZPoint BS_vtx(BSx,BSy,BSz);
	  

	  double dzbest1 = reco_d1->pseudoTrack().dz(bestvtx);	  // today change
	  double dxybest1 = reco_d1->pseudoTrack().dxy(bestvtx); // today change
	  double dzerror1 = std::sqrt(pseudoTrk1.dzError() * pseudoTrk1.dzError() + bestvzError * bestvzError);
	  double dxyerror1 = std::sqrt(pseudoTrk1.dxyError() * pseudoTrk1.dxyError() + bestvxError * bestvyError);
	  
	  dzos1[it] = dzbest1/dzerror1;
	  dxyos1[it] = dxybest1/dxyerror1;
	  Dtrk1Dz1[it] = 1.0*dzbest1;
	  Dtrk1Dxy1[it] = dxybest1;
	  Dtrk1DzError1[it] = dzerror1;
	  Dtrk1DxyError1[it] = dxyerror1;
	 	  

	  //const pat::PackedCandidate* reco_d2 = dynamic_cast<const pat::PackedCandidate*>(reco_d2);
	  const reco::Track& pseudoTrk2 = reco_d2->pseudoTrack();
	  trkChi2[it] = pseudoTrk2.normalizedChi2();
	  ptErr2[it] = pseudoTrk2.ptError();
	  //nhit2[it] = pseudoTrk2.hitPattern().numberOfValidHits();
	  //trkquality2[it] = pseudoTrk2.quality(reco::TrackBase::highPurity);
		  
	  //track Chi2
	  //trkChi2[it] = reco_d2->normalizedChi2();
	  
	  //track pT error
	  //ptErr2[it] = reco_d2->ptError();
	  
	  //vertexCovariance 00-xError 11-y 22-z
	  //secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
	  //std::cout<<"CHECK-11"<<" for event "<<iEvent.id().event()<<std::endl;
		  

	  double dzbest2 = reco_d2->pseudoTrack().dz(bestvtx);	  // today change
	  double dxybest2 = reco_d2->pseudoTrack().dxy(bestvtx); // today change
	  double dzerror2 = TMath::Sqrt(pseudoTrk2.dzError() * pseudoTrk2.dzError() + bestvzError * bestvzError);
	  double dxyerror2 = TMath::Sqrt(pseudoTrk2.dxyError() * pseudoTrk2.dxyError() + bestvxError * bestvyError);

	  dzos2[it] = dzbest2/dzerror2;
	  dxyos2[it] = dxybest2/dxyerror2;
	  Dtrk2Dz1[it] = 1.0*dzbest2;
	  Dtrk2Dxy1[it] = dxybest2;
	  Dtrk2DzError1[it] = dzerror2;
	  Dtrk2DxyError1[it] = dxyerror2;



	  if (threeProngDecay_){
	    //const pat::PackedCandidate* dau3 = dynamic_cast<const pat::PackedCandidate*>(reco_d3);
	    const reco::Track& pseudoTrk3 = reco_d3->pseudoTrack();
	    trkChi3[it] = pseudoTrk3.normalizedChi2();
	    ptErr3[it] = pseudoTrk3.ptError();
	    //nhit2[it] = pseudoTrk2.hitPattern().numberOfValidHits();
	    //trkquality2[it] = pseudoTrk2.quality(reco::TrackBase::highPurity);
	    
	    //track Chi2
	    //trkChi2[it] = reco_d2->normalizedChi2();
	    
	    //track pT error
	    //ptErr2[it] = reco_d2->ptError();
	    
	    //vertexCovariance 00-xError 11-y 22-z
	    //secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
	    
		  

	  double dzbest3 = reco_d3->pseudoTrack().dz(bestvtx);	  // today change
	  double dxybest3 = reco_d3->pseudoTrack().dxy(bestvtx); // today change
	  double dzerror3 = TMath::Sqrt(pseudoTrk3.dzError() * pseudoTrk3.dzError() + bestvzError * bestvzError);
	  double dxyerror3 = TMath::Sqrt(pseudoTrk3.dxyError() * pseudoTrk3.dxyError() + bestvxError * bestvyError);

	  dzos3[it] = dzbest3/dzerror3;
	  dxyos3[it] = dxybest3/dxyerror3;
	  Dtrk3Dz1[it] = 1.0*dzbest3;
	  Dtrk3Dxy1[it] = dxybest3;
	  Dtrk3DzError1[it] = dzerror3;
	  Dtrk3DxyError1[it] = dxyerror3;

	  //std::cout<<"CHECK-12"<<" for event "<<iEvent.id().event()<<std::endl;
	  }

	  
	  
       	  
#ifdef DEBUG
	  cout << "Done reco single iter" << endl;
#endif
	  	  
	}// Candidate loop

	//std::cout<<"Finished Tree producer loop"<<std::endl;
#ifdef DEBUG
	cout << "Fill reco done" << endl;
#endif
}




// ------------ method called once each job just before starting event loop  ------------
void
VertexCompositeTreeProducer::beginJob()
{
	TH1D::SetDefaultSumw2();

	if(!doRecoNtuple_ && !doGenNtuple_)
	{
		cout<<"No output for either RECO or GEN!! Fix config!!"<<endl; return;
	}


	if(saveTree_) initTree();
}

void 
VertexCompositeTreeProducer::initTree()
{ 
	VertexCompositeNtuple = fs->make< TTree>("VertexCompositeNtuple","VertexCompositeNtuple");

	if(doRecoNtuple_) 
	{ 

		// Event info
		VertexCompositeNtuple->Branch("Ntrkoffline",&Ntrkoffline,"Ntrkoffline/I");
		VertexCompositeNtuple->Branch("Npixel",&Npixel,"Npixel/I");
		VertexCompositeNtuple->Branch("HFsumET",&HFsumET,"HFsumET/F");
		VertexCompositeNtuple->Branch("HFsumETPlus",&HFsumETPlus,"HFsumETPlus/F");
		VertexCompositeNtuple->Branch("HFsumETMinus",&HFsumETMinus,"HFsumETMinus/F");
		VertexCompositeNtuple->Branch("ZDCPlus",&ZDCPlus,"ZDCPlus/F");
		VertexCompositeNtuple->Branch("ZDCMinus",&ZDCMinus,"ZDCMinus/F");
		VertexCompositeNtuple->Branch("PvtxX",&bestvx,"PvtxX/F");
		VertexCompositeNtuple->Branch("PvtxY",&bestvy,"PvtxY/F");
		VertexCompositeNtuple->Branch("PvtxZ",&bestvz,"PvtxZ/F");
		VertexCompositeNtuple->Branch("BSx",&BSx,"BSx/F");
		VertexCompositeNtuple->Branch("BSy",&BSy,"BSy/F");
		VertexCompositeNtuple->Branch("BSz",&BSz,"BSz/F");
		VertexCompositeNtuple->Branch("PvtxXErr",&bestvxError,"PvtxXErr/F");
		VertexCompositeNtuple->Branch("PvtxYErr",&bestvyError,"PvtxYErr/F");
		VertexCompositeNtuple->Branch("PvtxZErr",&bestvzError,"PvtxZErr/F");
		VertexCompositeNtuple->Branch("BSxErr",&BSxerror,"BSxErr/F");
		VertexCompositeNtuple->Branch("BSyErr",&BSyerror,"BSyErr/F");
		VertexCompositeNtuple->Branch("BSzErr",&BSzerror,"BSzErr/F");
		VertexCompositeNtuple->Branch("candSize",&candSize,"candSize/I");
		if(isCentrality_) VertexCompositeNtuple->Branch("centrality",&centrality,"centrality/I");
		// particle info
		VertexCompositeNtuple->Branch("pT",&pt,"pT[candSize]/F");
		VertexCompositeNtuple->Branch("y",&y,"y[candSize]/F");
		VertexCompositeNtuple->Branch("phi",&phi,"phi[candSize]/F");
		VertexCompositeNtuple->Branch("mass",&mass,"mass[candSize]/F");
		//if(useAnyMVA_) VertexCompositeNtuple->Branch("mva",&mva,"mva[candSize]/F");

		if(!isSkimMVA_)  
		{
			//Composite candidate info RECO
			VertexCompositeNtuple->Branch("flavor",&flavor,"flavor[candSize]/F");
			VertexCompositeNtuple->Branch("eta",&eta,"eta[candSize]/F");
			VertexCompositeNtuple->Branch("VtxProb",&VtxProb,"VtxProb[candSize]/F");
			VertexCompositeNtuple->Branch("VtxChi2",&vtxChi2,"VtxChi2[candSize]/F");
			VertexCompositeNtuple->Branch("3DPointingAngle",&agl_abs,"3DPointingAngle[candSize]/F");
			VertexCompositeNtuple->Branch("Ddca",&Ddca,"Ddca[candSize]/F");
			VertexCompositeNtuple->Branch("2DPointingAngle",&agl2D_abs,"2DPointingAngle[candSize]/F");
			VertexCompositeNtuple->Branch("3DDecayLengthSignificance",&dlos,"3DDecayLengthSignificance[candSize]/F");
			VertexCompositeNtuple->Branch("3DDecayLength",&dl,"3DDecayLength[candSize]/F");
			VertexCompositeNtuple->Branch("3DDecayLengthError",&dlerror,"3DDecayLengthError[candSize]/F");
			VertexCompositeNtuple->Branch("2DDecayLengthSignificance",&dlos2D,"2DDecayLengthSignificance[candSize]/F");
			VertexCompositeNtuple->Branch("2DDecayLength",&dl2D,"2DDecayLength[candSize]/F");
			VertexCompositeNtuple->Branch("2DDecayLengthError",&dl2Derror,"2DDecayLengthError[candSize]/F");
			VertexCompositeNtuple->Branch("DlxyBS",DlxyBS,"DlxyBS[candSize]/F");
			VertexCompositeNtuple->Branch("DlxyBSErr",DlxyBSErr,"DlxyBSErr[candSize]/F");
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther1",&dzos1,"zDCASignificanceDaugther1[candSize]/F");
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther1",&dxyos1,"xyDCASignificanceDaugther1[candSize]/F");
			VertexCompositeNtuple->Branch("pTD1",&pt1,"pTD1[candSize]/F");
			VertexCompositeNtuple->Branch("pTerrD1",&ptErr1,"pTerrD1[candSize]/F");
			VertexCompositeNtuple->Branch("EtaD1",&eta1,"EtaD1[candSize]/F");
			VertexCompositeNtuple->Branch("PhiD1",&phi1,"PhiD1[candSize]/F");
			VertexCompositeNtuple->Branch("dedxHarmonic2D1",&H2dedx1,"dedxHarmonic2D1[candSize]/F");
			VertexCompositeNtuple->Branch("zDCASignificanceDaugther2",&dzos2,"zDCASignificanceDaugther2[candSize]/F");
			VertexCompositeNtuple->Branch("xyDCASignificanceDaugther2",&dxyos2,"xyDCASignificanceDaugther2[candSize]/F");
			VertexCompositeNtuple->Branch("pTD2",&pt2,"pTD2[candSize]/F");
			VertexCompositeNtuple->Branch("pTerrD2",&ptErr2,"pTerrD2[candSize]/F");
			VertexCompositeNtuple->Branch("EtaD2",&eta2,"EtaD2[candSize]/F");
			VertexCompositeNtuple->Branch("PhiD2",&phi2,"PhiD2[candSize]/F");
			VertexCompositeNtuple->Branch("dedxHarmonic2D2",&H2dedx2,"dedxHarmonic2D2[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk1Dz1",&Dtrk1Dz1,"Dtrk1Dz1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk2Dz1",&Dtrk2Dz1,"Dtrk2Dz1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk1Dxy1",&Dtrk1Dxy1,"Dtrk1Dxy1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk2Dxy1",&Dtrk2Dxy1,"Dtrk2Dxy1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk1DzError1",&Dtrk1DzError1,"Dtrk1DzError1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk2DzError1",&Dtrk2DzError1,"Dtrk2DzError1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk1DxyError1",&Dtrk1DxyError1,"Dtrk1DxyError1[candSize]/F");
			VertexCompositeNtuple->Branch("Dtrk2DxyError1",&Dtrk2DxyError1,"Dtrk2DxyError1[candSize]/F");

			if(threeProngDecay_){
			  VertexCompositeNtuple->Branch("pTD3",&pt3,"pTD3[candSize]/F");
			  VertexCompositeNtuple->Branch("pTerrD3",&ptErr3,"pTerrD3[candSize]/F");
			  VertexCompositeNtuple->Branch("EtaD3",&eta3,"EtaD3[candSize]/F");
			  VertexCompositeNtuple->Branch("PhiD3",&phi3,"PhiD3[candSize]/F");
			  //VertexCompositeNtuple->Branch("dedxHarmonic2D3",&H2dedx3,"dedxHarmonic2D3[candSize]/F");
			  VertexCompositeNtuple->Branch("Dtrk3Dz1",&Dtrk3Dz1,"Dtrk3Dz1[candSize]/F");
			  VertexCompositeNtuple->Branch("Dtrk3Dxy1",&Dtrk3Dxy1,"Dtrk3Dxy1[candSize]/F");
			  VertexCompositeNtuple->Branch("Dtrk3DzError1",&Dtrk3DzError1,"Dtrk3DzError1[candSize]/F");
			  VertexCompositeNtuple->Branch("Dtrk3DxyError1",&Dtrk3DxyError1,"Dtrk3DxyError1[candSize]/F");
			  
			}

			if(doGenMatching_)
			{
				VertexCompositeNtuple->Branch("isSwap",&isSwap,"isSwap[candSize]/O");
				VertexCompositeNtuple->Branch("idmom_reco",&idmom_reco,"idmom_reco[candSize]/I");
				VertexCompositeNtuple->Branch("idD1_reco",&idd1_reco,"idD1_reco[candSize]/I");
				VertexCompositeNtuple->Branch("idD2_reco",&idd2_reco,"idD2_reco[candSize]/I");
				VertexCompositeNtuple->Branch("idD3_reco",&idd3_reco,"idD3_reco[candSize]/I");
				VertexCompositeNtuple->Branch("matchGEN",&matchGEN,"matchGEN[candSize]/O");
				VertexCompositeNtuple->Branch("gen3DPointingAngle",&gen_agl_abs,"gen3DPointingAngle[candSize]/F");
				VertexCompositeNtuple->Branch("gen2DPointingAngle",&gen_agl2D_abs,"gen2DPointingAngle[candSize]/F");
				VertexCompositeNtuple->Branch("gen3DDecayLength",&gen_dl,"gen3DDecayLength[candSize]/F");
				VertexCompositeNtuple->Branch("gen2DDecayLength",&gen_dl2D,"gen2DDecayLength[candSize]/F");
			}

		}

	} // doRecoNtuple_

	if(doGenNtuple_)
	{
	  VertexCompositeNtuple->Branch("pT_gen",&pt_gen,"pT_gen[candSize]/F");
	  VertexCompositeNtuple->Branch("eta_gen",&eta_gen,"eta_gen[candSize]/F");
	  VertexCompositeNtuple->Branch("y_gen",&y_gen,"y_gen[candSize]/F");
	  VertexCompositeNtuple->Branch("phi_gen",&phi_gen,"phi_gen[candSize]/F");
	  VertexCompositeNtuple->Branch("MotherID_gen",&idmom,"MotherID_gen[candSize]/I");
	  
	  if(decayInGen_)
	    {
	      
	      VertexCompositeNtuple->Branch("DauID1_gen",&iddau1,"DauID1_gen[candSize]/I");
	      VertexCompositeNtuple->Branch("DauID2_gen",&iddau2,"DauID2_gen[candSize]/I");
	      VertexCompositeNtuple->Branch("DauID3_gen",&iddau3,"DauID3_gen[candSize]/I");
	    }
	}
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
VertexCompositeTreeProducer::endJob() {

}

void
VertexCompositeTreeProducer::genDecayLength(const uint& it, const reco::GenParticle& gCand) {
	gen_dl[it] = -99.; gen_agl_abs[it] = -99.; gen_dl2D[it] = -99.; gen_agl2D_abs[it] = -99.;
	if(gCand.numberOfDaughters()==0 || !gCand.daughter(0)) return;
	const auto& dauVtx = gCand.daughter(0)->vertex();
	TVector3 ptosvec(dauVtx.X(), dauVtx.Y(), dauVtx.Z());
	TVector3 secvec(gCand.px(), gCand.py(), gCand.pz());
	gen_agl_abs[it] = secvec.Angle(ptosvec);
	gen_dl[it]  = ptosvec.Mag();
	TVector3 ptosvec2D(dauVtx.X(), dauVtx.Y(), 0.0);
	TVector3 secvec2D(gCand.px(), gCand.py(), 0.0);
	gen_agl2D_abs[it] = secvec2D.Angle(ptosvec2D);
	gen_dl2D[it]  = ptosvec2D.Mag();
}

//define this as a plug-in
DEFINE_FWK_MODULE(VertexCompositeTreeProducer);
