// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
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



#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

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
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

#include "CondFormats/GBRForest/interface/GBRForest.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416

using namespace std;

class VertexCompositeSelector : public edm::one::EDProducer<> {
public:
  explicit VertexCompositeSelector(const edm::ParameterSet&);
  ~VertexCompositeSelector();

  using MVACollection = std::vector<float>;

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void fillRECO(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;

  double GetMVACut(double y, double pt);

  // ----------member data ---------------------------
    
    //options

    int PID_;
    int PID_dau1_;
    int PID_dau2_;
    
    //cut variables
    double multMax_;
    double multMin_;
    double deltaR_; //deltaR for Gen matching
    bool   trkHighPurity_;
    double trkPMin_;
    double trkPtMin_;
    double trkEtaMax_;
    double trkPSumMin_;
    double trkPtSumMin_;
    double trkPtAsymMin_;
    double trkEtaDiffMax_;
    double trkPtErrMax_;
    int    trkNHitMin_;
    double candpTMin_;
    double candpTMax_;
    double candYMin_;
    double candYMax_;
    double cand3DDecayLengthSigMin_;
    double cand2DDecayLengthSigMin_;
    double cand3DPointingAngleMax_;
    double cand2DPointingAngleMax_;
    double cand3DDCAMin_;
    double cand3DDCAMax_;
    double cand2DDCAMin_;
    double cand2DDCAMax_;
    double candVtxProbMin_;

    //tree branches
    //event info
    int centrality;
    int Ntrkoffline;
    float bestvx;
    float bestvy;
    float bestvz;
    
    //Composite candidate info
    float mva;
    float pt;
    float eta;
    float flavor;
    float y;
    float mass;
    float VtxProb;
    float dlos;
    float dl;
    float dlerror;
    float agl;
    float vtxChi2;
    float ndf;
    float agl_abs;
    float agl2D;
    float agl2D_abs;
    float dlos2D;
    bool isSwap;
    bool matchGEN;
    int idmom_reco;
    


    //dau info
    float dzos1;
    float dzos2;
    float dxyos1;
    float dxyos2;
    float nhit1;
    float nhit2;
    bool trkquality1;
    bool trkquality2;
    float pt1;
    float pt2;
    float ptErr1;
    float ptErr2;
    float p1;
    float p2;
    float eta1;
    float eta2;
    float phi1;
    float phi2;
    float H2dedx1;
    float H2dedx2;
    float T4dedx1;
    float T4dedx2;
    float trkChi1;
    float trkChi2;
    bool  isPionD1;
    bool  isPionD2;
    bool  isKaonD1;
    bool  isKaonD2;
    
    int  selectFlavor_;
    bool usePID_;
    bool useAnyMVA_;
    bool useExistingMVA_;

    std::string mvaType_;
    std::string forestLabel_;
    GBRForest * forest_;
    bool useForestFromDB_;
    std::vector<float> mvaVals_;
    std::string dbFileName_;

    TF2* func_mva;
    std::vector<double> mvaCuts_;

    TH2D* hist_bdtcut;

    float mvaMin_;
    float mvaMax_;

    int   centMin_;
    int   centMax_;
    bool isCentrality_;
  bool isMC_cent_;

    edm::Handle<int> cbin_;

    //tokens
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<MVACollection> MVAValues_Token_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
    edm::EDGetTokenT<int> tok_centBinLabel_;
    edm::EDGetTokenT<reco::Centrality> tok_centSrc_;

    std::string v0IDName_;

    reco::VertexCompositeCandidateCollection theVertexComps;
    MVACollection theMVANew;
    edm::ESGetToken<GBRForest, GBRWrapperRcd> mvaToken_;
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

VertexCompositeSelector::VertexCompositeSelector(const edm::ParameterSet& iConfig)
{
    //options

    PID_ = iConfig.getUntrackedParameter<int>("PID");
    PID_dau1_ = iConfig.getUntrackedParameter<int>("PID_dau1");
    PID_dau2_ = iConfig.getUntrackedParameter<int>("PID_dau2");
    
    //cut variables
    centMin_ = iConfig.getUntrackedParameter<int>("centMin", 0);
    centMax_ = iConfig.getUntrackedParameter<int>("centMax", 10000);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", -1);
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", -1);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.03);
    mvaMax_ = iConfig.getUntrackedParameter<double>("mvaMax", 999.9);
    mvaMin_ = iConfig.getUntrackedParameter<double>("mvaMin", -999.9);

    trkHighPurity_ = iConfig.getUntrackedParameter<bool>("trkHighPurity");
    trkPMin_ = iConfig.getUntrackedParameter<double>("trkPMin");
    trkPtMin_ = iConfig.getUntrackedParameter<double>("trkPtMin");
    trkEtaMax_ = iConfig.getUntrackedParameter<double>("trkEtaMax");
    trkPSumMin_ = iConfig.getUntrackedParameter<double>("trkPSumMin");
    trkPtSumMin_ = iConfig.getUntrackedParameter<double>("trkPtSumMin");
    trkPtAsymMin_ = iConfig.getUntrackedParameter<double>("trkPtAsymMin");
    trkEtaDiffMax_ = iConfig.getUntrackedParameter<double>("trkEtaDiffMax");
    trkPtErrMax_ = iConfig.getUntrackedParameter<double>("trkPtErrMax");
    trkNHitMin_ = iConfig.getUntrackedParameter<int>("trkNHitMin");
    candpTMin_ = iConfig.getUntrackedParameter<double>("candpTMin");
    candpTMax_ = iConfig.getUntrackedParameter<double>("candpTMax");
    candYMin_ = iConfig.getUntrackedParameter<double>("candYMin");
    candYMax_ = iConfig.getUntrackedParameter<double>("candYMax");
    cand3DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand3DDecayLengthSigMin");
    cand2DDecayLengthSigMin_ = iConfig.getUntrackedParameter<double>("cand2DDecayLengthSigMin");
    cand3DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand3DPointingAngleMax");
    cand2DPointingAngleMax_ = iConfig.getUntrackedParameter<double>("cand2DPointingAngleMax");
    cand3DDCAMin_ = iConfig.getUntrackedParameter<double>("cand3DDCAMin");
    cand3DDCAMax_ = iConfig.getUntrackedParameter<double>("cand3DDCAMax");
    cand2DDCAMin_ = iConfig.getUntrackedParameter<double>("cand2DDCAMin");
    cand2DDCAMax_ = iConfig.getUntrackedParameter<double>("cand2DDCAMax");
    candVtxProbMin_ = iConfig.getUntrackedParameter<double>("candVtxProbMin");

    //input tokens
    tok_offlinePV_ = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCollection"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("TrackCollection"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag(iConfig.getUntrackedParameter<edm::InputTag>("GenParticleCollection")));

    usePID_ = false;
    selectFlavor_ = 0;
    if(iConfig.exists("usePID")) usePID_ = iConfig.getParameter<bool>("usePID");
    if(iConfig.exists("useFlavor")) selectFlavor_ = iConfig.getUntrackedParameter<int>("selectFlavor");
 
    // Loading TMVA
    useAnyMVA_ = false;
    useExistingMVA_ = false;

    forestLabel_ = "D0InpPb";
    std::string type = "BDT";
    useForestFromDB_ = true;
    dbFileName_ = "";

    forest_ = nullptr;

    isCentrality_ = false;
    if(iConfig.exists("isCentrality")) isCentrality_ = iConfig.getParameter<bool>("isCentrality");
    if(isCentrality_)
    {
      tok_centBinLabel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("centralityBinLabel"));
      tok_centSrc_ = consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("centralitySrc"));
    }

    if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
    if(iConfig.exists("useExistingMVA")) useExistingMVA_ = iConfig.getParameter<bool>("useExistingMVA");

    if(useAnyMVA_){
      if(useExistingMVA_ && iConfig.exists("MVACollection")){ 
        MVAValues_Token_ = consumes<MVACollection>(iConfig.getParameter<edm::InputTag>("MVACollection"));
      }
      else{
        if(iConfig.exists("mvaType"))type = iConfig.getParameter<std::string>("mvaType");
        if(iConfig.exists("GBRForestLabel"))forestLabel_ = iConfig.getParameter<std::string>("GBRForestLabel");
        if(iConfig.exists("GBRForestFileName")){
          dbFileName_ = iConfig.getParameter<std::string>("GBRForestFileName");
          useForestFromDB_ = false;
        }
      }
   
      mvaCuts_ = iConfig.getParameter< std::vector<double> >("mvaCuts");

      TString bdtcut_filename;
      if(iConfig.exists("BDTCutFileName")) bdtcut_filename = iConfig.getParameter<string>("BDTCutFileName"); 
      hist_bdtcut = 0;
      if(!bdtcut_filename.IsNull()) 
      {
        edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeAnalyzer/data/%s",bdtcut_filename.Data()));
        TFile ff(fip.fullPath().c_str(),"READ");
        hist_bdtcut = (TH2D*)ff.Get("hist_bdtcut");
        ff.Close();
      }

      func_mva = new TF2("func_mva","[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)*(1+[4]*y+[5]*y*y+[6]*y*y*y+[7]*y*y*y*y)",0,5.0,0,100);
      func_mva->SetParameters(mvaCuts_[0],mvaCuts_[1],mvaCuts_[2],mvaCuts_[3],mvaCuts_[4],mvaCuts_[5],mvaCuts_[6],mvaCuts_[7]);
    }

    if(!useForestFromDB_){
      edm::FileInPath fip(Form("VertexCompositeAnalysis/VertexCompositeAnalyzer/data/%s",dbFileName_.c_str()));
      TFile gbrfile(fip.fullPath().c_str(),"READ");
      forest_ = (GBRForest*)gbrfile.Get(forestLabel_.c_str());
      gbrfile.Close();
    }

    mvaType_ = type;

    v0IDName_ = (iConfig.getUntrackedParameter<edm::InputTag>("VertexCompositeCollection")).instance();

    produces< reco::VertexCompositeCandidateCollection >(v0IDName_);
    produces<MVACollection>(Form("MVAValuesNew%s",v0IDName_.c_str()));

    isPionD1 = true;
    isPionD2 = true;
    isKaonD1 = false;
    isKaonD2 = false;

    mvaToken_ = esConsumes<GBRForest, GBRWrapperRcd>(edm::ESInputTag("", forestLabel_));
}


VertexCompositeSelector::~VertexCompositeSelector()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VertexCompositeSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using std::vector;
    using namespace edm;

    fillRECO(iEvent,iSetup);

    auto theNewV0Cands = std::make_unique<reco::VertexCompositeCandidateCollection>();

    theNewV0Cands->reserve( theVertexComps.size() );

    std::copy( theVertexComps.begin(),
               theVertexComps.end(),
               std::back_inserter(*theNewV0Cands) );

    iEvent.put(std::move(theNewV0Cands), v0IDName_);
    theVertexComps.clear();

    if(useAnyMVA_)
    {
      auto mvas = std::make_unique<MVACollection>(theMVANew.begin(),theMVANew.end());
      iEvent.put(std::move(mvas), Form("MVAValuesNew%s",v0IDName_.c_str()));
      theMVANew.clear();
    }
}

void
VertexCompositeSelector::fillRECO(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //get collections
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    const reco::VertexCompositeCandidateCollection * v0candidates_ = v0candidates.product();

    edm::Handle<MVACollection> mvavalues;
    if(useAnyMVA_ && useExistingMVA_)
    {
      iEvent.getByToken(MVAValues_Token_,mvavalues);
      assert( (*mvavalues).size() == v0candidates->size() );
    }

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    if(usePID_) iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    if(usePID_) iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    centrality=-1;
    if(isCentrality_)
    {
      edm::Handle<reco::Centrality> cent;
      iEvent.getByToken(tok_centSrc_, cent);

      iEvent.getByToken(tok_centBinLabel_,cbin_);
      centrality = *cbin_;

      //if(isMC_cent_){centrality = getHiBinFromhiHF(HFsumET);}

    }
    if(centrality!=-1 && (centrality >= centMax_ || centrality < centMin_)) return;

    //best vertex
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

    
    //Ntrkoffline
    Ntrkoffline = 0;
    if(multMax_!=-1 && multMin_!=-1) 
    {
      for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];

        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();


        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        Ntrkoffline++;
      }

      if(Ntrkoffline >= multMax_ || Ntrkoffline < multMin_) return;
    }

    //Gen info for matching

    //RECO Candidate info
    for(unsigned it=0; it<v0candidates_->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_)[it];
        

        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();

        eta = trk.eta();
        y = trk.rapidity();
        pt = trk.pt();
        flavor = trk.pdgId()/421;


        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);

        //Gen match
        // select particle vs antiparticle  
        if(usePID_ && selectFlavor_ && (int)flavor!=selectFlavor_) continue; 

        // select on pT and y
        if(pt<candpTMin_ || pt>candpTMax_) continue;
        if(y<candYMin_ || y>candYMax_) continue;
        //pt
        pt1 = d1->pt();
        pt2 = d2->pt();

        if(pt1 < trkPtMin_ || pt2 < trkPtMin_) continue;
        if((pt1+pt2) < trkPtSumMin_) continue;
                
        if(pt2/pt1 < trkPtAsymMin_ || pt1/pt2 < trkPtAsymMin_) continue;

        //momentum
        p1 = d1->p();
        p2 = d2->p();
        
        if(p1 < trkPMin_ || p2 < trkPMin_) continue;
        if((p1+p2) < trkPSumMin_) continue;

        //eta
        eta1 = d1->eta();
        eta2 = d2->eta();
        
        if(fabs(eta1) > trkEtaMax_ || fabs(eta2) > trkEtaMax_) continue;
        if(fabs(eta1-eta2) > trkEtaDiffMax_) continue;

        //phi
        phi1 = d1->phi();
        phi2 = d2->phi();
        
        
        //vtxChi2
        vtxChi2 = trk.vertexChi2();
        ndf = trk.vertexNdof();
        VtxProb = TMath::Prob(vtxChi2,ndf);

        if(VtxProb < candVtxProbMin_) continue;       


        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
        
        agl = cos(secvec.Angle(ptosvec));
        agl_abs = secvec.Angle(ptosvec);
        if(agl_abs > cand3DPointingAngleMax_) continue;

        agl2D = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs = secvec2D.Angle(ptosvec2D);
        if(agl2D_abs > cand2DPointingAngleMax_) continue;
        
        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl = ROOT::Math::Mag(distanceVector);
        dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        dlos = dl/dlerror;
        if(dlos < cand3DDecayLengthSigMin_ || dlos > 1000.) continue;
 
        //Decay length 2D
        SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
        SVector6 v2(trk.vertexCovariance(0,0), trk.vertexCovariance(0,1),trk.vertexCovariance(1,1),0,0,0);
        
        SMatrixSym3D sv1(v1);
        SMatrixSym3D sv2(v2);
        
        SMatrixSym3D totalCov2D = sv1 + sv2;
        SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);
        
        double dl2D = ROOT::Math::Mag(distanceVector2D);
        double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D;
        
        dlos2D = dl2D/dl2Derror;
        if(dlos2D < cand2DDecayLengthSigMin_ || dlos2D > 1000.) continue;

        double dca3D = dl*sin(agl_abs);
        if(dca3D < cand3DDCAMin_ || dca3D > cand3DDCAMax_) continue;

        double dca2D = dl2D*sin(agl2D_abs);
        if(dca2D < cand2DDCAMin_ || dca2D > cand2DDCAMax_) continue;

        //trk info
        auto dau1 = d1->get<reco::TrackRef>();
            //trk quality
            trkquality1 = dau1->quality(reco::TrackBase::highPurity);
            if(trkHighPurity_ && !trkquality1) continue;
            
            //trk dEdx
            H2dedx1 = -999.9;
            T4dedx1 = -999.9;
            if(usePID_)
            {
               if(dEdxHandle1.isValid())
               {
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
                  H2dedx1 = dEdxTrack[dau1].dEdx();
                  if(H2dedx1 > (2.8/pow(pt1*cosh(eta1),0.4)+0.2) && H2dedx1 < (2.8/pow(pt1*cosh(eta1),0.9)+1.8) && H2dedx1> (2.8/pow(0.75,0.4)+0.2)) { isKaonD1 = true; isPionD1 = false; }
                  if((H2dedx1 < (2.8/pow(pt1*cosh(eta1),0.4)+0.2) || H2dedx1< (2.8/pow(0.75,0.4)+0.2)) && H2dedx1>0) { isPionD1 = true; isKaonD1 = false;}
               }
            
               if(dEdxHandle2.isValid())
               {
                  const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
                  T4dedx1 = dEdxTrack[dau1].dEdx();
               }
            }
            
            //track Chi2
            trkChi1 = dau1->normalizedChi2();
            
            //track pT error
            ptErr1 = dau1->ptError();
             if(ptErr1/dau1->pt() > trkPtErrMax_) continue;
    
            secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
           
            //trkNHits
            nhit1 = dau1->numberOfValidHits();
            if(nhit1 < trkNHitMin_) continue;
            //DCA
            math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
            
            double dzbest1 = dau1->dz(bestvtx);
            double dxybest1 = dau1->dxy(bestvtx);
            double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
            double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
            
            dzos1 = dzbest1/dzerror1;
            dxyos1 = dxybest1/dxyerror1;
        
        auto dau2 = d2->get<reco::TrackRef>();
       
        trkquality2 = dau2->quality(reco::TrackBase::highPurity);
        if(trkHighPurity_ && !trkquality2) continue;

        //trk dEdx
        H2dedx2 = -999.9;
        T4dedx2 = -999.9;
        
        if(usePID_)
        {
          if(dEdxHandle1.isValid())
          {
             const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
             H2dedx2 = dEdxTrack[dau2].dEdx();
        
             if(H2dedx2 > (2.8/pow(pt2*cosh(eta2),0.4)+0.2) && H2dedx2 < (2.8/pow(pt2*cosh(eta2),0.9)+1.8) && H2dedx2> (2.8/pow(0.75,0.4)+0.2)) { isKaonD2 = true; isPionD2 = false; }
             if((H2dedx2 < (2.8/pow(pt2*cosh(eta2),0.4)+0.2) || H2dedx2< (2.8/pow(0.75,0.4)+0.2)) && H2dedx2>0) { isPionD2 = true; isKaonD2 = false; }
          }

          if(dEdxHandle2.isValid()){
             const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
             T4dedx2 = dEdxTrack[dau2].dEdx();
          }

          if(flavor>0 && (!isPionD1 || !isKaonD2)) continue;
          if(flavor<0 && (!isPionD2 || !isKaonD1)) continue;
        }

        //track Chi2
        trkChi2 = dau2->normalizedChi2();
        
        //track pT error
        ptErr2 = dau2->ptError();
        if(ptErr2/dau2->pt() > trkPtErrMax_) continue;

        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        
        //trkNHits
        nhit2 = dau2->numberOfValidHits();
        if(nhit2 < trkNHitMin_) continue;
        
        //DCA
        //math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2 = dzbest2/dzerror2;
        dxyos2 = dxybest2/dxyerror2;
        
        

        // select MVA value
        mva=0;
        if(useAnyMVA_ && useExistingMVA_)
        {
          mva = (*mvavalues)[it];
          if(mva < mvaMin_ || mva > mvaMax_) continue;

          if(mva<GetMVACut(y,pt)) continue;

          theVertexComps.push_back( trk );
          theMVANew.push_back( mva );
          continue;
        }
        else if(useAnyMVA_ && !useExistingMVA_)
        {
          float gbrVals_[50];
          if(forestLabel_ == "D0InpPb" || forestLabel_ == "D0Inpp" || forestLabel_ == "D0InPbPb")
          { 
            gbrVals_[0] = pt;
            gbrVals_[1] = y;
            gbrVals_[2] = VtxProb;
            gbrVals_[3] = dlos;
            gbrVals_[4] = dlos2D;
            gbrVals_[5] = dl;
            gbrVals_[6] = agl_abs;
            gbrVals_[7] = agl2D_abs;
            gbrVals_[8] = dzos1;
            gbrVals_[9] = dzos2;
            gbrVals_[10] = dxyos1;
            gbrVals_[11] = dxyos2;
            gbrVals_[12] = pt1;
            gbrVals_[13] = pt2;
            gbrVals_[14] = eta1;
            gbrVals_[15] = eta2;
            gbrVals_[16] = nhit1;
            gbrVals_[17] = nhit2;
            gbrVals_[18] = ptErr1;
            gbrVals_[19] = ptErr2;
          }

          GBRForest const * forest = forest_;
          if(useForestFromDB_){
            const auto& forestHandle = iSetup.getHandle<GBRForest, GBRWrapperRcd>(mvaToken_);
            forest = forestHandle.product();
          }

          auto gbrVal = forest->GetClassifier(gbrVals_);

          if(gbrVal < mvaMin_ || gbrVal > mvaMax_) continue;
          if(gbrVal < GetMVACut(y,pt)) continue;

          theMVANew.push_back( gbrVal );
        } 
        theVertexComps.push_back( trk );
    }
}

double
VertexCompositeSelector::GetMVACut(double y, double pt)
{
  double mvacut = -1.0;
  if(fabs(y)>2.4) return mvacut;
  if(!hist_bdtcut) return mvacut;

  mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(pt));
  if(pt>7.4) mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(7.4));
  if(pt<1.37) mvacut = hist_bdtcut->GetBinContent(hist_bdtcut->GetXaxis()->FindBin(y),hist_bdtcut->GetYaxis()->FindBin(1.37));

  return mvacut;
}


// ------------ method called once each job just before starting event
//loop  ------------
void
VertexCompositeSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
VertexCompositeSelector::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"
DEFINE_FWK_MODULE(VertexCompositeSelector);
