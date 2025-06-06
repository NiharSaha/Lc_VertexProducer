// -*- C++ -*-
//
// Package:    PrimaryVertexProducer
// Class:      PrimaryVertexProducer
// 
/**\class PrimaryVertexProducer PrimaryVertexProducer.cc RecoVertex/PrimaryVertexProducer/src/PrimaryVertexProducer.cc

 Description: steers tracker primary vertex reconstruction and storage

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Pascal Vanlaer
//         Created:  Tue Feb 28 11:06:34 CET 2006
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFindingBase.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZT_vect.h"


#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/HITrackFilterForPVFinding.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GapClusterizerInZ.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include <algorithm>
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//
// class declaration
//

class PrimaryVertexRecoveryProducer : public edm::stream::EDProducer<> {
public:
  PrimaryVertexRecoveryProducer(const edm::ParameterSet&);
  ~PrimaryVertexRecoveryProducer() override;
  
  void produce(edm::Event&, const edm::EventSetup&) override;
  void skimTracks(reco::VertexCollection&);

  // access to config
  edm::ParameterSet config() const { return theConfig; }
  
private:
  // ----------member data ---------------------------
  TrackFilterForPVFindingBase* theTrackFilter; 
  TrackClusterizerInZ* theTrackClusterizer;

  // vtx fitting algorithms
  struct algo {
    VertexFitter<5> * fitter;
    VertexCompatibleWithBeam * vertexSelector;
    std::string  label;
    bool useBeamConstraint;
    double minNdof;
  };

  std::vector< algo > algorithms;

  edm::ParameterSet theConfig;
  bool fVerbose;

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transTrackBuilderToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  bool redoAllVertices;
  edm::EDGetTokenT<reco::VertexCollection> oldVtxToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimesToken;
  edm::EDGetTokenT<edm::ValueMap<float> > trkTimeResosToken;

  bool f4D;
};
