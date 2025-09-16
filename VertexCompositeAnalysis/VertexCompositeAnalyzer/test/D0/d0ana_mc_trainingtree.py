import FWCore.ParameterSet.Config as cms

process = cms.Process("d0ana")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options   = cms.untracked.PSet( wantSummary = 
cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) 
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_V0Cascade_v1/170301_201909/0000/pPb_HM_100.root',
'root://cms-xrd-global.cern.ch//store/user/davidlw/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/RecoSkim2016_pPb_D0Both_v2/170922_191620/0000/pPb_HM_1.root',
#'root://cms-xrd-global.cern.ch//store/user/davidlw/PAHighMultiplicity1/RecoSkim2016_pPb_D0_v2/170323_023918/0000/pPb_HM_121.root'
                ),
secondaryFileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/F8B77681-2DB0-E611-926A-FA163EB6A60D.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/ECFD5388-2DB0-E611-845F-FA163E6DC769.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/98D5D280-2DB0-E611-832E-02163E01439F.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/7261827F-2DB0-E611-850D-FA163E8AFBCE.root',
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/517/00000/62856385-2DB0-E611-B60A-02163E014682.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/206ADD2E-3A9B-E711-9561-0242AC11000B.root',
'root://cms-xrd-global.cern.ch//store/himc/pPb816Summer16DR/PromptD0_D0pT-1p2_pPb-EmbEPOS_8p16_Pythia8/AODSIM/pPbEmb_80X_mcRun2_pA_v4-v1/10000/28D799E3-5D9B-E711-BEE3-0242AC110008.root'
#'root://cms-xrd-global.cern.ch//store/hidata/PARun2016C/PAHighMultiplicity1/AOD/PromptReco-v1/000/285/480/00000/0C73EB2B-22AF-E611-85F1-02163E013657.root'
)
                            )

process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0selector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.d0analyzer_ntp_cff")

process.TFileService = cms.Service("TFileService",
                                       fileName = 
cms.string('d0ana_mc.root')
                                   )

process.d0ana_mc_genmatch = process.d0ana_mc.clone()
process.d0ana_mc_genunmatch = process.d0ana_mc.clone()
process.d0ana_mc_genmatchswap = process.d0ana_mc.clone()
process.d0ana_mc_genmatchunswap = process.d0ana_mc.clone()

process.d0ana_mc_genmatch.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatch:D0")
process.d0ana_mc_genunmatch.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenUnMatch:D0")
process.d0ana_mc_genmatchswap.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatchSwap:D0")
process.d0ana_mc_genmatchunswap.VertexCompositeCollection = cms.untracked.InputTag("d0selectorMCGenMatchUnSwap:D0")
process.d0ana_wrongsign_mc.VertexCompositeCollection = cms.untracked.InputTag("d0selectorWSMC:D0")

process.d0ana_genmatch_seq = cms.Sequence(process.d0selectorMCGenMatch * process.d0ana_mc_genmatch)
process.d0ana_genunmatch_seq = cms.Sequence(process.d0selectorMCGenUnMatch * process.d0ana_mc_genunmatch)
process.d0ana_genmatchswap_seq = cms.Sequence(process.d0selectorMCGenMatchSwap * process.d0ana_mc_genmatchswap)
process.d0ana_genmatchunswap_seq = cms.Sequence(process.d0selectorMCGenMatchUnSwap * process.d0ana_mc_genmatchunswap)
#process.d0ana_wrongsign_seq = cms.Sequence(process.d0selectorWSMC * process.d0ana_wrongsign)

#process.p = cms.Path(process.d0ana_genmatch_seq)
process.p1 = cms.Path(process.d0ana_genunmatch_seq)
process.p2 = cms.Path(process.d0ana_genmatchswap_seq)
process.p3 = cms.Path(process.d0ana_genmatchunswap_seq)
#process.p4 = cms.Path(process.d0ana_wrongsign_seq)
