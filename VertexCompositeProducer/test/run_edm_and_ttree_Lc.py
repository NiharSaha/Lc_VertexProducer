import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Era_Run3_pp_on_PbPb_2023_cff import Run3_pp_on_PbPb_2023

process = cms.Process('ANASKIM', eras.Run3_2023, Run3_pp_on_PbPb_2023)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')
process.load('HeavyIonsAnalysis.EventAnalysis.hievtanalyzer_data_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.hltobject_cfi')




process.load("HeavyIonsAnalysis.EventAnalysis.HiForestInfo_cfi")
process.HiForestInfo.info = cms.vstring("HiForest, miniAOD, 132X, data")


# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000)) 
process.TFileService = cms.Service("TFileService",
    fileName =cms.string('TTree_Lc.root'))


# Define the input source

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    #fileNames = cms.untracked.vstring('file:output_withMC.root'  # Use the EDM output file
    fileNames = cms.untracked.vstring(
       #'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/374/668/00000/06179488-b7e6-44f6-bec9-eb242a290ffd.root'
        'root://xrootd-cms.infn.it//store/hidata/HIRun2023A/HIPhysicsRawPrime0/MINIAOD/PromptReco-v2/000/375/055/00000/2d8cd07d-f92f-44df-8e0f-eb28dca3108b.root'


    ),
        #eventsToProcess = cms.untracked.VEventRange('1:1430:199505260')  # Replace with your specific run, lumi, event numbers
)



#GT and centrality calibration
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v7', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '132X_dataRun3_Prompt_v6', '')
process.HiForestInfo.GlobalTagLabel = process.GlobalTag.globaltag


# Define centrality binning
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")



# =============== Import Sequences =====================

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hffilter_cfi')

#from HeavyIonsAnalysis.EventAnalysis.hltobject_cfi import trigger_list_data_2023_skimmed
#process.hltobject.triggerNames = trigger_list_data_2023_skimmed

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone(
    HLTPaths = [
        "HLT_HIMinimumBias*"
    ]
)
process.eventFilter_HM = cms.Sequence(process.hltFilter)


process.event_filters = cms.Sequence(
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter  *
    process.phfCoincFilter2Th4
)

from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import changeToMiniAOD

# Define the analysis steps

########## LamC3P candidate rereco ###############################################################


VertexCollection_PAT = cms.InputTag("offlineSlimmedPrimaryVertices")
TrackCollection_PAT = cms.InputTag("packedPFCandidates")
GenParticleCollection_PAT = cms.untracked.InputTag("prunedGenParticles")

process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalLamC3PCandidates_cff")

process.generalLamC3PCandidatesNew = process.generalLamC3PCandidates.clone()
process.generalLamC3PCandidatesNew.vertexRecoAlgorithm = VertexCollection_PAT
process.generalLamC3PCandidatesNew.trackRecoAlgorithm = TrackCollection_PAT
process.generalLamC3PCandidatesNew.tkEtaDiffCut = cms.double(999.9)
process.generalLamC3PCandidatesNew.tkNhitsCut = cms.int32(11)
process.generalLamC3PCandidatesNew.tkPtErrCut = cms.double(0.1)
process.generalLamC3PCandidatesNew.tkPtCut = cms.double(2.0)
process.generalLamC3PCandidatesNew.alphaCut = cms.double(0.30)
process.generalLamC3PCandidatesNew.alpha2DCut = cms.double(999.9)
process.generalLamC3PCandidatesNew.dPtCut = cms.double(0.0)
#process.generalLamC3PCandidatesNew.mPiKCutMin = cms.double(1.74)
#process.generalLamC3PCandidatesNew.mPiKCutMax = cms.double(2.00)
#process.generalLamC3PCandidatesNew.d0MassCut = cms.double(0.125)
process.generalLamC3PCandidatesNew.VtxChiProbCut = cms.double(0.010)


# produce LamC3P trees
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3pselector_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.lamc3panalyzer_ntp_cff")
process.load("VertexCompositeAnalysis.VertexCompositeAnalyzer.eventinfotree_cff")



process.LamC3PselectorNewReduced = process.lamc3pselector.clone()
process.LamC3PselectorNewReduced.VertexCollection = VertexCollection_PAT
process.LamC3PselectorNewReduced.TrackCollection = TrackCollection_PAT
process.LamC3PselectorNewReduced.GenParticleCollection = GenParticleCollection_PAT
process.LamC3PselectorNewReduced.DCAValCollection = cms.InputTag("generalLamC3PCandidatesNew:DCAValuesLamC3P")
process.LamC3PselectorNewReduced.DCAErrCollection = cms.InputTag("generalLamC3PCandidatesNew:DCAErrorsLamC3P")
process.LamC3PselectorNewReduced.cand3DDecayLengthSigMin = cms.untracked.double(0.)
process.LamC3PselectorNewReduced.cand3DPointingAngleMax = cms.untracked.double(1.0)
#process.LamC3PselectorNewReduced.trkNHitMin = cms.untracked.int32(11)
process.LamC3PselectorNewReduced.trkNHitMin = cms.untracked.int32(-1)


process.LamC3PAna = process.lamc3pana.clone()
process.LamC3PAna.VertexCollection = VertexCollection_PAT
process.LamC3PAna.TrackCollection = TrackCollection_PAT
process.LamC3PAna.GenParticleCollection = GenParticleCollection_PAT
process.LamC3PAna.VertexCompositeCollection = cms.untracked.InputTag("LamC3PselectorNewReduced:LamC3P")
process.LamC3PAna.DCAValCollection = cms.InputTag("LamC3PselectorNewReduced:DCAValuesNewLamC3P")
process.LamC3PAna.DCAErrCollection = cms.InputTag("LamC3PselectorNewReduced:DCAErrorsNewLamC3P")
process.LamC3PAna.isCentrality = cms.bool(True) # Centrality
process.LamC3PAna.centralityBinLabel = cms.InputTag("centralityBin", "HFtowers")#centrality
process.LamC3PAna.centralitySrc = cms.InputTag("hiCentrality") #central
process.LamC3PAna.doGenNtuple = cms.untracked.bool(False) #MConly
process.LamC3PAna.doGenMatching = cms.untracked.bool(False) #MConly

process.LamC3Pana_seq2 = cms.Sequence(process.LamC3PselectorNewReduced * process.LamC3PAna)


# Define the process schedule
process.EventSelections = cms.Path(process.centralityBin * process.eventFilter_HM * process.event_filters * process.generalLamC3PCandidatesNew * process.LamC3Pana_seq2)
process.schedule = cms.Schedule(process.EventSelections)


changeToMiniAOD(process)
process.options.numberOfThreads = 1

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('edm.root'),
        outputCommands = cms.untracked.vstring( #which data to include and exclude 
        "drop *", #no data is kept unless explicitly specified
        'keep recoVertexCompositeCandidates_LamC3PselectorNewReduced_LamC3P_*',  # Keep the first collection
        'keep recoVertexCompositeCandidates_generalLamC3PCandidatesNew_LamC3P_*' 
        )
)

process.outputPath = cms.EndPath(process.output)
process.schedule.append(process.outputPath)



