if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'VertexCompositeAna'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
#    config.JobType.maxMemoryMB = 3000
#    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/d0ana_mc_trainingtree_signal.py'
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 10
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/'
    config.Data.publication = False
    config.Data.useParent = True
    config.Data.inputDBS = 'phys03'
#    config.Site.storageSite = 'T2_CH_CERN'
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T3_US_Rice'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'pPb2016_pPbMC_PromptJPsi_TrainingTree_signal_combined_v2'
    config.Data.outputDatasetTag = 'pPb_TrainingTree_signal_combined_v2'
    config.Data.inputDataset = '/Psi1SToMuMu_pTMu-2p5_pPb-Bst_8p16-Pythia8/davidlw-RecoSkim2016_pPb_JPsi_v1-0a66f66d7871bedbe3caf4f3ff1ba8c1/USER'
    config.JobType.psetName = '../test/jpsiana_mc_trainingtree_signal.py'
    submit(config)

    config.General.requestName = 'pPb2016_PbpMC_PromptJPsi_TrainingTree_signal_combined_v2'
    config.Data.outputDatasetTag = 'Pbp_TrainingTree_signal_combined_v2'
    config.Data.inputDataset = '/Psi1SToMuMu_pTMu-2p5_PbP-Bst_8p16-Pythia8/davidlw-RecoSkim2016_Pbp_JPsi_v1-0a66f66d7871bedbe3caf4f3ff1ba8c1/USER'
    submit(config)
